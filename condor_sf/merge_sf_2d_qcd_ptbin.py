"""Merge per-pT-bin QCD MuEnriched SF sub-job outputs into a single per-region SF.

Each sub-job produced by submit_sf_2d_qcd_ptbin.sh writes
    sf_2d_<YEAR>_Inclusive_<region>_qcd_muenriched_PT_<bin>_data.coffea
    sf_2d_<YEAR>_Inclusive_<region>_qcd_muenriched_PT_<bin>_mc.coffea
plus a per-bin CSV/figure tree under qcd_muenriched_PT_<bin>/.

This script:
  1. Picks one canonical 'data' coffea per region (data side is identical across
     pT-bin sub-jobs — they all run the same MuonData_<YEAR>_<era> keys). We use
     the lowest-pT bin's data file as the canonical 'data' so it survives.
  2. Adds the MC coffea histograms from all 12 pT bins into a single merged MC.
  3. Saves the merged outputs at the canonical (non-suffixed) location:
        sf_2d_<YEAR>_Inclusive_<region>_qcd_muenriched_{data,mc}.coffea
  4. Re-runs the plotting + per-trigger CSV from the merged coffea, by invoking
     scale_factor_2d.py with --skip-processing (TODO if not present) — or, if
     that's not wired up, simply prints the merged Combined-OR SF for hand-fill.

Usage:
    python3 merge_sf_2d_qcd_ptbin.py --year 2024 --region fail
    python3 merge_sf_2d_qcd_ptbin.py --year 2024 --region all
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np  # noqa: E402

from coffea import util  # noqa: E402
import hist  # noqa: F401  needed for histogram unpickling


EMPTY_MC_BYTES = 500  # coffea no-entries skeleton is ~470 B; anything bigger has content


def discover_subjob_files(out_dir, year, region):
    """Find all per-pT-bin coffea pairs under out_dir for the given region.

    Returns list of (mc_path, data_path, bin_label) tuples, sorted by bin label.
    Skips MC files smaller than EMPTY_MC_BYTES, which are skeleton outputs from
    pT bins whose events all fall below the SF probe baseline (pT > 450 GeV).
    """
    # Pattern produced by submit_sf_2d_qcd_ptbin.sh / --mc-output-suffix
    mc_glob = os.path.join(
        out_dir, f"sf_2d_{year}_Inclusive_{region}_qcd_muenriched_PT_*_mc.coffea"
    )
    pairs = []
    for mc_path in sorted(glob.glob(mc_glob)):
        data_path = mc_path.replace("_mc.coffea", "_data.coffea")
        if not os.path.exists(data_path):
            print(f"  WARN: no data side for {mc_path}; skipping")
            continue
        if os.path.getsize(mc_path) < EMPTY_MC_BYTES:
            print(f"  SKIP empty bin: {os.path.basename(mc_path)} "
                  f"({os.path.getsize(mc_path)} B < {EMPTY_MC_BYTES})")
            continue
        # Extract the bin tag (PT_<...>) from the filename
        base = os.path.basename(mc_path)
        # ...qcd_muenriched_PT_300to470_mc.coffea -> PT_300to470
        tag = base.split("_qcd_muenriched_")[1].rsplit("_mc.coffea", 1)[0]
        pairs.append((mc_path, data_path, tag))
    return pairs


def add_outputs(outputs):
    """Type-aware accumulator across coffea processor outputs.

    The output dict is a mix of:
      - boost-Hist histograms with fixed binning -> elementwise `+`
      - bare numpy arrays of per-event records (variable size per sub-job)
        -> np.concatenate
      - nested dicts -> recurse
      - scalars (ints, floats) -> `+`

    Coffea's official processor.accumulate fails here because the per-event
    arrays are bare numpy.ndarray, not column_accumulator-wrapped. My earlier
    naive `a + b` silently dropped shape-mismatched leaves, which biased
    per-region MC denominators (pass_bb showed MC baseline=2 vs ~thousands
    of truth, breaking the per-trigger CSV).
    """
    if not outputs:
        return None
    merged = outputs[0]
    for o in outputs[1:]:
        merged = _accumulate(merged, o)
    return merged


def _accumulate(a, b):
    """Recursive type-aware accumulator. See add_outputs for the leaf rules."""
    if isinstance(a, dict):
        out = dict(a)
        for k, v in b.items():
            if k in out:
                out[k] = _accumulate(out[k], v)
            else:
                out[k] = v
        return out
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        # 1D per-event records: concatenate along axis 0.
        # 0-D / scalar arrays: fall back to elementwise + (just below).
        if a.ndim == 1 and b.ndim == 1:
            return np.concatenate([a, b])
    # boost-Hist objects, ints, floats, and equal-shape arrays all support +.
    return a + b


def merge_region(out_dir, year, region):
    pairs = discover_subjob_files(out_dir, year, region)
    if not pairs:
        print(f"[{region}] no per-pT-bin coffea pairs found under {out_dir} — skip")
        return False
    print(f"[{region}] found {len(pairs)} pT-bin sub-jobs:")
    for _, _, t in pairs: print(f"           {t}")

    # The data side is identical across sub-jobs (same MuonData keys). Take the first.
    data = util.load(pairs[0][1])
    print(f"[{region}] using data side from {pairs[0][1]}")

    # Add MC histograms across all pT bins
    mc_objs = [util.load(p[0]) for p in pairs]
    mc_merged = add_outputs(mc_objs)
    print(f"[{region}] merged {len(mc_objs)} MC outputs")

    # Per-era coffea (added 2026-07-09): each sub-job saves a
    # {era: coffea_output_of_that_era_data_only} dict. Data is same across
    # sub-jobs, so we can take the first non-empty per_era file. Only used
    # for downstream per-era CSV regen.
    canonical_pt_tag = pairs[0][2]
    per_era_src = pairs[0][1].replace("_data.coffea", "_per_era.coffea").replace(
        f"_{canonical_pt_tag}_", "_"
    )
    # The above replacement is fragile; be more precise:
    first_mc = pairs[0][0]
    per_era_src = first_mc.replace("_mc.coffea", "_per_era.coffea").replace(
        f"_{canonical_pt_tag}_per_era", f"_{canonical_pt_tag}_per_era"
    )
    # Simpler: use pattern
    import glob as _glob
    per_era_candidates = sorted(_glob.glob(
        os.path.join(out_dir, f"sf_2d_{year}_Inclusive_{region}_qcd_muenriched_PT_*_per_era.coffea")
    ))
    per_era_data = None
    if per_era_candidates:
        per_era_data = util.load(per_era_candidates[0])
        print(f"[{region}] using per-era side from {per_era_candidates[0]}")

    out_data = os.path.join(out_dir, f"sf_2d_{year}_Inclusive_{region}_qcd_muenriched_data.coffea")
    out_mc = os.path.join(out_dir, f"sf_2d_{year}_Inclusive_{region}_qcd_muenriched_mc.coffea")
    out_per_era = os.path.join(out_dir, f"sf_2d_{year}_Inclusive_{region}_qcd_muenriched_per_era.coffea")
    if per_era_data is not None:
        util.save(per_era_data, out_per_era)
        print(f"[{region}] saved merged per-era to {out_per_era}")
    util.save(data, out_data)
    util.save(mc_merged, out_mc)
    print(f"[{region}] saved merged outputs:\n  {out_data}\n  {out_mc}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", required=True)
    ap.add_argument("--region", default="all",
                    help="fail | pass_bb | pass_cc | all")
    ap.add_argument("--out-dir",
                    default=None,
                    help="Path to sf_2d_<YEAR>/output (default inferred from --year)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = args.out_dir or os.path.join(here, f"sf_2d_{args.year}", "output")
    if not os.path.isdir(out_dir):
        sys.exit(f"ERROR: out dir not found: {out_dir}")

    regions = ["fail", "pass_bb", "pass_cc"] if args.region == "all" else [args.region]
    any_done = False
    for r in regions:
        if merge_region(out_dir, args.year, r):
            any_done = True

    if not any_done:
        sys.exit("no regions merged")

    print("\nMerge complete. To regenerate per-trigger CSVs and plots from the merged coffea,")
    print("invoke scale_factor_2d.py's plotting / CSV codepath on the merged files directly,")
    print("or use coffea.util.load + the existing helper functions in scale_factor_2d.py.")


if __name__ == "__main__":
    main()
