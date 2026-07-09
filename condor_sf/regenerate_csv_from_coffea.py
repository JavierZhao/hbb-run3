"""Regenerate per-trigger SF CSV from already-saved coffea files.

Used after merge_sf_2d_qcd_ptbin.py to derive the integrated Combined-OR SF
from the merged data+MC coffea outputs, without rerunning the heavy processor.

Usage:
    python3 regenerate_csv_from_coffea.py --year 2024 --region fail --mc-sample qcd_muenriched
    python3 regenerate_csv_from_coffea.py --year 2024 --region all --mc-sample qcd_muenriched
"""
from __future__ import annotations
import argparse
import os
import sys

# scale_factor_2d.py is in the same dir; import the helpers we need.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import scale_factor_2d as sf  # noqa: E402
from coffea import util  # noqa: E402


def regen(year, region, mc_sample, out_dir):
    data_path = os.path.join(
        out_dir, f"sf_2d_{year}_Inclusive_{region}_{mc_sample}_data.coffea"
    )
    mc_path = os.path.join(
        out_dir, f"sf_2d_{year}_Inclusive_{region}_{mc_sample}_mc.coffea"
    )
    if not os.path.exists(data_path) or not os.path.exists(mc_path):
        print(f"[{region}] missing {data_path} or {mc_path}; skip")
        return

    print(f"[{region}] loading merged coffea ...")
    output_data = util.load(data_path)
    output_mc = util.load(mc_path)
    print(f"[{region}] data keys: {list(output_data.keys())[:5]}...")
    print(f"[{region}] mc   keys: {list(output_mc.keys())[:5]}...")

    csv_output_dir = os.path.join(
        out_dir, year, "Inclusive", region, mc_sample
    )
    os.makedirs(csv_output_dir, exist_ok=True)

    sf.generate_per_trigger_efficiency_csv(
        output_data, output_mc,
        sf.trig_vars_2d,
        triggers=sf.trigger_dict_periods[year],
        save_dir=csv_output_dir,
        year=f"{year}_Inclusive_{region}_{mc_sample}",
        pt_min=sf.PT_MIN_AN,
        pt_max=sf.PT_MAX_AN,
        msd_min=sf.MSD_MIN_AN,
        msd_max=sf.MSD_MAX_AN,
    )

    # Per-era CSV: needs the per-era dict from scale_factor_2d.py's saved
    # _per_era.coffea (added 2026-07-09). Skipped for older productions
    # that didn't emit it.
    per_era_path = os.path.join(
        out_dir, f"sf_2d_{year}_Inclusive_{region}_{mc_sample}_per_era.coffea"
    )
    if os.path.exists(per_era_path):
        per_era_data = util.load(per_era_path)
        eras = sorted(per_era_data.keys())
        print(f"[{region}] per-era file present; regenerating CSV over eras: {eras}")
        sf.generate_per_era_efficiency_csv(
            per_era_data, output_mc, sf.trig_vars_2d,
            triggers=sf.trigger_dict_periods[year],
            save_dir=csv_output_dir,
            year=f"{year}_Inclusive_{region}_{mc_sample}",
            pt_min=sf.PT_MIN_AN,
            pt_max=sf.PT_MAX_AN,
            msd_min=sf.MSD_MIN_AN,
            msd_max=sf.MSD_MAX_AN,
        )
    else:
        print(f"[{region}] no per-era coffea at {per_era_path}; skip per-era CSV")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", required=True)
    ap.add_argument("--region", default="all",
                    help="fail | pass_bb | pass_cc | all")
    ap.add_argument("--mc-sample", default="qcd_muenriched")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = args.out_dir or os.path.join(here, f"sf_2d_{args.year}", "output")
    regions = ["fail", "pass_bb", "pass_cc"] if args.region == "all" else [args.region]
    for r in regions:
        regen(args.year, r, args.mc_sample, out_dir)


if __name__ == "__main__":
    main()
