"""Preliminary study: do PFJet500, AK8PFJet425_SoftDropMass40, PFHT1050 help?

For each (year, MC sample) we measure trigger-OR efficiency in the same
muon-enriched probe as scale_factor_2d.py and compare:

    current_OR     = OR of the HBB triggers already in scale_factor_2d.py
    augmented_OR   = current_OR ∪ (PFJet500 ∪ AK8PFJet425_SoftDropMass40 ∪ PFHT1050)
    gain           = augmented_OR ∧ ¬current_OR     (events ONLY the new triggers catch)

Outputs (per year, per MC sample):
    figures_extra_trig/<year>/<mc>/eff_vs_pt.png
    figures_extra_trig/<year>/<mc>/per_trigger_counts.csv

Usage
-----
    # Branch-availability check (no processing, ~10s/year):
    python3 study_extra_triggers.py --year 2023 --check-only

    # Quick local test (10 files/dataset, ~10 min):
    python3 study_extra_triggers.py --year 2023 --test

    # Full local run (slow; prefer condor):
    python3 study_extra_triggers.py --year 2023 --mc-sample qcd_muenriched
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import awkward as ak
import uproot
import hist
from coffea import processor
from coffea.nanoevents import NanoAODSchema
from coffea.processor import dict_accumulator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

# Reuse the analysis probe selection so SFs measured here are directly
# comparable to scale_factor_2d.py
from scale_factor_2d import (
    create_baseline_selection_mask_new,
    binomial_eff_err,
    sf_err_from_counts,
    trigger_dict_periods as CURRENT_TRIGGERS_BY_YEAR,
    muon_triggers as MUON_REF_BY_YEAR,
    _era_from_key,
    round_robin_sample_lists,
)

# Candidate triggers we are evaluating
CANDIDATE_TRIGGERS = [
    "PFJet500",                       # single PF jet pT > 500
    "AK8PFJet425_SoftDropMass40",     # AK8 pT > 425 + mSD > 40
    "PFHT1050",                       # event HT > 1050
]

# Muon dataset keys per year (mirror scale_factor_2d.py post-fix).
# 2023 is restricted to v3/v4 — v1/v2 NanoAODs lack
# HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06 (see scan_2023 study).
MUON_KEYS_BY_PERIOD = {
    '2022':     ["Muon_Run2022C", "Muon_Run2022D"],
    '2022EE':   ["Muon_Run2022E", "Muon_Run2022F", "Muon_Run2022G"],
    '2023':     [
        "Muon0_Run2023C-v3", "Muon0_Run2023C-v4",
        "Muon1_Run2023C-v3", "Muon1_Run2023C-v4",
    ],
    '2023BPix': [
        "Muon0_Run2023D-v1", "Muon0_Run2023D-v2",
        "Muon1_Run2023D-v1", "Muon1_Run2023D-v2",
    ],
}

PT_EDGES = [200, 280, 360, 450, 550, 700, 1000]  # match scale_factor_2d.py


# ---------------------------------------------------------------------------
# Branch-availability check (uproot only — no coffea/processor needed)
# ---------------------------------------------------------------------------
def check_trigger_branches(year: str, mc_sample: str = "qcd_muenriched") -> dict:
    """Open one data file and one MC file; report which trigger branches exist."""
    print(f"\n=== Branch-availability check for {year} ===")
    current_trigs = CURRENT_TRIGGERS_BY_YEAR.get(year, [])
    all_trigs = current_trigs + CANDIDATE_TRIGGERS

    # one data file
    data_json = os.path.join(f"infiles/{year}", f"{year}_MuonData.json")
    with open(data_json) as f:
        data_dict = json.load(f)
    data_file = None
    for k in MUON_KEYS_BY_PERIOD[year]:
        if data_dict.get(k):
            data_file = data_dict[k][0]
            break
    # one MC file
    if mc_sample == "ttbar":
        mc_json = os.path.join(f"infiles/{year}", f"{year}_ttbar.json")
        with open(mc_json) as f:
            mc_dict = json.load(f)
        mc_file = mc_dict.get("TTtoLNu2Q", [None])[0]
    else:
        mc_json = os.path.join(f"infiles/{year}", f"{year}_QCD_MuEnriched.json")
        with open(mc_json) as f:
            mc_dict = json.load(f)
        # take first non-empty key
        mc_file = next((v[0] for v in mc_dict.values() if v), None)

    results = {}
    for tag, fname in [("data", data_file), (f"mc ({mc_sample})", mc_file)]:
        if fname is None:
            print(f"  {tag}: no file in JSON")
            continue
        try:
            with uproot.open(fname) as f:
                branches = set(f["Events"].keys())
        except Exception as exc:
            print(f"  {tag}: failed to open {fname}: {exc!r}")
            continue
        print(f"  {tag} file: {fname.split('/')[-1]}")
        for t in all_trigs:
            present = f"HLT_{t}" in branches
            mark = "OK " if present else "MISS"
            tag_marker = "(current)  " if t in current_trigs else "(candidate)"
            print(f"    [{mark}] {tag_marker} HLT_{t}")
            results[(tag, t)] = present
    return results


# ---------------------------------------------------------------------------
# Processor
# ---------------------------------------------------------------------------
class ExtraTriggerStudyProcessor(processor.ProcessorABC):
    """Fills baseline + per-trigger + OR counts and a 1D pT histogram per category."""

    def __init__(self, year, current_triggers, candidate_triggers,
                 use_jetid_correction=True):
        self.year = year
        self.current_triggers = list(current_triggers)
        self.candidate_triggers = list(candidate_triggers)
        self.all_triggers = self.current_triggers + self.candidate_triggers
        self.use_jetid_correction = use_jetid_correction

        # Categories we fill: baseline, each trigger individually, current_OR, augmented_OR, gain
        categories = (
            ["baseline"]
            + [f"trig_{t}" for t in self.all_triggers]
            + ["current_OR", "augmented_OR", "gain_new_only"]
        )
        pt_axis = hist.axis.Variable(PT_EDGES, name="pt", label="Leading AK8 jet $p_T$ [GeV]")
        cat_axis = hist.axis.StrCategory(categories, name="cat", label="Category")
        self._accumulator = dict_accumulator({
            "h_pt": hist.Hist(cat_axis, pt_axis, storage="double"),
            "counts": dict_accumulator({c: 0 for c in categories}),
            "missing_triggers": dict_accumulator(),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self.accumulator.copy()
        n0 = len(events)
        if n0 == 0:
            return out

        # ---- muon reference trigger ----
        ref_trigs = MUON_REF_BY_YEAR.get(self.year, [])
        ref_pass = np.zeros(n0, dtype=bool)
        for t in ref_trigs:
            tfield = t[4:] if t.startswith("HLT_") else t
            if tfield in events.HLT.fields:
                ref_pass |= ak.to_numpy(getattr(events.HLT, tfield)).astype(bool)

        # ---- probe (same as scale_factor_2d.py with txbb_region='inclusive') ----
        baseline_mask = create_baseline_selection_mask_new(
            events, tag="Inclusive", year=self.year,
            txbb_region="inclusive",
            use_jetid_correction=self.use_jetid_correction,
        )
        probe_mask = ref_pass & np.asarray(baseline_mask, dtype=bool)
        if probe_mask.sum() == 0:
            return out
        events_p = events[probe_mask]

        # Candidate-fatjet pT (matches the leading-score candidate logic in the main script;
        # for this 1D study we use the leading-pT AK8 as a simple proxy — sufficient
        # for *relative* eff comparison between current vs augmented OR.)
        # NOTE: ak.pad_none(..., 1, clip=True) already gracefully handles the
        # "no FatJet" case by yielding None, which ak.fill_none turns into -1.
        if "FatJet" in events_p.fields:
            cand_pt = ak.to_numpy(
                ak.fill_none(ak.firsts(ak.pad_none(events_p.FatJet.pt, 1, clip=True)), -1.0)
            )
        else:
            cand_pt = np.full(len(events_p), -1.0)

        # ---- per-trigger pass arrays ----
        trig_pass = {}
        missing = []
        for t in self.all_triggers:
            tfield = t[4:] if t.startswith("HLT_") else t
            if tfield in events_p.HLT.fields:
                trig_pass[t] = ak.to_numpy(getattr(events_p.HLT, tfield)).astype(bool)
            else:
                trig_pass[t] = np.zeros(len(events_p), dtype=bool)
                missing.append(t)
        for t in missing:
            out["missing_triggers"][t] = out["missing_triggers"].get(t, 0) + 1

        if self.current_triggers:
            current_or = np.logical_or.reduce([trig_pass[t] for t in self.current_triggers])
        else:
            current_or = np.zeros(len(events_p), dtype=bool)
        if self.candidate_triggers:
            cand_or = np.logical_or.reduce([trig_pass[t] for t in self.candidate_triggers])
        else:
            cand_or = np.zeros(len(events_p), dtype=bool)
        augmented_or = current_or | cand_or
        gain_only = cand_or & ~current_or

        # ---- fill ----
        out["h_pt"].fill(cat="baseline", pt=cand_pt)
        out["counts"]["baseline"] += len(events_p)
        for t, m in trig_pass.items():
            out["h_pt"].fill(cat=f"trig_{t}", pt=cand_pt[m])
            out["counts"][f"trig_{t}"] += int(m.sum())
        for name, m in [("current_OR", current_or),
                        ("augmented_OR", augmented_or),
                        ("gain_new_only", gain_only)]:
            out["h_pt"].fill(cat=name, pt=cand_pt[m])
            out["counts"][name] += int(m.sum())

        return out

    def postprocess(self, accumulator):
        return accumulator


# ---------------------------------------------------------------------------
# Fileset loading
# ---------------------------------------------------------------------------
def load_fileset(year: str, mc_sample: str, quick_files: int | None,
                 outdir_for_log: str) -> dict:
    fileset = {}
    # data
    data_json = os.path.join(f"infiles/{year}", f"{year}_MuonData.json")
    with open(data_json) as f:
        data_dict = json.load(f)
    data_sublists = []
    for k in MUON_KEYS_BY_PERIOD[year]:
        data_sublists.append(list(data_dict.get(k, [])))
    if quick_files:
        data_files = round_robin_sample_lists(data_sublists, quick_files)
    else:
        data_files = [f for sub in data_sublists for f in sub]
    fileset[f"MuonData_{year}"] = data_files
    print(f"  data files: {len(data_files)}")

    # MC
    if mc_sample == "ttbar":
        mc_json = os.path.join(f"infiles/{year}", f"{year}_ttbar.json")
        with open(mc_json) as f:
            mc_dict = json.load(f)
        mc_sublists = [list(mc_dict.get("TTtoLNu2Q", []))]
    else:
        mc_json = os.path.join(f"infiles/{year}", f"{year}_QCD_MuEnriched.json")
        if not os.path.exists(mc_json):
            raise FileNotFoundError(
                f"{mc_json} missing — run make_qcd_muenriched_infiles.py first."
            )
        with open(mc_json) as f:
            mc_dict = json.load(f)
        mc_sublists = [list(v) for v in mc_dict.values()]
    if quick_files:
        mc_files = round_robin_sample_lists(mc_sublists, quick_files)
    else:
        mc_files = [f for sub in mc_sublists for f in sub]
    fileset[f"{mc_sample}_{year}"] = mc_files
    print(f"  mc   files: {len(mc_files)} ({mc_sample})")
    return fileset


# ---------------------------------------------------------------------------
# Output: tables and plots
# ---------------------------------------------------------------------------
def write_outputs(out_data, out_mc, year, mc_sample, current_triggers,
                  candidate_triggers, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    all_triggers = current_triggers + candidate_triggers

    # ---- CSV ----
    rows = []
    base_d = out_data["counts"]["baseline"]
    base_m = out_mc["counts"]["baseline"]

    def _r(x):
        # Coerce numpy scalars / 0-d arrays to Python float before round() —
        # round(np.ndarray, n) raises TypeError.
        try:
            return round(float(x), 4)
        except Exception:
            return 0.0

    def row(label, k_d, k_m):
        e_d, ue_d = binomial_eff_err(k_d, base_d) if base_d else (0.0, 0.0)
        e_m, ue_m = binomial_eff_err(k_m, base_m) if base_m else (0.0, 0.0)
        if float(e_m) > 0:
            sf = float(e_d) / float(e_m)
            _, sf_err = sf_err_from_counts(k_d, base_d, k_m, base_m)
        else:
            sf, sf_err = 0.0, 0.0
        return {
            "label": label,
            "data_pass": int(k_d), "data_total": int(base_d),
            "data_eff": _r(e_d), "data_eff_err": _r(ue_d),
            "mc_pass":   int(k_m), "mc_total":   int(base_m),
            "mc_eff":    _r(e_m), "mc_eff_err":   _r(ue_m),
            "sf":        _r(sf),  "sf_err":       _r(sf_err),
        }

    for t in all_triggers:
        tag = "current" if t in current_triggers else "candidate"
        rows.append({"group": tag, **row(t, out_data["counts"][f"trig_{t}"],
                                         out_mc["counts"][f"trig_{t}"])})
    for name in ["current_OR", "augmented_OR", "gain_new_only"]:
        rows.append({"group": "summary",
                     **row(name, out_data["counts"][name], out_mc["counts"][name])})

    csv_path = os.path.join(save_dir, "per_trigger_counts.csv")
    headers = ["group", "label", "data_pass", "data_total", "data_eff", "data_eff_err",
               "mc_pass", "mc_total", "mc_eff", "mc_eff_err", "sf", "sf_err"]
    with open(csv_path, "w") as f:
        f.write(",".join(headers) + "\n")
        for r in rows:
            f.write(",".join(str(r[h]) for h in headers) + "\n")
    print(f"  wrote {csv_path}")

    # ---- Headline printout ----
    print("\n--- Headline ---")
    print(f"{'category':24s}  {'data_eff':>14s}  {'mc_eff':>14s}  {'SF':>12s}")
    for r in rows:
        if r["group"] != "summary":
            continue
        de = f"{r['data_eff']*100:6.2f}% +- {r['data_eff_err']*100:5.2f}%"
        me = f"{r['mc_eff']*100:6.2f}% +- {r['mc_eff_err']*100:5.2f}%"
        sf = f"{r['sf']:.3f} +- {r['sf_err']:.3f}"
        print(f"  {r['label']:24s}  {de:>14s}  {me:>14s}  {sf:>12s}")

    cur_d  = next(r for r in rows if r["label"] == "current_OR")
    aug_d  = next(r for r in rows if r["label"] == "augmented_OR")
    gain_d = next(r for r in rows if r["label"] == "gain_new_only")
    print(f"\n  Data eff:   {cur_d['data_eff']*100:.2f}% -> {aug_d['data_eff']*100:.2f}% "
          f"(+{(aug_d['data_eff']-cur_d['data_eff'])*100:.2f}%, unique gain {gain_d['data_eff']*100:.2f}%)")
    print(f"  MC eff:     {cur_d['mc_eff']*100:.2f}% -> {aug_d['mc_eff']*100:.2f}% "
          f"(+{(aug_d['mc_eff']-cur_d['mc_eff'])*100:.2f}%)")
    print(f"  SF:         {cur_d['sf']:.3f} -> {aug_d['sf']:.3f} "
          f"(Δ={aug_d['sf']-cur_d['sf']:+.3f})")

    # ---- Plot: eff vs pT for current_OR, augmented_OR ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, (label, accum, total_k) in zip(
        axes,
        [("Data",      out_data, "data"),
         (f"MC ({mc_sample})", out_mc,   "mc")],
    ):
        base = accum["h_pt"][{"cat": "baseline"}].values()
        edges = accum["h_pt"].axes[1].edges
        centers = 0.5 * (edges[:-1] + edges[1:])

        for cat, color, ls in [
            ("current_OR",   "C0", "-"),
            ("augmented_OR", "C3", "-"),
            ("gain_new_only","C2", "--"),
        ]:
            num = accum["h_pt"][{"cat": cat}].values()
            with np.errstate(invalid="ignore", divide="ignore"):
                eff = np.where(base > 0, num / base, 0.0)
                err = np.where(base > 0, np.sqrt(eff * (1 - eff) / np.maximum(base, 1)), 0.0)
            ax.errorbar(centers, eff, yerr=err, fmt="o-" if ls == "-" else "o--",
                        color=color, label=cat, capsize=2)
        ax.set_xlabel("Leading AK8 jet $p_T$ [GeV]")
        ax.set_title(f"{label}  ({year})")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=10)
    axes[0].set_ylabel("Trigger efficiency")
    fig.tight_layout()
    png_path = os.path.join(save_dir, "eff_vs_pt.png")
    fig.savefig(png_path, dpi=110)
    plt.close(fig)
    print(f"  wrote {png_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--year", required=True,
                   choices=list(CURRENT_TRIGGERS_BY_YEAR.keys()))
    p.add_argument("--mc-sample", default="qcd_muenriched",
                   choices=["ttbar", "qcd_muenriched"])
    p.add_argument("--check-only", action="store_true",
                   help="Only check branch availability and exit.")
    p.add_argument("--test", action="store_true",
                   help="Quick mode: 10 files/dataset, 2 chunks max.")
    p.add_argument("--quick-files", type=int, default=None,
                   help="Override files/dataset (default: 10 with --test, None otherwise).")
    p.add_argument("--chunksize", type=int, default=100_000)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--out-base", default="figures_extra_trig")
    args = p.parse_args()

    if args.check_only:
        check_trigger_branches(args.year, args.mc_sample)
        return

    quick_files = args.quick_files
    maxchunks = None
    if args.test and quick_files is None:
        quick_files = 10
        maxchunks = 2

    current_triggers = CURRENT_TRIGGERS_BY_YEAR[args.year]
    save_dir = os.path.join(args.out_base, args.year, args.mc_sample)

    print(f"=== Extra-trigger study: year={args.year}, mc={args.mc_sample} ===")
    print(f"  current OR: {current_triggers}")
    print(f"  candidates: {CANDIDATE_TRIGGERS}")
    print(f"  save dir:   {save_dir}")
    print(f"  test mode:  {bool(args.test)} (quick_files={quick_files}, maxchunks={maxchunks})")

    fileset = load_fileset(args.year, args.mc_sample, quick_files, save_dir)
    if any(len(v) == 0 for v in fileset.values()):
        print("ERROR: empty fileset for one of the keys; aborting.")
        sys.exit(2)

    proc = ExtraTriggerStudyProcessor(
        year=args.year,
        current_triggers=current_triggers,
        candidate_triggers=CANDIDATE_TRIGGERS,
        use_jetid_correction=True,
    )
    runner = processor.Runner(
        executor=processor.FuturesExecutor(
            compression=None, workers=args.workers, recoverable=True,
        ),
        schema=NanoAODSchema,
        skipbadfiles=True,
        savemetrics=True,
        chunksize=args.chunksize,
        maxchunks=maxchunks,
        xrootdtimeout=300,
    )

    outputs = {}
    for key, files in fileset.items():
        if not files:
            continue
        print(f"\nRunning on {key} ({len(files)} files)...")
        try:
            out = runner({key: files}, treename="Events", processor_instance=proc)
        except Exception as exc:
            print(f"  ERROR for {key}: {exc!r}")
            continue
        if not isinstance(out, tuple) or out[0] is None:
            print(f"  WARNING: empty output for {key} (all files failed metadata?)")
            continue
        outputs[key] = out[0]

    data_key = next((k for k in outputs if k.startswith("MuonData_")), None)
    mc_key   = next((k for k in outputs if not k.startswith("MuonData_")), None)
    if data_key is None or mc_key is None:
        print("ERROR: missing data or MC output; nothing to compare.")
        sys.exit(3)

    write_outputs(
        outputs[data_key], outputs[mc_key],
        args.year, args.mc_sample,
        current_triggers, CANDIDATE_TRIGGERS,
        save_dir=save_dir,
    )


if __name__ == "__main__":
    main()
