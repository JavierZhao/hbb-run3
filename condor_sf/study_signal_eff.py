"""Signal-MC trigger efficiency: ggF / VBF (and any other signal sample)
under the current_OR vs augmented_OR trigger soups, with the analysis baseline
(pT > 450, mSD > 40) but WITHOUT the Z(μμ)-style muon veto.

Why this exists: scale_factor_2d.py measures SF on data + QCD/ttbar MC using a
Z(μμ)-orthogonal probe. For *signal* MC (ggF→bb, VBF→bb, …) we just want to
know what fraction of analysis-region events fire the trigger soup — no probe,
no orthogonal trigger. Adds a per-sample efficiency table that's directly
usable as the "trigger acceptance" number in the analysis cutflow.

Output per (year, sample):
    figures_signal_eff/<year>/<sample>/per_trigger_efficiency.csv
    figures_signal_eff/<year>/<sample>/eff_vs_pt.png

Usage
-----
    # Quick check (10 files per sample, ~5 min):
    python3 study_signal_eff.py --year 2022 --sample both --quick-files 10

    # Full run:
    python3 study_signal_eff.py --year 2022 --sample both
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
import hist
from coffea import processor
from coffea.nanoevents import NanoAODSchema
from coffea.processor import dict_accumulator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

# Reuse helpers + trigger lists from scale_factor_2d.py
from scale_factor_2d import (
    create_baseline_selection_mask_new,
    binomial_eff_err,
    trigger_dict_periods as CURRENT_TRIGGERS_BY_YEAR,
    round_robin_sample_lists,
)
from study_extra_triggers import CANDIDATE_TRIGGERS, PT_EDGES


# JSON keys for signal MC per year. Conventionally the JSON has 1 key whose
# name we don't try to enforce — we just take all files across all keys.
SIGNAL_KEYS_BY_SAMPLE = ("ggF", "VBF", "WH", "ZH", "ttH")


# ---------------------------------------------------------------------------
# Processor
# ---------------------------------------------------------------------------
class SignalEffProcessor(processor.ProcessorABC):
    """Counts events that pass the analysis baseline (pT>450, mSD>40, no muon
    veto) and then check which trigger(s) fire. Fills counts + 1D pT histogram
    per category (baseline, each trigger, current_OR, augmented_OR, gain_only).
    """

    def __init__(self, year, current_triggers, candidate_triggers,
                 use_jetid_correction=True,
                 probe_pt_min=None, probe_ht_min=None):
        self.year = year
        self.current_triggers = list(current_triggers)
        self.candidate_triggers = list(candidate_triggers)
        self.all_triggers = self.current_triggers + self.candidate_triggers
        self.use_jetid_correction = use_jetid_correction
        self.probe_pt_min = probe_pt_min
        self.probe_ht_min = probe_ht_min

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
        if len(events) == 0:
            return out

        # ---- baseline WITHOUT muon veto, WITHOUT reference trigger ----
        baseline_mask = create_baseline_selection_mask_new(
            events, tag="Inclusive", year=self.year,
            txbb_region="inclusive",
            use_jetid_correction=self.use_jetid_correction,
            require_muon=False,
        )
        baseline_mask = np.asarray(baseline_mask, dtype=bool)
        if baseline_mask.sum() == 0:
            return out
        events_p = events[baseline_mask]

        # ---- optional further-tightening probe cuts (rarely needed for signal) ----
        if self.probe_pt_min is not None or self.probe_ht_min is not None:
            tight_mask = np.ones(len(events_p), dtype=bool)
            if self.probe_pt_min is not None:
                cand_pt_pre = ak.fill_none(
                    ak.firsts(ak.pad_none(events_p.FatJet.pt, 1, clip=True)), -1.0
                )
                tight_mask &= ak.to_numpy(cand_pt_pre) > float(self.probe_pt_min)
            if self.probe_ht_min is not None:
                ak4_for_ht = (events_p.Jet.pt > 30) & (abs(events_p.Jet.eta) < 2.5)
                ht_pre = ak.sum(events_p.Jet.pt[ak4_for_ht], axis=1)
                tight_mask &= ak.to_numpy(ht_pre) > float(self.probe_ht_min)
            events_p = events_p[tight_mask]
            if len(events_p) == 0:
                return out

        # Leading AK8 pT (-1 if no AK8 in event)
        if "FatJet" in events_p.fields:
            cand_pt = ak.to_numpy(
                ak.fill_none(ak.firsts(ak.pad_none(events_p.FatJet.pt, 1, clip=True)), -1.0)
            )
        else:
            cand_pt = np.full(len(events_p), -1.0)

        # ---- per-trigger pass arrays ----
        trig_pass = {}
        for t in self.all_triggers:
            tfield = t[4:] if t.startswith("HLT_") else t
            if tfield in events_p.HLT.fields:
                trig_pass[t] = ak.to_numpy(getattr(events_p.HLT, tfield)).astype(bool)
            else:
                trig_pass[t] = np.zeros(len(events_p), dtype=bool)
                out["missing_triggers"][t] = out["missing_triggers"].get(t, 0) + 1

        current_or = np.logical_or.reduce([trig_pass[t] for t in self.current_triggers]) \
            if self.current_triggers else np.zeros(len(events_p), dtype=bool)
        cand_or = np.logical_or.reduce([trig_pass[t] for t in self.candidate_triggers]) \
            if self.candidate_triggers else np.zeros(len(events_p), dtype=bool)
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
def load_signal_fileset(year: str, sample: str, quick_files):
    json_path = os.path.join(f"infiles/{year}", f"{year}_{sample}.json")
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"{json_path} missing — check signal-sample name.")
    with open(json_path) as f:
        sample_dict = json.load(f)
    # Sample JSONs may have multiple dataset keys (e.g. different mass points).
    # Round-robin so quick-mode sees all of them.
    sublists = [list(v) for v in sample_dict.values() if isinstance(v, list)]
    if quick_files:
        files = round_robin_sample_lists(sublists, quick_files)
    else:
        files = [f for sub in sublists for f in sub]
    return files


# ---------------------------------------------------------------------------
# Output (CSV + plot)
# ---------------------------------------------------------------------------
def write_outputs(out_mc, year, sample, current_triggers, candidate_triggers, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    all_triggers = current_triggers + candidate_triggers

    base = int(out_mc["counts"]["baseline"])
    rows = []

    def _r(x):
        try: return round(float(x), 4)
        except Exception: return 0.0

    def row(label, k):
        e, ue = binomial_eff_err(k, base) if base else (0.0, 0.0)
        return {
            "label": label,
            "pass":  int(k),
            "total": int(base),
            "eff":   _r(e),
            "eff_err": _r(ue),
        }

    for t in all_triggers:
        tag = "current" if t in current_triggers else "candidate"
        rows.append({"group": tag, **row(t, out_mc["counts"][f"trig_{t}"])})
    for name in ["current_OR", "augmented_OR", "gain_new_only"]:
        rows.append({"group": "summary", **row(name, out_mc["counts"][name])})

    csv_path = os.path.join(save_dir, "per_trigger_efficiency.csv")
    headers = ["group", "label", "pass", "total", "eff", "eff_err"]
    with open(csv_path, "w") as f:
        f.write(",".join(headers) + "\n")
        for r in rows:
            f.write(",".join(str(r[h]) for h in headers) + "\n")
    print(f"  wrote {csv_path}")

    # Headline
    print("\n--- Headline ---")
    print(f"{'category':24s}  {'pass/total':>16s}  {'eff':>14s}")
    for r in rows:
        if r["group"] == "summary":
            print(f"  {r['label']:24s}  {r['pass']:>7d}/{r['total']:<8d}  {r['eff']*100:6.2f}% +- {r['eff_err']*100:5.2f}%")

    cur = next(r for r in rows if r["label"] == "current_OR")
    aug = next(r for r in rows if r["label"] == "augmented_OR")
    print(f"\n  {sample} {year}: current_OR -> augmented_OR  "
          f"{cur['eff']*100:.2f}% -> {aug['eff']*100:.2f}%  "
          f"(+{(aug['eff']-cur['eff'])*100:.2f}%)")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    base_vals = out_mc["h_pt"][{"cat": "baseline"}].values()
    edges = out_mc["h_pt"].axes[1].edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    for cat, color, ls in [
        ("current_OR",    "C0", "-"),
        ("augmented_OR",  "C3", "-"),
        ("gain_new_only", "C2", "--"),
    ]:
        num = out_mc["h_pt"][{"cat": cat}].values()
        with np.errstate(invalid="ignore", divide="ignore"):
            eff = np.where(base_vals > 0, num / base_vals, 0.0)
            err = np.where(base_vals > 0,
                           np.sqrt(eff * (1 - eff) / np.maximum(base_vals, 1)),
                           0.0)
        ax.errorbar(centers, eff, yerr=err, fmt="o-" if ls == "-" else "o--",
                    color=color, label=cat, capsize=2)
    ax.set_xlabel("Leading AK8 jet $p_T$ [GeV]")
    ax.set_ylabel("Trigger efficiency")
    ax.set_title(f"{sample} signal MC ({year})  ·  baseline pT > 450, mSD > 40")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=10)
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
    p.add_argument("--sample", default="both",
                   choices=["ggF", "VBF", "WH", "ZH", "ttH", "both"],
                   help="Signal sample. 'both' = ggF + VBF in one job.")
    p.add_argument("--quick-files", type=int, default=None,
                   help="Limit to N files/sample (round-robin across JSON keys).")
    p.add_argument("--chunksize", type=int, default=100_000)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--out-base", default="figures_signal_eff")
    p.add_argument("--probe-pt-min", type=float, default=None)
    p.add_argument("--probe-ht-min", type=float, default=None)
    args = p.parse_args()

    samples = ["ggF", "VBF"] if args.sample == "both" else [args.sample]
    current_triggers = CURRENT_TRIGGERS_BY_YEAR[args.year]

    print(f"=== Signal-MC trigger eff: year={args.year}, samples={samples} ===")
    print(f"  current OR: {current_triggers}")
    print(f"  candidates: {CANDIDATE_TRIGGERS}")
    print(f"  baseline:   pT > 450, mSD > 40 (no muon veto, no ref trigger)")
    if args.probe_pt_min or args.probe_ht_min:
        print(f"  extra cuts: pT > {args.probe_pt_min}, HT > {args.probe_ht_min}")

    runner = processor.Runner(
        executor=processor.FuturesExecutor(
            compression=None, workers=args.workers, recoverable=True,
        ),
        schema=NanoAODSchema,
        skipbadfiles=True,
        savemetrics=True,
        chunksize=args.chunksize,
        xrootdtimeout=300,
    )

    for sample in samples:
        save_dir = os.path.join(args.out_base, args.year, sample)
        files = load_signal_fileset(args.year, sample, args.quick_files)
        print(f"\n--- {sample}: {len(files)} files ---")
        if not files:
            print(f"  (no files; skipping)")
            continue

        proc = SignalEffProcessor(
            year=args.year,
            current_triggers=current_triggers,
            candidate_triggers=CANDIDATE_TRIGGERS,
            use_jetid_correction=True,
            probe_pt_min=args.probe_pt_min,
            probe_ht_min=args.probe_ht_min,
        )

        try:
            out = runner({f"{sample}_{args.year}": files}, treename="Events",
                         processor_instance=proc)
        except Exception as exc:
            print(f"  ERROR running on {sample}: {exc!r}")
            continue
        if not isinstance(out, tuple) or out[0] is None:
            print(f"  WARNING: empty output for {sample} (all files failed metadata?)")
            continue

        write_outputs(out[0], args.year, sample, current_triggers,
                      CANDIDATE_TRIGGERS, save_dir=save_dir)


if __name__ == "__main__":
    main()
