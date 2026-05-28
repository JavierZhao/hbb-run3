"""Signal-MC trigger efficiency: ggF / VBF (and any other signal sample)
under the current_OR vs augmented_OR trigger soups, with the analysis baseline
(pT > 450, mSD > 40) but WITHOUT the Z(μμ)-style muon veto.

Why this exists: scale_factor_2d.py measures SF on data + QCD/ttbar MC using a
Z(μμ)-orthogonal probe. For *signal* MC (ggF→bb, VBF→bb, …) we just want to
know what fraction of analysis-region events fire the trigger soup — no probe,
no orthogonal trigger. This gives the "trigger acceptance" numbers + plots that
parallel the SF pipeline for the AN.

Mirrors scale_factor_2d.py coverage:
  - efficiency vs leading-jet pT, mSD, and Xbb+Xcc score
  - split by AN Sec 5.1 regions: inclusive, fail, pass_bb, pass_cc

Output per (year, sample, region):
    figures_signal_eff/<year>/<sample>/<region>/per_trigger_efficiency.csv
    figures_signal_eff/<year>/<sample>/<region>/eff_vs_{pt,msd,score}.png

Usage
-----
    python3 study_signal_eff.py --year 2022 --sample both --quick-files 10
    python3 study_signal_eff.py --year 2022 --sample both
    python3 study_signal_eff.py --year 2023 --sample VBF \
        --candidates VBF_DiPFJet125_45_Mjj720_Detajj3p0
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
    get_txbb_txcc_score,
    corrected_jetid_tight_mask,
)
from study_extra_triggers import CANDIDATE_TRIGGERS, PT_EDGES

# AN Sec 5.1 regions (inclusive + the bb/cc/fail split)
REGIONS = ["inclusive", "fail", "pass_bb", "pass_cc"]

# Variable binning (pT matches scale_factor_2d post-pT>450; mSD matches its
# _msd_edges; score is the Xbb+Xcc tagger axis 0..1).
MSD_EDGES = [0, 30, 60, 90, 120, 160, 220, 300]
VAR_SPECS = {
    "pt":    (hist.axis.Variable(PT_EDGES,  name="val"), "Leading AK8 jet $p_T$ [GeV]"),
    "msd":   (hist.axis.Variable(MSD_EDGES, name="val"), "Leading AK8 jet $m_{SD}$ [GeV]"),
    "score": (hist.axis.Regular(20, 0.0, 1.0, name="val"), "Leading AK8 jet $X_{bb}+X_{cc}$"),
}


# ---------------------------------------------------------------------------
# Processor
# ---------------------------------------------------------------------------
class SignalEffProcessor(processor.ProcessorABC):
    """Trigger efficiency on signal MC vs pT / mSD / score, in 4 regions.

    Baseline = analysis preselection (pT>450, mSD>40, FatJet preselection,
    JetID) WITHOUT the muon veto and WITHOUT a reference trigger. The region
    split (fail/pass_bb/pass_cc) reuses create_baseline_selection_mask_new so
    the definitions are identical to the SF pipeline.
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

        self.cats = (
            ["baseline"]
            + [f"trig_{t}" for t in self.all_triggers]
            + ["current_OR", "augmented_OR", "gain_new_only"]
        )
        region_axis = hist.axis.StrCategory(REGIONS, name="region")
        cat_axis = hist.axis.StrCategory(self.cats, name="cat")
        hists = {}
        for vname, (vaxis, _label) in VAR_SPECS.items():
            hists[f"h_{vname}"] = hist.Hist(region_axis, cat_axis, vaxis, storage="double")
        hists["counts"] = dict_accumulator(
            {f"{r}/{c}": 0 for r in REGIONS for c in self.cats}
        )
        hists["missing_triggers"] = dict_accumulator()
        self._accumulator = dict_accumulator(hists)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self.accumulator.copy()
        n0 = len(events)
        if n0 == 0:
            return out

        # ---- per-event observables (leading-pT FatJet for pT/mSD; leading-score
        #      candidate for the tagger-score axis, mirroring scale_factor_2d) ----
        fj_pt = ak.to_numpy(ak.fill_none(
            ak.firsts(ak.pad_none(events.FatJet.pt, 1, clip=True)), -1.0))
        fj_msd = ak.to_numpy(ak.fill_none(
            ak.firsts(ak.pad_none(events.FatJet.msoftdrop, 1, clip=True)), -1.0))

        fatjets = events.FatJet
        fj_tight = corrected_jetid_tight_mask(
            fatjets, jet_type="AK8", year=self.year,
            use_jetid_correction=self.use_jetid_correction)
        cand = fatjets[(fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fj_tight][:, :2]
        cand_scores = get_txbb_txcc_score(cand)
        cand_lead = ak.firsts(cand[ak.argmax(cand_scores, axis=1, keepdims=True)])
        score = ak.to_numpy(ak.fill_none(get_txbb_txcc_score(cand_lead), -1.0))

        val_by_var = {"pt": fj_pt, "msd": fj_msd, "score": score}

        # ---- per-trigger pass arrays (full events) ----
        trig_pass = {}
        for t in self.all_triggers:
            tfield = t[4:] if t.startswith("HLT_") else t
            if tfield in events.HLT.fields:
                trig_pass[t] = ak.to_numpy(getattr(events.HLT, tfield)).astype(bool)
            else:
                trig_pass[t] = np.zeros(n0, dtype=bool)
                out["missing_triggers"][t] = out["missing_triggers"].get(t, 0) + 1

        cur = (np.logical_or.reduce([trig_pass[t] for t in self.current_triggers])
               if self.current_triggers else np.zeros(n0, dtype=bool))
        cnd = (np.logical_or.reduce([trig_pass[t] for t in self.candidate_triggers])
               if self.candidate_triggers else np.zeros(n0, dtype=bool))
        aug = cur | cnd
        gain = cnd & ~cur
        cat_masks = {"baseline": np.ones(n0, dtype=bool)}
        for t in self.all_triggers:
            cat_masks[f"trig_{t}"] = trig_pass[t]
        cat_masks["current_OR"] = cur
        cat_masks["augmented_OR"] = aug
        cat_masks["gain_new_only"] = gain

        # ---- optional probe tightening (same cut for all regions) ----
        tight = np.ones(n0, dtype=bool)
        if self.probe_pt_min is not None:
            tight &= fj_pt > float(self.probe_pt_min)
        if self.probe_ht_min is not None:
            ak4 = (events.Jet.pt > 30) & (abs(events.Jet.eta) < 2.5)
            ht = ak.to_numpy(ak.sum(events.Jet.pt[ak4], axis=1))
            tight &= ht > float(self.probe_ht_min)

        # ---- per region: baseline mask (no muon veto) then fill ----
        for region in REGIONS:
            rmask = np.asarray(create_baseline_selection_mask_new(
                events, tag="Inclusive", year=self.year,
                txbb_region=region,
                use_jetid_correction=self.use_jetid_correction,
                require_muon=False,
            ), dtype=bool)
            rmask &= tight
            if rmask.sum() == 0:
                continue
            for cat, cmask in cat_masks.items():
                sel = rmask & cmask
                ns = int(sel.sum())
                out["counts"][f"{region}/{cat}"] += ns
                if ns == 0:
                    continue
                for vname in VAR_SPECS:
                    out[f"h_{vname}"].fill(region=region, cat=cat,
                                           val=val_by_var[vname][sel])
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
    sublists = [list(v) for v in sample_dict.values() if isinstance(v, list)]
    if quick_files:
        files = round_robin_sample_lists(sublists, quick_files)
    else:
        files = [f for sub in sublists for f in sub]
    return files


# ---------------------------------------------------------------------------
# Output (CSV + plots), per region
# ---------------------------------------------------------------------------
def _r(x):
    try:
        return round(float(x), 4)
    except Exception:
        return 0.0


def write_outputs(out, year, sample, current_triggers, candidate_triggers, save_base):
    all_triggers = current_triggers + candidate_triggers

    for region in REGIONS:
        save_dir = os.path.join(save_base, region)
        os.makedirs(save_dir, exist_ok=True)
        base = int(out["counts"].get(f"{region}/baseline", 0))

        # ---- CSV ----
        def row(label, k):
            e, ue = binomial_eff_err(k, base) if base else (0.0, 0.0)
            return {"label": label, "pass": int(k), "total": int(base),
                    "eff": _r(e), "eff_err": _r(ue)}

        rows = []
        for t in all_triggers:
            tag = "current" if t in current_triggers else "candidate"
            rows.append({"group": tag, **row(t, out["counts"].get(f"{region}/trig_{t}", 0))})
        for name in ["current_OR", "augmented_OR", "gain_new_only"]:
            rows.append({"group": "summary", **row(name, out["counts"].get(f"{region}/{name}", 0))})

        csv_path = os.path.join(save_dir, "per_trigger_efficiency.csv")
        headers = ["group", "label", "pass", "total", "eff", "eff_err"]
        with open(csv_path, "w") as f:
            f.write(",".join(headers) + "\n")
            for r in rows:
                f.write(",".join(str(r[h]) for h in headers) + "\n")

        # Headline to stdout
        cur = next(r for r in rows if r["label"] == "current_OR")
        aug = next(r for r in rows if r["label"] == "augmented_OR")
        print(f"  [{region:9s}] base={base:6d}  current_OR={cur['eff']*100:6.2f}%  "
              f"augmented_OR={aug['eff']*100:6.2f}%  (+{(aug['eff']-cur['eff'])*100:.2f}%)")

        if base == 0:
            continue  # no plots for an empty region

        # ---- plots: eff vs each variable ----
        for vname, (vaxis, vlabel) in VAR_SPECS.items():
            h = out[f"h_{vname}"]
            edges = h.axes[2].edges
            centers = 0.5 * (edges[:-1] + edges[1:])
            base_vals = h[{"region": region, "cat": "baseline"}].values()
            fig, ax = plt.subplots(figsize=(8, 6))
            for cat, color, ls in [("current_OR", "C0", "-"),
                                    ("augmented_OR", "C3", "-"),
                                    ("gain_new_only", "C2", "--")]:
                num = h[{"region": region, "cat": cat}].values()
                with np.errstate(invalid="ignore", divide="ignore"):
                    eff = np.where(base_vals > 0, num / base_vals, 0.0)
                    err = np.where(base_vals > 0,
                                   np.sqrt(np.clip(eff * (1 - eff), 0, None) / np.maximum(base_vals, 1)),
                                   0.0)
                ax.errorbar(centers, eff, yerr=err,
                            fmt="o-" if ls == "-" else "o--",
                            color=color, label=cat, capsize=2)
            ax.set_xlabel(vlabel)
            ax.set_ylabel("Trigger efficiency")
            ax.set_title(f"{sample} {year} · {region} · baseline pT>450, mSD>40")
            ax.set_ylim(0, 1.05)
            ax.grid(alpha=0.3)
            ax.legend(loc="lower right", fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(save_dir, f"eff_vs_{vname}.png"), dpi=110)
            plt.close(fig)
        print(f"    wrote {save_dir}/{{per_trigger_efficiency.csv, eff_vs_pt.png, eff_vs_msd.png, eff_vs_score.png}}")


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
    p.add_argument("--candidates", default=None,
                   help="Comma-separated trigger names (no HLT_ prefix) to use as the "
                        "candidate-OR overlay, replacing the default from study_extra_triggers.")
    args = p.parse_args()

    samples = ["ggF", "VBF"] if args.sample == "both" else [args.sample]
    current_triggers = CURRENT_TRIGGERS_BY_YEAR[args.year]
    if args.candidates:
        candidate_triggers = [t.strip() for t in args.candidates.split(",") if t.strip()]
    else:
        candidate_triggers = CANDIDATE_TRIGGERS

    print(f"=== Signal-MC trigger eff: year={args.year}, samples={samples} ===")
    print(f"  current OR: {current_triggers}")
    print(f"  candidates: {candidate_triggers}")
    print(f"  regions:    {REGIONS}")
    print(f"  variables:  pT, mSD, Xbb+Xcc score")
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
        save_base = os.path.join(args.out_base, args.year, sample)
        files = load_signal_fileset(args.year, sample, args.quick_files)
        print(f"\n--- {sample}: {len(files)} files ---")
        if not files:
            print("  (no files; skipping)")
            continue

        proc = SignalEffProcessor(
            year=args.year,
            current_triggers=current_triggers,
            candidate_triggers=candidate_triggers,
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
                      candidate_triggers, save_base=save_base)


if __name__ == "__main__":
    main()
