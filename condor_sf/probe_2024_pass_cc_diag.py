"""Diagnose why 2024 QCD MuEnriched pass_cc came up empty.

For a handful of 2024 QCD MuEnriched MC files (PT_300to470 + higher, the bins
that contribute at jet pT > 450 GeV), apply the SF baseline (candidatejet
pT > 450, mSD > 40) and report the Xbb / Xcc score distributions under
globalParT3:

  - N(baseline)
  - N with Xbb+Xcc >= 0.82                 -> total in pass region
  - of which N with Xbb > Xcc              -> pass_bb count
  - of which N with Xcc > Xbb              -> pass_cc count
  - fraction with Xcc > Xbb                -> answers the "is it physics" question
  - median Xbb and Xcc for baseline events

If that fraction is ~0, the N/A is real physics: globalParT3 essentially
never commits to a cc hypothesis for QCD jets.
If it's >> 0 but the original SF job still found nothing, the bug is upstream
(probably the 2024C-only data-side bug, plus the JSON contamination we fixed).

Also probes a 2023 QCD MuEnriched file for comparison (uses ParticleNet) to
illustrate the difference between the two tagger semantics.
"""
from __future__ import annotations
import json
import os
import sys

import numpy as np
import uproot
import awkward as ak

print = lambda *a, _p=print, **kw: _p(*a, **{**kw, "flush": True})

INFILES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "infiles")

PT_BASELINE = 450.0
MSD_BASELINE = 40.0
TXBB_WP = 0.82

# Probe targets: 2024 high-pT bins (the ones that contribute at pT>450)
TARGETS_2024 = [
    ("2024", "QCD_MuEnriched_PT_300to470", 2),
    ("2024", "QCD_MuEnriched_PT_470to600", 2),
    ("2024", "QCD_MuEnriched_PT_600to800", 2),
]
TARGETS_2023 = [
    ("2023", "QCD_MuEnriched_PT_300to470", 2),
    ("2023", "QCD_MuEnriched_PT_470to600", 2),
]


def fetch_files(year, key):
    path = os.path.join(INFILES, year, f"{year}_QCD_MuEnriched.json")
    with open(path) as f:
        d = json.load(f)
    return d.get(key, [])


def probe_one_file(url, year):
    """Open a single MC file, apply baseline, extract Xbb/Xcc."""
    try:
        f = uproot.open(url, timeout=60)
        events = f["Events"]
    except Exception as exc:
        print(f"  open failed for {url}: {exc}")
        return None

    fj_branches = [
        "FatJet_pt", "FatJet_msoftdrop", "FatJet_eta",
        "FatJet_globalParT3_Xbb", "FatJet_globalParT3_Xcc", "FatJet_globalParT3_QCD",
        "FatJet_particleNet_XbbVsQCD", "FatJet_particleNet_XccVsQCD",
        "FatJet_particleNet_QCD",
    ]
    available = set(events.keys())
    needed = [b for b in fj_branches if b in available]
    arr = events.arrays(needed, library="ak")

    # The candidate jet is jet 0 (leading by pT)
    jets = ak.zip({b.replace("FatJet_", ""): arr[b] for b in needed})

    # Some events have no jets; keep only those with at least one
    has_jet = ak.num(jets) > 0
    jets = jets[has_jet]
    cand = ak.firsts(jets)

    baseline = (cand.pt > PT_BASELINE) & (cand.msoftdrop > MSD_BASELINE)
    cand = cand[baseline]
    n_base = ak.sum(baseline)

    if n_base == 0:
        return {"n_base": 0}

    # Show BOTH tagger semantics for 2024 so we can compare ParT3 vs ParticleNet.
    # The actual SF code prefers ParT3 (raw probabilities) when both are present.
    results = {}
    for tag_label, key_bb, key_cc, key_qcd in [
        ("globalParT3 (raw probabilities, used by SF code)",
         "globalParT3_Xbb", "globalParT3_Xcc", "globalParT3_QCD"),
        ("particleNet (X / (X+QCD) scores, legacy fallback)",
         "particleNet_XbbVsQCD", "particleNet_XccVsQCD", "particleNet_QCD"),
    ]:
        if key_bb not in cand.fields or key_cc not in cand.fields:
            results[tag_label] = None
            continue
        xbb = np.asarray(ak.fill_none(getattr(cand, key_bb), 0.0))
        xcc = np.asarray(ak.fill_none(getattr(cand, key_cc), 0.0))
        sum_score = xbb + xcc
        pass_mask = sum_score >= TXBB_WP
        pass_bb = pass_mask & (xbb > xcc)
        pass_cc = pass_mask & (xcc > xbb)
        results[tag_label] = {
            "n_pass": int(pass_mask.sum()),
            "n_pass_bb": int(pass_bb.sum()),
            "n_pass_cc": int(pass_cc.sum()),
            "frac_pass_cc_in_pass": (pass_cc.sum() / max(pass_mask.sum(), 1)),
            "med_xbb_base": float(np.median(xbb)),
            "med_xcc_base": float(np.median(xcc)),
            "med_xbb_pass": float(np.median(xbb[pass_mask])) if pass_mask.sum() > 0 else None,
            "med_xcc_pass": float(np.median(xcc[pass_mask])) if pass_mask.sum() > 0 else None,
        }
    return {"n_base": int(n_base), "per_tagger": results}


def main():
    print(f"baseline: pt > {PT_BASELINE}, msoftdrop > {MSD_BASELINE}")
    print(f"pass region: Xbb + Xcc >= {TXBB_WP}")
    print(f"pass_cc subset: pass AND Xcc > Xbb")
    print()
    # Accumulate totals per (year, tagger) so we can compare ParT3 vs ParticleNet.
    yearly = {}  # year -> tagger -> tot dict
    for label, targets in [("2024", TARGETS_2024), ("2023", TARGETS_2023)]:
        print("=" * 70)
        print(f"  {label}")
        print("=" * 70)
        yearly[label] = {}
        for year, key, n_files in targets:
            files = fetch_files(year, key)
            if not files:
                print(f"  no files in {year}/{key}; skip")
                continue
            print(f"\n  -- {key} ({n_files} of {len(files)} files) --")
            for url in files[:n_files]:
                r = probe_one_file(url, year)
                if r is None or r["n_base"] == 0:
                    print(f"     (no baseline events in {os.path.basename(url)})")
                    continue
                print(f"     {os.path.basename(url)[:60]}  baseline={r['n_base']}")
                for tag_lbl, td in r["per_tagger"].items():
                    if td is None:
                        print(f"       [{tag_lbl}] BRANCH MISSING")
                        continue
                    print(f"       [{tag_lbl}]")
                    print(f"         pass={td['n_pass']}  pass_bb={td['n_pass_bb']} "
                          f"pass_cc={td['n_pass_cc']}  "
                          f"frac_pass_cc/pass={td['frac_pass_cc_in_pass']:.3%}")
                    print(f"         median Xbb base/pass = {td['med_xbb_base']:.4f} / "
                          f"{td['med_xbb_pass']:.4f}"
                          if td['med_xbb_pass'] is not None else
                          f"         median Xbb base = {td['med_xbb_base']:.4f}")
                    print(f"         median Xcc base/pass = {td['med_xcc_base']:.4f} / "
                          f"{td['med_xcc_pass']:.4f}"
                          if td['med_xcc_pass'] is not None else
                          f"         median Xcc base = {td['med_xcc_base']:.4f}")
                    tot = yearly[label].setdefault(tag_lbl,
                        {"n_base": 0, "n_pass": 0, "n_pass_bb": 0, "n_pass_cc": 0})
                    tot["n_base"] += r["n_base"]
                    tot["n_pass"] += td["n_pass"]
                    tot["n_pass_bb"] += td["n_pass_bb"]
                    tot["n_pass_cc"] += td["n_pass_cc"]
    # Final summary
    print()
    print("=" * 70)
    print("  TOTALS PER (YEAR, TAGGER)")
    print("=" * 70)
    for year_lbl, tagger_totals in yearly.items():
        for tag_lbl, tot in tagger_totals.items():
            if tot["n_pass"] == 0:
                print(f"  {year_lbl} | {tag_lbl}: baseline={tot['n_base']}, pass=0")
                continue
            print(f"  {year_lbl} | {tag_lbl}:")
            print(f"    baseline={tot['n_base']}, pass={tot['n_pass']}, "
                  f"pass_bb={tot['n_pass_bb']} ({tot['n_pass_bb']/tot['n_pass']:.1%}), "
                  f"pass_cc={tot['n_pass_cc']} ({tot['n_pass_cc']/tot['n_pass']:.1%})")
    print()
    print("Verdict:")
    print("  If 2024 has pass_cc fraction << 2023, it's the tagger semantics")
    print("  (globalParT3 raw probabilities vs particleNet normalized scores).")
    print("  If 2024 pass_cc fraction is comparable to 2023, the original empty")
    print("  pass_cc job result was a code bug elsewhere (probably the same")
    print("  2024C-only data-side merge bug that bit ttbar).")


if __name__ == "__main__":
    main()
