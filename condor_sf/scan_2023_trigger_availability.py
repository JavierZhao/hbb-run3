"""Scan every 2023 Muon data file and report whether it carries
HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06 (and the other 4 triggers in
the study) as a TBranch.

Why this exists: the --check-only run on a single 2023 data file showed
HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06 missing, while MC has it. If a
non-trivial fraction of 2023 data files lacks this branch, the scale_factor_2d
trigger OR is silently degraded for those files (the code at L539 skips
missing branches), biasing data eff down and producing the observed
~30-point data/MC gap in 2023.

Usage:
    python3 scan_2023_trigger_availability.py [--year 2023] \
        [--max-files 0] [--report scan_2023_report.csv]

By default scans every file listed in infiles/2023/2023_MuonData.json.
Outputs:
  - stdout: per-key summary (files have / files miss)
  - CSV with one row per file: file_path,key,trigger_a_present,trigger_b_present,...
"""
from __future__ import annotations
import argparse
import json
import os
import signal
import sys
import time

import uproot

# Force every print to flush so condor's streamed .out updates in real time
# (tee + python normally line-buffers but flush=True is the seatbelt).
print = lambda *a, _p=print, **kw: _p(*a, **{**kw, "flush": True})


class _FileTimeout(Exception):
    pass


def _open_with_timeout(path, seconds):
    """Open `path` with uproot under an alarm-based wall clock.

    uproot.open() can hang indefinitely on a slow xrootd server (we observed
    >100 min on a single file in the previous run). SIGALRM gives us a hard
    cap. Linux-only, but that's fine for LPC.
    """
    def _h(signum, frame):
        raise _FileTimeout(f"hung > {seconds}s")
    prev = signal.signal(signal.SIGALRM, _h)
    signal.alarm(seconds)
    try:
        return uproot.open(path)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prev)

# Same 5 triggers as study_extra_triggers.py, with HLT_ prefix
TRIGGERS_BY_YEAR = {
    "2022": [
        "HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35",
        "HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65",
        "HLT_PFJet500",
        "HLT_AK8PFJet425_SoftDropMass40",
        "HLT_PFHT1050",
    ],
    "2022EE": [
        "HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35",
        "HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65",
        "HLT_PFJet500",
        "HLT_AK8PFJet425_SoftDropMass40",
        "HLT_PFHT1050",
    ],
    "2023": [
        "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06",
        "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70",
        "HLT_PFJet500",
        "HLT_AK8PFJet425_SoftDropMass40",
        "HLT_PFHT1050",
    ],
    "2023BPix": [
        "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06",
        "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70",
        "HLT_PFJet500",
        "HLT_AK8PFJet425_SoftDropMass40",
        "HLT_PFHT1050",
    ],
}


def scan(year: str, max_files: int, report: str) -> None:
    json_path = os.path.join("infiles", year, f"{year}_MuonData.json")
    with open(json_path) as f:
        muon_dict = json.load(f)
    triggers = TRIGGERS_BY_YEAR[year]

    # Flatten files preserving key tag for per-key reporting
    work = []
    for key, files in muon_dict.items():
        for fpath in files:
            work.append((key, fpath))

    if max_files > 0:
        work = work[:max_files]

    print(f"Scanning {len(work)} files across {len(set(k for k,_ in work))} keys for year {year}")
    print(f"Looking for {len(triggers)} trigger branches:")
    for t in triggers:
        print(f"  - {t}")

    # Per-trigger, per-key tallies
    per_key_present = {}   # key -> trigger -> count
    per_key_total   = {}   # key -> count
    missing_files   = []   # list of (file, key, [missing_triggers])
    errors          = []   # list of (file, err)

    rows = []
    headers = ["file", "key"] + [t.replace("HLT_", "") for t in triggers]
    rows.append(headers)

    t0 = time.time()
    PER_FILE_TIMEOUT = 60  # seconds; the prior run hung > 100 min on one file
    for i, (key, fpath) in enumerate(work, 1):
        per_key_total[key] = per_key_total.get(key, 0) + 1
        per_key_present.setdefault(key, {t: 0 for t in triggers})
        try:
            f = _open_with_timeout(fpath, PER_FILE_TIMEOUT)
            try:
                # Just the Events tree key listing; cheaper than a full keys() walk.
                branches = set(f["Events"].keys())
            finally:
                f.close()
        except Exception as exc:
            errors.append((fpath, type(exc).__name__ + ": " + str(exc)))
            rows.append([fpath, key] + ["ERR"] * len(triggers))
            # Print right away so we see slow/stuck patterns live
            print(f"  [{i:4d}/{len(work)}] ERR  key={key} {fpath.split('/')[-1]}: {type(exc).__name__}: {str(exc)[:80]}")
            continue
        present_flags = [t in branches for t in triggers]
        for t, p in zip(triggers, present_flags):
            if p:
                per_key_present[key][t] += 1
        miss = [t for t, p in zip(triggers, present_flags) if not p]
        if miss:
            missing_files.append((fpath, key, miss))
        rows.append([fpath, key] + ["1" if p else "0" for p in present_flags])

        # Tighter cadence (every 10 files) so the .out shows progress live
        if i % 10 == 0 or i == len(work):
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-6)
            print(f"  [{i:4d}/{len(work)}] elapsed={elapsed:.0f}s rate={rate:.2f} file/s")

    # ---- Per-key summary ----
    print("\n=== Per-MuonData-key summary ===")
    print(f"{'key':30s} {'total':>6s}  " + "  ".join(f"{t.replace('HLT_',''):42s}" for t in triggers))
    for key in sorted(per_key_total):
        tot = per_key_total[key]
        if tot == 0:
            continue
        cells = []
        for t in triggers:
            p = per_key_present[key][t]
            frac = 100.0 * p / tot
            cells.append(f"{p:5d}/{tot:<5d} ({frac:5.1f}%)            ")
        print(f"{key:30s} {tot:>6d}  " + "  ".join(c[:44] for c in cells))

    if errors:
        print(f"\n=== {len(errors)} file open errors (showing first 5) ===")
        for fpath, msg in errors[:5]:
            print(f"  {fpath}: {msg}")

    if missing_files:
        print(f"\n=== {len(missing_files)} file(s) MISS at least one trigger (showing first 10) ===")
        for fpath, key, miss in missing_files[:10]:
            short = fpath.split("/")[-1]
            print(f"  [{key}] {short}: missing {[m.replace('HLT_','') for m in miss]}")

    # ---- CSV ----
    with open(report, "w") as f:
        for row in rows:
            f.write(",".join(str(c) for c in row) + "\n")
    print(f"\nWrote per-file report: {report}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--year", default="2023", choices=list(TRIGGERS_BY_YEAR.keys()))
    p.add_argument("--max-files", type=int, default=0,
                   help="If >0, only scan first N files (debugging).")
    p.add_argument("--report", default=None,
                   help="Per-file CSV report (default: scan_<year>_trig_report.csv)")
    args = p.parse_args()

    report = args.report or f"scan_{args.year}_trig_report.csv"
    scan(args.year, args.max_files, report)


if __name__ == "__main__":
    main()
