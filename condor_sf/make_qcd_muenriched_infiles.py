"""
Generate infiles JSON for QCD muon-enriched MC samples using the DAS REST API.

These are the standard CMS Run 3 QCD_PT-<bin>_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8
NanoAOD samples, used as a more representative MC for the trigger-SF measurement
(jet flavour mix closer to the single-muon data than TTtoLNu2Q).

Usage (requires a valid grid proxy):
    voms-proxy-init --voms cms --valid 168:00
    python3 make_qcd_muenriched_infiles.py

Output (one per period):
    infiles/2022/2022_QCD_MuEnriched.json
    infiles/2022EE/2022EE_QCD_MuEnriched.json
    infiles/2023/2023_QCD_MuEnriched.json
    infiles/2023BPix/2023BPix_QCD_MuEnriched.json
    infiles/2024/2024_QCD_MuEnriched.json

Note: the boosted Hbb analysis baseline is jet pT > 300 GeV, so only the
PT-300to470 and higher bins matter in practice. Lower bins are still queried
so the JSON can be reused for other studies.

2024 caveat: the MC dataset naming convention changed for the
RunIII2024Summer24NanoAODv15 campaign. Old (2022/2023) DAS names follow
    /QCD_PT-<bin>_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8/<campaign>-*/NANOAODSIM
New (2024) names follow
    /QCD_Bin-PT-<bin>_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/<campaign>-*/NANOAODSIM
We branch on the period to pick the right pattern.

Pass `--periods 2024` (or any space-separated list) to limit the run to specific
periods; default is all four (legacy) + 2024.
"""

import argparse
import json
import os
import subprocess
import sys

XROOTD_PREFIX = "root://cmsxrootd.fnal.gov/"

# Standard pT bin labels for the MuEnriched samples.
# The highest bin uses 'PT-1000' (no upper bound).
PT_BINS = [
    "PT-15to20",   "PT-20to30",   "PT-30to50",   "PT-50to80",
    "PT-80to120",  "PT-120to170", "PT-170to300", "PT-300to470",
    "PT-470to600", "PT-600to800", "PT-800to1000", "PT-1000",
]

# Map analysis period -> (campaign string, mc-condition string used in DAS).
# 2024 added on 2026-05-29 — uses new "Bin-PT_..._Fil-MuEnriched" naming convention
# (see DATASET_PATTERN_BY_PERIOD below).
CAMPAIGNS = {
    '2022':     ('Run3Summer22NanoAODv12',         '130X_mcRun3_2022_realistic_v5'),
    '2022EE':   ('Run3Summer22EENanoAODv12',       '130X_mcRun3_2022_realistic_postEE_v6'),
    '2023':     ('Run3Summer23NanoAODv12',         '130X_mcRun3_2023_realistic_v14'),
    '2023BPix': ('Run3Summer23BPixNanoAODv12',     '130X_mcRun3_2023_realistic_postBPix_v2'),
    '2024':     ('RunIII2024Summer24NanoAODv15',   '150X_mcRun3_2024_realistic_v2'),
}

# Two distinct dataset-name conventions are in use across the Run 3 campaigns.
# 2024 (NanoAODv15) re-labeled QCD MC files as `QCD_Bin-PT-..._Fil-MuEnriched`,
# while 2022/2023 (NanoAODv12) still use the legacy `QCD_PT-..._MuEnrichedPt5`.
LEGACY_DATASET = "/QCD_{pt}_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8/{campaign}-*/NANOAODSIM"
V15_DATASET    = "/QCD_Bin-{pt}_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/{campaign}-*/NANOAODSIM"

DATASET_PATTERN_BY_PERIOD = {
    '2022':     LEGACY_DATASET,
    '2022EE':   LEGACY_DATASET,
    '2023':     LEGACY_DATASET,
    '2023BPix': LEGACY_DATASET,
    '2024':     V15_DATASET,
}

# Variant tags to skip when multiple matches exist for a single pt bin.
# We always prefer the standard reco over JEC/BTV variants.
_VARIANT_SKIP_TAGS = ('JMENano', 'BTVNano')


def query_das_dataset(pattern):
    """Query DAS for datasets matching a pattern. Returns list of full dataset names."""
    cmd = ['dasgoclient', '--query', f'dataset={pattern}', '--limit', '0']
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            print(f"  DAS query failed for {pattern}: {r.stderr.strip()}")
            return []
        return [d.strip() for d in r.stdout.strip().split('\n') if d.strip()]
    except FileNotFoundError:
        print("ERROR: dasgoclient not found. Run:")
        print("  source /cvmfs/cms.cern.ch/cmsset_default.sh")
        sys.exit(1)


def query_das_files(dataset):
    """List files for one resolved dataset."""
    cmd = ['dasgoclient', '--query', f'file dataset={dataset}', '--limit', '0']
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        print(f"  file query failed for {dataset}: {r.stderr.strip()}")
        return []
    return [XROOTD_PREFIX + f.strip() for f in r.stdout.strip().split('\n') if f.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Generate QCD MuEnriched infiles JSON via DAS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--periods", nargs="+", choices=list(CAMPAIGNS.keys()),
        default=list(CAMPAIGNS.keys()),
        help="Subset of analysis periods to query (default: all).",
    )
    args = parser.parse_args()

    for period in args.periods:
        campaign, _ = CAMPAIGNS[period]
        ds_template = DATASET_PATTERN_BY_PERIOD[period]
        print(f"\n{'='*60}")
        print(f"Generating QCD_MuEnriched infiles for {period} (campaign {campaign})")
        print(f"  dataset pattern: {ds_template.format(pt='<bin>', campaign=campaign)}")
        print(f"{'='*60}")

        result = {}

        for pt_bin in PT_BINS:
            # Build a wildcarded pattern; DAS resolves to the latest version.
            pattern = ds_template.format(pt=pt_bin, campaign=campaign)
            print(f"  Resolving {pt_bin} ...", flush=True)
            datasets = query_das_dataset(pattern)
            if not datasets:
                print(f"    no dataset matched {pattern}")
                continue

            # Prefer the standard reco over JEC/BTV variants.
            standard = [d for d in datasets if not any(t in d for t in _VARIANT_SKIP_TAGS)]
            datasets = standard if standard else datasets
            # If multiple -vN versions remain, take the highest.
            datasets.sort(key=lambda s: s.rsplit('-v', 1)[-1])
            ds = datasets[-1]
            print(f"    using {ds}")
            files = query_das_files(ds)
            if not files:
                continue

            key = f"QCD_MuEnriched_{pt_bin.replace('-', '_')}"
            result[key] = files
            print(f"    -> {len(files)} files")

        out_path = f"infiles/{period}/{period}_QCD_MuEnriched.json"
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\nSaved {len(result)} pT bins -> {out_path}")
        print(f"Total files: {sum(len(v) for v in result.values())}")


if __name__ == '__main__':
    main()
