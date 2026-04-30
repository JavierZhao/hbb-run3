"""
Generate infiles JSON for ParkingVBF NanoAOD datasets using the DAS REST API.

Usage (requires a valid grid proxy):
    voms-proxy-init --voms cms --valid 168:00
    python3 make_parking_vbf_infiles.py

Output:
    infiles/2023/2023_ParkingVBFData.json
    infiles/2023BPix/2023BPix_ParkingVBFData.json
"""

import json
import os
import subprocess
import sys

XROOTD_PREFIX = "root://cmsxrootd.fnal.gov/"

# ParkingVBF streams: 0-9 (CMS uses 10 parallel parking streams)
STREAMS = [f"ParkingVBF{i}" for i in range(10)]

# Dataset patterns to query
DATASETS = {
    '2023': {
        'datasets': [
            f"/{stream}/Run2023C-22Sep2023_v3-v1/NANOAOD"
            for stream in STREAMS
        ] + [
            f"/{stream}/Run2023C-22Sep2023_v4-v1/NANOAOD"
            for stream in STREAMS
        ],
        'output': 'infiles/2023/2023_ParkingVBFData.json',
        'key_prefix': 'ParkingVBF_Run2023C',
    },
    '2023BPix': {
        'datasets': [
            f"/{stream}/Run2023D-22Sep2023_v1-v1/NANOAOD"
            for stream in STREAMS
        ] + [
            f"/{stream}/Run2023D-22Sep2023_v2-v1/NANOAOD"
            for stream in STREAMS
        ],
        'output': 'infiles/2023BPix/2023BPix_ParkingVBFData.json',
        'key_prefix': 'ParkingVBF_Run2023D',
    }
}


def query_das_files(dataset):
    """Query DAS for files in a dataset using dasgoclient."""
    cmd = ['dasgoclient', '--query', f'file dataset={dataset}', '--limit', '0']
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"  DAS query failed for {dataset}: {result.stderr.strip()}")
            return []
        files = [f.strip() for f in result.stdout.strip().split('\n') if f.strip()]
        return [XROOTD_PREFIX + f for f in files]
    except FileNotFoundError:
        print("ERROR: dasgoclient not found. Please set up CMSSW environment:")
        print("  source /cvmfs/cms.cern.ch/cmsset_default.sh && cmsenv")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(f"  Timeout querying {dataset}")
        return []


def main():
    for period, config in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"Generating ParkingVBF infiles for {period}")
        print(f"{'='*60}")

        result = {}

        for dataset in config['datasets']:
            # Extract stream and version from dataset name
            # e.g. /ParkingVBF0/Run2023C-22Sep2023_v3-v1/NANOAOD
            parts = dataset.strip('/').split('/')
            stream = parts[0]       # e.g. ParkingVBF0
            era_version = parts[1]  # e.g. Run2023C-22Sep2023_v3-v1
            # Make a clean key: ParkingVBF0_Run2023C-v3
            run_part = era_version.split('-')[0]  # Run2023C
            ver_part = era_version.split('_v')[1].split('-')[0]  # 3
            key = f"{stream}_{run_part}-v{ver_part}"

            print(f"  Querying {dataset} ...", flush=True)
            files = query_das_files(dataset)

            if files:
                print(f"  -> Found {len(files)} files")
                result[key] = files
            else:
                print(f"  -> No files found (dataset may not exist)")

        # Save output
        os.makedirs(os.path.dirname(config['output']), exist_ok=True)
        with open(config['output'], 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\nSaved {len(result)} datasets to {config['output']}")

        # Summary
        total_files = sum(len(v) for v in result.values())
        print(f"Total files: {total_files}")


if __name__ == '__main__':
    main()
