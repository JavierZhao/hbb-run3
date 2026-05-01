#!/usr/bin/env bash
# Run ParkingVBF trigger efficiency study
# Usage: ./run_parking_vbf_eff.sh [--year 2023|2023BPix] [--test]

set -euo pipefail

YEAR="2023"
METHOD="zmumu_vbf"
TEST_FLAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --year)    YEAR="$2";   shift 2 ;;
        --method)  METHOD="$2"; shift 2 ;;
        --test)    TEST_FLAG="--test"; shift ;;
        *)         echo "Unknown argument: $1"; exit 1 ;;
    esac
done

echo "Running VBF trigger study: year=${YEAR} method=${METHOD} ${TEST_FLAG}"
python3 parking_vbf_eff.py --year "${YEAR}" --method "${METHOD}" ${TEST_FLAG}
