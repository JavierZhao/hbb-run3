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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

echo "Running VBF trigger study: year=${YEAR} method=${METHOD} ${TEST_FLAG}"
python3 parking_vbf_eff.py --year "${YEAR}" --method "${METHOD}" ${TEST_FLAG}
PYRC=$?

# Build the summary doc as the last step of the job.
if [[ "${PYRC}" -eq 0 ]] && [[ -d "figures_parking_vbf" ]]; then
  echo ">>> Building summary index ..."
  python3 -m pip install --user --quiet markdown 2>/dev/null || \
    echo ">>> (python-markdown unavailable; will emit summary.md only, no HTML)"
  python3 "${SCRIPT_DIR}/summarize_outputs.py" \
    --base-outdir . \
    --out figures_parking_vbf/summary.md \
    --pdf
fi
