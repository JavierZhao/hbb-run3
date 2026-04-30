#!/usr/bin/env bash
# run_year_combined.sh — Wrapper to run both 2D and 1D SF analyses for a single year
#
# Usage:
#   ./run_year_combined.sh --year 2022 [--test]

set -euo pipefail

# Parse arguments
YEAR=""
TEST_FLAG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -y|--year)
            YEAR="$2"
            shift 2
            ;;
        --test)
            TEST_FLAG="--test"
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$YEAR" ]]; then
    echo "Error: --year is required"
    exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

echo "========================================"
echo "Running Combined Analysis for ${YEAR}"
echo "========================================"
echo ""

# Run 2D SF analysis
echo ">>> Step 1/2: Running 2D Scale Factor Analysis..."
"${SCRIPT_DIR}/run_sf_2d.sh" --year "${YEAR}" ${TEST_FLAG}

echo ""
echo ">>> Step 2/2: Running 1D Scale Factor Analysis..."
"${SCRIPT_DIR}/run_sf.sh" --year "${YEAR}"

echo ""
echo "========================================"
echo "Completed Combined Analysis for ${YEAR}"
echo "========================================"
