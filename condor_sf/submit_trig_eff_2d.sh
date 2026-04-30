#!/usr/bin/env bash
# submit_trig_eff_2d.sh — wrapper script to submit 2D trigger efficiency Condor job
# Creates necessary directories before submission

set -euo pipefail

# Default year
YEAR="${1:-2023}"

# Validate year
if [[ ! "$YEAR" =~ ^(2022|2022EE|2023|2023BPix)$ ]]; then
  echo "Error: Invalid year '$YEAR'"
  echo "Usage: $0 [YEAR]"
  echo "  YEAR: 2022, 2022EE, 2023, or 2023BPix (default: 2023)"
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/trig_eff_2d_${YEAR}"

echo "=========================================="
echo "Submitting 2D trigger efficiency job"
echo "Year: ${YEAR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "=========================================="

# Create output directories
echo "Creating output directories..."
mkdir -p "${OUTPUT_DIR}/figures_2d"
mkdir -p "${OUTPUT_DIR}/output"
mkdir -p "${OUTPUT_DIR}"

echo "Directories created successfully."
echo ""

# Submit the job
echo "Submitting Condor job..."
cd "${SCRIPT_DIR}"
condor_submit trig_eff_2d-condor.jdl YEAR="${YEAR}"

echo ""
echo "Job submitted successfully!"
echo "Monitor with: condor_q"
echo "Check logs in: ${OUTPUT_DIR}/"
