#!/usr/bin/env bash
# submit_sf_2d.sh — wrapper script to submit 2D scale factor Condor job
# Creates necessary directories before submission

set -euo pipefail

# Default year
YEAR="${1:-2023}"

# Validate year
if [[ ! "$YEAR" =~ ^(2022|2022EE|2023|2023BPix|2024)$ ]]; then
  echo "Error: Invalid year '$YEAR'"
  echo "Usage: $0 [YEAR]"
  echo "  YEAR: 2022, 2022EE, 2023, 2023BPix, or 2024 (default: 2023)"
  exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/sf_2d_${YEAR}"

echo "=========================================="
echo "Submitting 2D scale factor job"
echo "Year: ${YEAR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "=========================================="

# Move any previous output aside BEFORE condor's file-transfer back-fills
# the destination. condor merges new files into the existing dir, so an old
# .coffea/.png from a previous submission survives if the new run fails to
# produce that specific (region, mc) combination. This bit us in 2023BPix on
# 2026-05-22 — the QCD_MuEnriched[fail] step crashed via try/except and the
# May-6 stale file silently fed back into summarize_outputs.py.
TS="$(date +%Y%m%d_%H%M%S)"
for sub in output figures_sf_2d; do
  if [[ -d "${OUTPUT_DIR}/${sub}" ]]; then
    BAK="${OUTPUT_DIR}/${sub}_prev_${TS}"
    echo "Stashing previous ${sub}: ${OUTPUT_DIR}/${sub} -> ${BAK}"
    mv "${OUTPUT_DIR}/${sub}" "${BAK}"
  fi
done

# Create fresh output directories
echo "Creating output directories..."
mkdir -p "${OUTPUT_DIR}/figures_sf_2d"
mkdir -p "${OUTPUT_DIR}/output"
mkdir -p "${OUTPUT_DIR}"

echo "Directories created successfully."
echo ""

# Submit the job (pass-through any extra args like MC_SAMPLE=ttbar)
echo "Submitting Condor job..."
cd "${SCRIPT_DIR}"
condor_submit sf_2d-condor.jdl YEAR="${YEAR}" "${@:2}"

echo ""
echo "Job submitted successfully!"
echo "Monitor with: condor_q"
echo "Check logs in: ${OUTPUT_DIR}/"
