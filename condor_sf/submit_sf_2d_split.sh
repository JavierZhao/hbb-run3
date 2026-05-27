#!/usr/bin/env bash
# Submit the SF pipeline as 6 sub-jobs per year (3 regions × 2 MCs) instead of
# one monolithic job. Each sub-job processes exactly one (txbb_region, mc_sample)
# combination, so its wall time is ~1/6 of the all-in-one job and well within
# condor's 48h limit.
#
# On 2026-05-20–22, the monolithic submission produced 48h timeouts for the 3
# largest years (2022, 2022EE, 2023) — outputs lost on ABORT. Splitting fixes
# this; if one combo fails the others still produce clean outputs.
#
# Stashes existing output/ and figures_sf_2d/ to *_prev_<TS>/ BEFORE submitting
# so condor's file transfer can't merge new files into a partially-stale tree.
#
# Usage:
#   bash submit_sf_2d_split.sh 2022
#   bash submit_sf_2d_split.sh 2023BPix
#   bash submit_sf_2d_split.sh 2022 EXTRA_FLAGS=--no-recoverable
#
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <year> [extra condor_submit args, e.g. EXTRA_FLAGS=--no-recoverable]"
  exit 2
fi

YEAR="$1"
shift

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
RUN_DIR="${SCRIPT_DIR}/sf_2d_${YEAR}"
mkdir -p "${RUN_DIR}"

# --- Stash any previous outputs ---
if [[ -d "${RUN_DIR}" ]]; then
  TS="$(date +%Y%m%d_%H%M%S)"
  for sub in output figures_sf_2d; do
    if [[ -d "${RUN_DIR}/${sub}" ]]; then
      BAK="${RUN_DIR}/${sub}_prev_${TS}"
      echo "[clean] mv ${RUN_DIR}/${sub} -> ${BAK}"
      mv "${RUN_DIR}/${sub}" "${BAK}"
    fi
  done
fi

cd "${SCRIPT_DIR}"

# --- Submit 6 sub-jobs (3 regions × 2 MCs) ---
SUBMITTED=()
for MC in ttbar qcd_muenriched; do
  for REG in fail pass_bb pass_cc; do
    echo "[submit] YEAR=${YEAR}  MC_SAMPLE=${MC}  TXBB_REGION=${REG}"
    CL=$(condor_submit sf_2d-condor.jdl \
        YEAR="${YEAR}" \
        MC_SAMPLE="${MC}" \
        TXBB_REGION="${REG}" \
        "$@" 2>&1 | tee /dev/stderr | awk '/submitted to cluster/ {gsub(/\./, ""); print $NF}')
    SUBMITTED+=("${YEAR}/${MC}/${REG}: ${CL}")
  done
done

echo ""
echo "[done] Submitted ${#SUBMITTED[@]} sub-jobs for ${YEAR}:"
for s in "${SUBMITTED[@]}"; do echo "  $s"; done
