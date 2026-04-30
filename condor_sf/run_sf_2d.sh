#!/usr/bin/env bash
# run_sf_2d.sh — wrapper to run scale_factor_2d.py (2D scale factor calculation)

set -euo pipefail

# Defaults (can be overridden by CLI)
YEAR="2022"
TEST_MODE=""
PYTHON_BIN="${PYTHON_BIN:-python3}"   # optional env override

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--year) YEAR="$2"; shift 2 ;;
    --test) TEST_MODE="--test"; shift ;;
    -p|--python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [-y|--year 2022] [--test] [-p|--python python3]"; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# Per-run directory: sf_2d_YEAR
RUN_DIR="${SCRIPT_DIR}/sf_2d_${YEAR}"
mkdir -p "${RUN_DIR}"

# Headless-safe plotting
export MPLBACKEND=Agg
export MPLCONFIGDIR="${RUN_DIR}/.mplconfig"
mkdir -p "${MPLCONFIGDIR}" "${RUN_DIR}"

# Keep threads modest on shared workers
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TS="$(date +%Y%m%d_%H%M%S)"
LOGFILE="${RUN_DIR}/run_${YEAR}_${TS}.log"

echo ">>> Running 2D scale factor calculation: year=${YEAR}"
if [[ -n "${TEST_MODE}" ]]; then
  echo ">>> TEST MODE enabled (10 files only)"
fi
echo ">>> Run dir: ${RUN_DIR}"
echo ">>> Log: ${LOGFILE}"

set -x
"${PYTHON_BIN}" "${SCRIPT_DIR}/scale_factor_2d.py" --year "${YEAR}" ${TEST_MODE} 2>&1 | tee "${LOGFILE}"
set +x

echo ">>> Done. Processed both VBF and ggF production modes for data and MC."
