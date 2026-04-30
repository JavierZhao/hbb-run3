#!/usr/bin/env bash
# run_sf.sh — wrapper to run scale_factor.py with CLI args

set -euo pipefail

# Defaults (can be overridden by CLI)
YEAR="2022"
PYTHON_BIN="${PYTHON_BIN:-python3}"   # optional env override

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--year) YEAR="$2"; shift 2 ;;
  # trigger flag removed; kept for backward compatibility placeholder
  -t|--trigger) echo "Warning: trigger flag removed; ignoring value $2"; shift 2 ;;
    -p|--python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [-y|--year 2022] [-p|--python python3]"; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# Per-run directory (requested): sf_YEAR_OR
RUN_DIR="${SCRIPT_DIR}/sf_${YEAR}_OR"
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
LOGFILE="${RUN_DIR}/run_${YEAR}_OR_${TS}.log"

echo ">>> Running: year=${YEAR}, trigger=OR"
echo ">>> Run dir: ${RUN_DIR}"
echo ">>> Log: ${LOGFILE}"

set -x
"${PYTHON_BIN}" "${SCRIPT_DIR}/trig_eff.py" --year "${YEAR}" 2>&1 | tee "${LOGFILE}"
set +x

echo ">>> Done."
