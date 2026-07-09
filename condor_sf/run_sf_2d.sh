#!/usr/bin/env bash
# run_sf_2d.sh — wrapper to run scale_factor_2d.py (2D scale factor calculation)

set -euo pipefail

# Defaults (can be overridden by CLI)
YEAR="2022"
TEST_MODE=""
MC_SAMPLE=""        # ttbar | qcd_muenriched (empty = use script default)
TXBB_REGION=""      # fail | pass_bb | pass_cc | pass | inclusive | all
EXTRA_FLAGS=""      # space-separated extra flags forwarded to python (e.g. --no-recoverable)
PYTHON_BIN="${PYTHON_BIN:-python3}"   # optional env override

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--year) YEAR="$2"; shift 2 ;;
    --test) TEST_MODE="--test"; shift ;;
    --quick-test) TEST_MODE="--quick-test"; shift ;;
    --mc-sample) MC_SAMPLE="--mc-sample $2"; shift 2 ;;
    --txbb-region) TXBB_REGION="--txbb-region $2"; shift 2 ;;
    --no-recoverable) EXTRA_FLAGS="${EXTRA_FLAGS} --no-recoverable"; shift ;;
    --workers) EXTRA_FLAGS="${EXTRA_FLAGS} --workers $2"; shift 2 ;;
    --mc-keys-filter) EXTRA_FLAGS="${EXTRA_FLAGS} --mc-keys-filter $2"; shift 2 ;;
    --mc-output-suffix) EXTRA_FLAGS="${EXTRA_FLAGS} --mc-output-suffix $2"; shift 2 ;;
    -p|--python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [-y|--year 2022] [--test] [--mc-sample ttbar|qcd_muenriched] [--txbb-region all|fail|pass_bb|pass_cc] [-p|--python python3]"
      exit 0 ;;
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
"${PYTHON_BIN}" "${SCRIPT_DIR}/scale_factor_2d.py" --year "${YEAR}" ${TEST_MODE} ${MC_SAMPLE} ${TXBB_REGION} ${EXTRA_FLAGS} 2>&1 | tee "${LOGFILE}"
PYRC=${PIPESTATUS[0]}
set +x

# Build the summary doc as the last step of the job.
# `pip install --user markdown` is a small pure-Python wheel (no system deps);
# if no network is available the script still emits summary.md without HTML.
#
# Only meaningful when the job covers the FULL scope (all 3 regions × both
# MC samples) — i.e. when MC_SAMPLE is empty/"both" AND TXBB_REGION is
# empty/"all". For per-region split sub-jobs we'd be summarizing just one
# combo (and six concurrent sub-jobs would race to overwrite each other's
# summary.md), so skip and let the user run summarize_year.sh post-merge.
_is_full_scope=true
[[ -n "${MC_SAMPLE}"   && ! "${MC_SAMPLE}"   =~ both ]] && _is_full_scope=false
[[ -n "${TXBB_REGION}" && ! "${TXBB_REGION}" =~ all  ]] && _is_full_scope=false

if [[ "${PYRC}" -eq 0 ]] && [[ -d "figures_sf_2d" ]] && ${_is_full_scope}; then
  echo ">>> Building summary index ..."
  "${PYTHON_BIN}" -m pip install --user --quiet markdown 2>/dev/null || \
    echo ">>> (python-markdown unavailable; will emit summary.md only, no HTML)"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/summarize_outputs.py" \
    --base-outdir . \
    --out figures_sf_2d/summary.md \
    --pdf 2>&1 | tee -a "${LOGFILE}"
elif ! ${_is_full_scope}; then
  echo ">>> Partial-scope run (MC_SAMPLE=${MC_SAMPLE} TXBB_REGION=${TXBB_REGION})."
  echo ">>> Skipping summary build — run summarize_year.sh ${YEAR} after all sub-jobs complete."
fi

echo ">>> Done. Processed both VBF and ggF production modes for data and MC."
