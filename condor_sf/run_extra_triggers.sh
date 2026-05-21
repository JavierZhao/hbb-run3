#!/usr/bin/env bash
# run_extra_triggers.sh — wrapper for study_extra_triggers.py.
# Usage examples:
#   bash run_extra_triggers.sh --year 2022 --check-only
#   bash run_extra_triggers.sh --year 2022 --test
#   bash run_extra_triggers.sh --year 2022 --mc-sample qcd_muenriched
set -euo pipefail

YEAR="2022"
MC_SAMPLE="qcd_muenriched"
OUT_BASE="figures_extra_trig"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--year)        YEAR="$2"; shift 2 ;;
    --mc-sample)      MC_SAMPLE="$2"; shift 2 ;;
    --out-base)       OUT_BASE="$2"; shift 2 ;;
    --check-only)     EXTRA_ARGS+=("--check-only"); shift ;;
    --test)           EXTRA_ARGS+=("--test"); shift ;;
    --chunksize)      EXTRA_ARGS+=("--chunksize" "$2"); shift 2 ;;
    --workers)        EXTRA_ARGS+=("--workers" "$2"); shift 2 ;;
    --quick-files)    EXTRA_ARGS+=("--quick-files" "$2"); shift 2 ;;
    -h|--help)
      cat <<EOF
Usage: $0 [-y YEAR] [--mc-sample ttbar|qcd_muenriched]
          [--check-only] [--test] [--chunksize N] [--workers N]

Year choices: 2022, 2022EE, 2023, 2023BPix
EOF
      exit 0 ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
RUN_DIR="${SCRIPT_DIR}/extra_trig_${YEAR}"
mkdir -p "${RUN_DIR}"

export MPLBACKEND=Agg
export MPLCONFIGDIR="${RUN_DIR}/.mplconfig"
mkdir -p "${MPLCONFIGDIR}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TS="$(date +%Y%m%d_%H%M%S)"
LOGFILE="${RUN_DIR}/run_${YEAR}_${MC_SAMPLE}_${TS}.log"

echo ">>> extra-trigger study: year=${YEAR}, mc=${MC_SAMPLE}"
echo ">>> log: ${LOGFILE}"

# Always make sure the out-base directory exists so condor's transfer_output_files
# never sees a missing path (e.g. in --check-only mode where the python script
# doesn't write any figures).
mkdir -p "${OUT_BASE}"

set -x
python3 "${SCRIPT_DIR}/study_extra_triggers.py" \
    --year "${YEAR}" \
    --mc-sample "${MC_SAMPLE}" \
    --out-base "${OUT_BASE}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee "${LOGFILE}"
