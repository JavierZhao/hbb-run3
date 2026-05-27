#!/usr/bin/env bash
# Wrapper for study_signal_eff.py.
# Usage examples:
#   bash run_signal_eff.sh --year 2022 --sample both --quick-files 10
#   bash run_signal_eff.sh --year 2023 --sample ggF
set -euo pipefail

YEAR="2022"
SAMPLE="both"
OUT_BASE="figures_signal_eff"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -y|--year)        YEAR="$2"; shift 2 ;;
    --sample)         SAMPLE="$2"; shift 2 ;;
    --out-base)       OUT_BASE="$2"; shift 2 ;;
    --quick-files)    EXTRA_ARGS+=("--quick-files" "$2"); shift 2 ;;
    --probe-pt-min)   EXTRA_ARGS+=("--probe-pt-min" "$2"); shift 2 ;;
    --probe-ht-min)   EXTRA_ARGS+=("--probe-ht-min" "$2"); shift 2 ;;
    --candidates)     EXTRA_ARGS+=("--candidates" "$2"); shift 2 ;;
    --chunksize)      EXTRA_ARGS+=("--chunksize" "$2"); shift 2 ;;
    --workers)        EXTRA_ARGS+=("--workers" "$2"); shift 2 ;;
    -h|--help)
      echo "Usage: $0 [-y YEAR] [--sample ggF|VBF|both] [--out-base DIR] [--quick-files N]"
      exit 0 ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
RUN_DIR="${SCRIPT_DIR}/signal_eff_${YEAR}"
mkdir -p "${RUN_DIR}"

export MPLBACKEND=Agg
export MPLCONFIGDIR="${RUN_DIR}/.mplconfig"
mkdir -p "${MPLCONFIGDIR}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${OUT_BASE}"

TS="$(date +%Y%m%d_%H%M%S)"
LOGFILE="${RUN_DIR}/run_${YEAR}_${SAMPLE}_${TS}.log"

echo ">>> signal-eff study: year=${YEAR}, sample=${SAMPLE}, out=${OUT_BASE}"
echo ">>> log: ${LOGFILE}"

set -x
python3 -u "${SCRIPT_DIR}/study_signal_eff.py" \
    --year "${YEAR}" \
    --sample "${SAMPLE}" \
    --out-base "${OUT_BASE}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee "${LOGFILE}"
