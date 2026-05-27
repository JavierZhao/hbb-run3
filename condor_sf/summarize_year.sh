#!/usr/bin/env bash
# Post-merge summary builder for the SF pipeline.
#
# Each per-region sub-job from submit_sf_2d_split.sh writes its output into the
# shared sf_2d_<YEAR>/{output,figures_sf_2d}/ trees but skips the
# summarize_outputs.py step (because its scope is partial). Once all 6 sub-jobs
# for a year have completed, run this to build the consolidated summary.
#
# Runs summarize_outputs.py inside the coffea container so it has the same
# environment as the sub-jobs. The `shell` helper binds the current working
# directory ($PWD = condor_sf/) to /srv, so both summarize_outputs.py and
# sf_2d_<YEAR>/ are accessible from inside.
#
# Usage:
#   bash summarize_year.sh 2022
#   bash summarize_year.sh 2023BPix
#   bash summarize_year.sh 2022 --no-pdf       # forward args to summarize_outputs.py
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <year> [extra args forwarded to summarize_outputs.py]"
  exit 2
fi

YEAR="$1"
shift

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
RUN_DIR="${SCRIPT_DIR}/sf_2d_${YEAR}"

if [[ ! -d "${RUN_DIR}" ]]; then
  echo "ERROR: ${RUN_DIR} does not exist."
  exit 3
fi
if [[ ! -d "${RUN_DIR}/figures_sf_2d" ]] && [[ ! -d "${RUN_DIR}/output" ]]; then
  echo "ERROR: ${RUN_DIR} has no figures_sf_2d/ or output/ — has the batch finished?"
  exit 4
fi

# Default extra args: build PDF if pandoc/weasyprint/wkhtmltopdf available.
EXTRA_ARGS="--pdf"
if [[ $# -gt 0 ]]; then EXTRA_ARGS="$*"; fi

export MPLBACKEND=Agg
export MPLCONFIGDIR="${RUN_DIR}/.mplconfig"
mkdir -p "${MPLCONFIGDIR}"

# `shell` binds $PWD -> /srv. cd to condor_sf/ so both summarize_outputs.py
# (at /srv/summarize_outputs.py) and the per-year output tree (at
# /srv/sf_2d_<YEAR>/) are reachable from inside the container.
cd "${SCRIPT_DIR}"
echo ">>> summarize_year ${YEAR}: condor_sf -> /srv inside container"

CONTAINER="coffeateam/coffea-dask-almalinux8:2025.10.2-py3.10"
"${SCRIPT_DIR}/shell" "${CONTAINER}" \
  "python3 summarize_outputs.py \
     --base-outdir sf_2d_${YEAR} \
     --out sf_2d_${YEAR}/figures_sf_2d/summary.md \
     ${EXTRA_ARGS}"

echo ">>> Wrote ${RUN_DIR}/figures_sf_2d/summary.{md,html,pdf}"
