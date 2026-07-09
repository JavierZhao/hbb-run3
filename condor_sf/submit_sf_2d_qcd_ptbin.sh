#!/usr/bin/env bash
# Fan out the QCD muon-enriched SF job into per-pT-bin sub-jobs (12 bins × N regions).
#
# Why this exists: the monolithic qcd_muenriched job for 2024 hit LPC's 2-day
# SYSTEM_PERIODIC_REMOVE cap with the data side still 2 eras short. The QCD
# MuEnriched MC is ~10x the size of TTtoLNu2Q (8338 vs 790 files) and the
# bottleneck is sequential per-file preprocessing on a single worker.
#
# This script reads the year's QCD_MuEnriched.json, queries its top-level keys
# (one per pT bin: QCD_MuEnriched_PT_15to20, ..., QCD_MuEnriched_PT_1000), and
# submits one job per (pT bin, txbb_region). Each sub-job uses:
#   --mc-keys-filter '^<key>$'         (limit MC to that single pT bin)
#   --mc-output-suffix _<bin>           (so coffea/CSV/figures don't clobber)
#
# After all sub-jobs land, merge with merge_sf_2d_qcd_ptbin.py to produce the
# canonical per-region coffea + CSV + figures.
#
# Usage:
#   bash submit_sf_2d_qcd_ptbin.sh 2024 fail,pass_bb
#   bash submit_sf_2d_qcd_ptbin.sh 2024 all                # = fail,pass_bb,pass_cc
#   bash submit_sf_2d_qcd_ptbin.sh 2024 fail               # one region
#   bash submit_sf_2d_qcd_ptbin.sh 2024 fail PT_300to470   # one region, one bin
#
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <year> [regions=all] [pt_bin_grep=*]"
  echo "  regions: comma-separated subset of {fail,pass_bb,pass_cc} or 'all'"
  echo "  pt_bin_grep: substring filter on bin label (e.g. PT_300to470). Default: all 12 bins."
  exit 2
fi

YEAR="$1"
REGIONS_ARG="${2:-all}"
BIN_FILTER="${3:-}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
JSON="${SCRIPT_DIR}/infiles/${YEAR}/${YEAR}_QCD_MuEnriched.json"
if [[ ! -f "$JSON" ]]; then
  echo "ERROR: $JSON not found"; exit 2
fi

# Resolve regions
if [[ "$REGIONS_ARG" == "all" ]]; then
  REGIONS=(fail pass_bb pass_cc)
else
  IFS=',' read -ra REGIONS <<< "$REGIONS_ARG"
fi

# Resolve pT-bin keys from the JSON (one key per bin).
mapfile -t ALL_BINS < <(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('\n'.join(d.keys()))" "$JSON")
if [[ -n "$BIN_FILTER" ]]; then
  FILTERED=()
  for b in "${ALL_BINS[@]}"; do
    [[ "$b" == *"$BIN_FILTER"* ]] && FILTERED+=("$b")
  done
  BINS=("${FILTERED[@]}")
else
  BINS=("${ALL_BINS[@]}")
fi

if [[ ${#BINS[@]} -eq 0 ]]; then
  echo "ERROR: no pT bins matched filter '$BIN_FILTER'"; exit 2
fi

echo "[plan] year=$YEAR  regions=${REGIONS[*]}  bins=${#BINS[@]}"
for b in "${BINS[@]}"; do echo "       $b"; done

cd "$SCRIPT_DIR"

# Submit
SUBMITTED=()
for REG in "${REGIONS[@]}"; do
  for KEY in "${BINS[@]}"; do
    # suffix: replace QCD_MuEnriched_ prefix with leading _ so paths read sf_2d_2024_..._qcd_muenriched_PT_300to470_*
    SUFFIX="_${KEY#QCD_MuEnriched_}"
    REGEX="^${KEY}$"
    echo "[submit] YEAR=$YEAR  TXBB_REGION=$REG  --mc-keys-filter='$REGEX'  --mc-output-suffix='$SUFFIX'"
    CL=$(condor_submit sf_2d-condor.jdl \
        YEAR="$YEAR" \
        MC_SAMPLE="qcd_muenriched" \
        TXBB_REGION="$REG" \
        EXTRA_FLAGS="--mc-keys-filter $REGEX --mc-output-suffix $SUFFIX" \
        2>&1 | tee /dev/stderr | awk '/submitted to cluster/ {gsub(/\./, ""); print $NF}')
    SUBMITTED+=("$YEAR/$REG/$KEY: $CL")
  done
done

echo
echo "[done] Submitted ${#SUBMITTED[@]} pT-bin sub-jobs for $YEAR:"
for s in "${SUBMITTED[@]}"; do echo "  $s"; done
