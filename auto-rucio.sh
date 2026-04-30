#!/usr/bin/env bash
# resolve_and_add.sh
set -euo pipefail

RSE="${RSE:-T1_US_FNAL_Disk}"
LIFETIME="${LIFETIME:-14000000}"
ACTIVITY="${ACTIVITY:-User AutoApprove}"
ASK="${ASK:---ask-approval}"
COMMENT="${COMMENT:-}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing $1 in PATH"; exit 1; }; }
need dasgoclient; need rucio

# Provide *partials* below (wildcards OK). The script resolves them to exact DIDs.
PARTIALS=(
#   '/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-*/NANOAODSIM'
#   '/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-*/NANOAODSIM'
#   '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-*/NANOAODSIM'
)

for P in "${PARTIALS[@]}"; do
  echo "Resolving: $P"
  DID=$(dasgoclient -query="dataset=$P" | head -n 1 || true)
  if [[ -z "$DID" ]]; then
    echo "  ! No match in DAS: $P"
    continue
  fi
  echo "  -> $DID"
  rucio add-rule "cms:${DID}" 1 "$RSE" --activity "$ACTIVITY" --lifetime "$LIFETIME" ${ASK:+$ASK} --comment "$COMMENT"
done
