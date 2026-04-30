#!/usr/bin/env bash
# rucio_bulk_add_rules.sh

set -euo pipefail

echo "Using bash: ${BASH_VERSION:-unknown}"
echo "Script: $0"
echo

# ---- Configurable knobs (overridable via env) ----
RSE="${RSE:-T1_US_FNAL_Disk}"
LIFETIME="${LIFETIME:-14000000}"
ACTIVITY="${ACTIVITY:-User AutoApprove}"
ASK_FLAG="${ASK_FLAG:---ask-approval}"   # set to "" to disable
COMMENT="${COMMENT:-}"                    # optional
DRYRUN="${DRYRUN:-0}"                     # 1 = don't execute, just print

# ---- Check rucio client exists ----
if ! command -v rucio >/dev/null 2>&1; then
  echo "ERROR: rucio CLI not found in PATH." >&2
  exit 1
fi

# ---- Dataset DIDs ----
DIDS="$(cat <<'EOF'
cms:/TTto4Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22NanoAODv12/NANOAODSIM
cms:/TTto2L2Nu_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22NanoAODv12/NANOAODSIM
cms:/TTtoLNu2Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22NanoAODv12/NANOAODSIM
cms:/TTto4Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22EENanoAODv12/NANOAODSIM
cms:/TTto2L2Nu_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22EENanoAODv12/NANOAODSIM
cms:/TTtoLNu2Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer22EENanoAODv12/NANOAODSIM
cms:/Muon0/Run2023C-22Sep2023_v1-v1/NANOAOD
cms:/Muon0/Run2023C-22Sep2023_v2-v1/NANOAOD
cms:/Muon0/Run2023C-22Sep2023_v3-v1/NANOAOD
cms:/Muon0/Run2023C-22Sep2023_v4-v1/NANOAOD
cms:/Muon1/Run2023C-22Sep2023_v1-v1/NANOAOD
cms:/Muon1/Run2023C-22Sep2023_v2-v1/NANOAOD
cms:/Muon1/Run2023C-22Sep2023_v3-v1/NANOAOD
cms:/Muon1/Run2023C-22Sep2023_v4-v1/NANOAOD
cms:/TTto4Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer23NanoAODv12/NANOAODSIM
cms:/TTto2L2Nu_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer23NanoAODv12/NANOAODSIM
cms:/TTtoLNu2Q_TuneCP5_13.6TeV_powheg-pythia8/Run3Summer23NanoAODv12/NANOAODSIM
cms:/Muon0/Run2023D-22Sep2023_v1-v1/NANOAOD
cms:/Muon0/Run2023D-22Sep2023_v2-v1/NANOAOD
cms:/Muon1/Run2023D-22Sep2023_v1-v1/NANOAOD
cms:/Muon1/Run2023D-22Sep2023_v2-v1/NANOAOD
EOF
)"

echo "Target RSE:     $RSE"
echo "Lifetime:       $LIFETIME"
echo "Activity:       $ACTIVITY"
echo "Ask approval:   ${ASK_FLAG:-<disabled>}"
echo "Comment:        ${COMMENT:-<empty>}"
echo "Dry-run:        $DRYRUN"
echo

# ---- Run submissions ----
while IFS= read -r DID; do
  [[ -z "$DID" ]] && continue
  echo ">>> Submitting: $DID"
  cmd=( rucio add-rule "$DID" 1 "$RSE" --activity "$ACTIVITY" --lifetime "$LIFETIME" )
  [[ -n "${ASK_FLAG}" ]] && cmd+=( "$ASK_FLAG" )
  cmd+=( --comment "$COMMENT" )
  echo "+ ${cmd[*]}"
  if [[ "$DRYRUN" != "1" ]]; then
    "${cmd[@]}"
  fi
  echo
done <<< "$DIDS"

echo "Done."
