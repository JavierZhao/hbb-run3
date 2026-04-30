#!/bin/bash
# Drop into a coffea-dask container shell from condor_sf/.
#
# Usage:
#   bash dask_shell.sh
#   ./dask_shell.sh

set -euo pipefail
cd "$(dirname "$0")"

# Refresh VOMS proxy if missing or expired
if ! voms-proxy-check 2>/dev/null; then
    echo "Proxy missing or expired — running voms-proxy-init ..."
    voms-proxy-init -voms cms
fi

exec ./shell coffeateam/coffea-dask:latest-py3.9 bash
