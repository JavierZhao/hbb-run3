#!/usr/bin/env bash
# submit_all_analyses.sh — Submit condor jobs for all trigger efficiency analyses
#
# This script provides options for submitting analyses:
#   1. Sequential: One job that runs all periods sequentially
#   2. Parallel: Four jobs, one per period (faster)
#
# Usage:
#   ./submit_all_analyses.sh [--parallel] [--test] [--year YEAR]
#
# Options:
#   --parallel    Submit separate jobs for each year (faster, default)
#   --sequential  Submit one job for all years (slower)
#   --test        Run in test mode (10 files only)
#   --year YEAR   Submit only for specific year (with --parallel)

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "${SCRIPT_DIR}"

# Default options
MODE="parallel"
TEST_FLAG=""
SPECIFIC_YEAR=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            MODE="parallel"
            shift
            ;;
        --sequential)
            MODE="sequential"
            shift
            ;;
        --test)
            TEST_FLAG="--test"
            shift
            ;;
        --year)
            SPECIFIC_YEAR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--parallel|--sequential] [--test] [--year YEAR]"
            echo ""
            echo "Options:"
            echo "  --parallel      Submit separate jobs for each year (default, faster)"
            echo "  --sequential    Submit one job for all years (slower)"
            echo "  --test          Run in test mode (10 files only)"
            echo "  --year YEAR     Submit only for specific year (with --parallel)"
            echo ""
            echo "Examples:"
            echo "  $0                           # Submit 4 parallel jobs (one per year)"
            echo "  $0 --sequential              # Submit 1 job for all years"
            echo "  $0 --year 2022               # Submit 1 job for 2022 only"
            echo "  $0 --test                    # Submit 4 parallel test jobs"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run with --help for usage"
            exit 1
            ;;
    esac
done

echo "========================================"
echo "Condor Job Submission"
echo "========================================"
echo "Mode: $MODE"
if [[ -n "$TEST_FLAG" ]]; then
    echo "Test mode: ENABLED (10 files only)"
fi
if [[ -n "$SPECIFIC_YEAR" ]]; then
    echo "Year: $SPECIFIC_YEAR"
fi
echo ""

# Create output directories
if [[ "$MODE" == "parallel" ]]; then
    if [[ -n "$SPECIFIC_YEAR" ]]; then
        mkdir -p "year_${SPECIFIC_YEAR}"
        echo "Output will go to: year_${SPECIFIC_YEAR}/"
    else
        mkdir -p year_2022 year_2022EE year_2023 year_2023BPix
        echo "Output will go to: year_<YEAR>/"
    fi
else
    mkdir -p all_analyses
    echo "Output will go to: all_analyses/"
fi

echo ""
read -p "Press Enter to submit jobs or Ctrl+C to abort..."

# Submit based on mode
if [[ "$MODE" == "sequential" ]]; then
    echo ""
    echo "Submitting sequential job (all years in one job)..."

    if [[ -n "$TEST_FLAG" ]]; then
        condor_submit run_all-condor.jdl TEST_FLAG="$TEST_FLAG"
    else
        condor_submit run_all-condor.jdl
    fi

    echo ""
    echo "Submitted 1 job"
    echo "This will run all 4 periods sequentially (slow but uses less resources)"

elif [[ "$MODE" == "parallel" ]]; then
    echo ""
    if [[ -n "$SPECIFIC_YEAR" ]]; then
        echo "Submitting parallel job for year: $SPECIFIC_YEAR..."

        if [[ -n "$TEST_FLAG" ]]; then
            condor_submit run_year-condor.jdl YEAR="$SPECIFIC_YEAR" TEST_FLAG="$TEST_FLAG"
        else
            condor_submit run_year-condor.jdl YEAR="$SPECIFIC_YEAR"
        fi

        echo ""
        echo "Submitted 1 job for $SPECIFIC_YEAR"
    else
        echo "Submitting parallel jobs (one per year)..."

        if [[ -n "$TEST_FLAG" ]]; then
            condor_submit run_year-condor.jdl TEST_FLAG="$TEST_FLAG"
        else
            condor_submit run_year-condor.jdl
        fi

        echo ""
        echo "Submitted 4 jobs (one per period: 2022, 2022EE, 2023, 2023BPix)"
        echo "This will run all periods in parallel (fast, recommended)"
    fi
fi

echo ""
echo "========================================"
echo "Monitor jobs with:"
echo "  condor_q"
echo "  condor_q -nobatch"
echo ""
echo "Check logs:"
if [[ "$MODE" == "sequential" ]]; then
    echo "  tail -f all_analyses/<cluster_id>_0.out"
elif [[ -n "$SPECIFIC_YEAR" ]]; then
    echo "  tail -f year_${SPECIFIC_YEAR}/<cluster_id>_0.out"
else
    echo "  tail -f year_2022/<cluster_id>_0.out"
    echo "  tail -f year_2022EE/<cluster_id>_0.out"
    echo "  tail -f year_2023/<cluster_id>_0.out"
    echo "  tail -f year_2023BPix/<cluster_id>_0.out"
fi
echo ""
echo "Done!"
