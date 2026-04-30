#!/usr/bin/env bash
# run_all_analyses.sh — Master script to generate all plots and tables for the Trigger AN
#
# This script runs all necessary analyses to produce:
#   1. 2D trigger scale factor maps (pT vs mSD) for VBF and ggF
#   2. 1D trigger scale factor curves (pT, mSD, HT, etc.)
#   3. Per-trigger efficiency tables (CSV format) for each period
#
# Usage:
#   ./run_all_analyses.sh [--test]
#
# Options:
#   --test    Run in test mode (process only 10 files per dataset)

set -euo pipefail

# Parse arguments
TEST_MODE=""
if [[ $# -gt 0 ]] && [[ "$1" == "--test" ]]; then
    TEST_MODE="--test"
    echo ">>> Running in TEST MODE (10 files per dataset)"
fi

# Get script directory
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "${SCRIPT_DIR}"

# All data-taking periods for Run 3
PERIODS=("2022" "2022EE" "2023" "2023BPix")

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Trigger Efficiency AN - Master Analysis${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "This script will generate all plots and tables for the trigger efficiency AN:"
echo "  - 2D scale factor maps (Data/MC) for pT vs mSD"
echo "  - 1D scale factor curves for pT, mSD, HT, etc."
echo "  - Per-trigger efficiency tables (CSV) for each period"
echo ""

# Summary of what will be produced
echo -e "${YELLOW}Output Summary:${NC}"
echo "  2D SF plots:  figures_sf_2d/{YEAR}/{MODE}/scale_factor_pt_vs_msd_2d.png"
echo "  1D SF plots:  sf_{YEAR}_OR/scale_factor_{VAR}_{YEAR}.png"
echo "  CSV tables:   output/{YEAR}/{MODE}/per_trigger_efficiency_{YEAR}_{MODE}.csv"
echo "  Coffea files: output/sf_2d_{YEAR}_{MODE}_{data,mc}.coffea"
echo ""
echo -e "${YELLOW}Data-taking periods:${NC} ${PERIODS[*]}"
echo ""

read -p "Press Enter to continue or Ctrl+C to abort..."

# Track timing
START_TIME=$(date +%s)

#==============================================================================
# Part 1: 2D Scale Factor Analysis
#==============================================================================
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Part 1: 2D Scale Factor Analysis${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "This produces:"
echo "  - 2D efficiency maps (pT vs mSD) for Data and MC"
echo "  - 2D scale factor maps (Data/MC ratio)"
echo "  - Per-trigger efficiency tables (CSV) with pT ≥ 450 GeV cut"
echo ""

for YEAR in "${PERIODS[@]}"; do
    echo ""
    echo -e "${BLUE}>>> Processing 2D SF for period: ${YEAR}${NC}"

    if [[ -x "${SCRIPT_DIR}/run_sf_2d.sh" ]]; then
        "${SCRIPT_DIR}/run_sf_2d.sh" --year "${YEAR}" ${TEST_MODE}
    else
        echo "Error: run_sf_2d.sh not found or not executable"
        exit 1
    fi

    echo -e "${GREEN}✓ Completed 2D SF analysis for ${YEAR}${NC}"
done

#==============================================================================
# Part 2: 1D Scale Factor Analysis
#==============================================================================
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Part 2: 1D Scale Factor Analysis${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "This produces:"
echo "  - 1D efficiency turn-on curves vs pT, mSD, HT, num_ak4, etc."
echo "  - 1D scale factor curves (Data/MC ratio)"
echo ""

for YEAR in "${PERIODS[@]}"; do
    echo ""
    echo -e "${BLUE}>>> Processing 1D SF for period: ${YEAR}${NC}"

    if [[ -x "${SCRIPT_DIR}/run_sf.sh" ]]; then
        "${SCRIPT_DIR}/run_sf.sh" --year "${YEAR}"
    else
        echo "Error: run_sf.sh not found or not executable"
        exit 1
    fi

    echo -e "${GREEN}✓ Completed 1D SF analysis for ${YEAR}${NC}"
done

#==============================================================================
# Summary
#==============================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All Analyses Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "Total time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo ""

# Print summary of outputs
echo -e "${YELLOW}Output Locations:${NC}"
echo ""
echo "2D Scale Factor Plots (for each period and production mode):"
for YEAR in "${PERIODS[@]}"; do
    echo "  - figures_sf_2d/${YEAR}/VBF/scale_factor_pt_vs_msd_2d.png"
    echo "  - figures_sf_2d/${YEAR}/ggF/scale_factor_pt_vs_msd_2d.png"
done
echo ""

echo "Per-Trigger Efficiency Tables (CSV):"
for YEAR in "${PERIODS[@]}"; do
    echo "  - output/${YEAR}/VBF/per_trigger_efficiency_${YEAR}_VBF.csv"
    echo "  - output/${YEAR}/ggF/per_trigger_efficiency_${YEAR}_ggF.csv"
done
echo ""

echo "1D Scale Factor Plots (combined OR of all triggers):"
for YEAR in "${PERIODS[@]}"; do
    echo "  - sf_${YEAR}_OR/*.png (pT, mSD, HT, num_ak4, particleNet_XbbVsQCD)"
done
echo ""

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Analysis Note Ready Outputs:${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "1. 2D efficiency maps:    Use figures_sf_2d/{YEAR}/{MODE}/"
echo "2. Per-trigger tables:    Use output/{YEAR}/{MODE}/per_trigger_efficiency_*.csv"
echo "3. 1D efficiency curves:  Use sf_{YEAR}_OR/"
echo ""
echo "For the AN, you can now:"
echo "  - Include 2D SF maps showing Data eff, MC eff, and SF"
echo "  - Include per-trigger efficiency tables (similar to Tables 2-4 in Run 2 AN)"
echo "  - Include 1D turn-on curves for key variables"
echo ""
echo -e "${GREEN}Done!${NC}"
