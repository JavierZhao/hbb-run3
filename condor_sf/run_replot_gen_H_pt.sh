#!/bin/bash
# Replot gen_H_pt trigger efficiency (pT >= 300 GeV) for all 4 periods.
# Run from condor_sf/:
#     bash run_replot_gen_H_pt.sh

set -e
cd "$(dirname "$0")"

echo "=============================================="
echo " Replotting gen_H_pt efficiency (300–1200 GeV)"
echo "=============================================="

# Map: year label -> expected coffea file (used to skip missing periods)
declare -A COFFEA_FILES
COFFEA_FILES[2022]="eff_22/output/trig_soup_eff_2022_ggF.coffea"
COFFEA_FILES[2022EE]="eff_22ee/output/trig_soup_eff_2022EE_ggF.coffea"
COFFEA_FILES[2023]="eff_23/output/trig_soup_eff_2023_ggF.coffea"
COFFEA_FILES[2023BPix]="eff_23bpix/output/trig_soup_eff_2023BPix_ggF.coffea"

for YEAR in 2022 2022EE 2023 2023BPix; do
    echo ""
    echo "=============================================="
    echo " ${YEAR}"
    echo "=============================================="

    if [ ! -f "${COFFEA_FILES[$YEAR]}" ]; then
        echo "Skipping ${YEAR} — ${COFFEA_FILES[$YEAR]} not found"
        continue
    fi

    python3 replot_gen_H_pt.py --year "${YEAR}" --mode both

    if [ $? -eq 0 ]; then
        echo "  -> done. Plots in figures_gen_H_pt_300/${YEAR}/"
    else
        echo "  -> ERROR during ${YEAR}. Continuing..."
    fi
done

echo ""
echo "=============================================="
echo " All periods complete"
echo "=============================================="
echo ""
echo "Output files:"
find figures_gen_H_pt_300 -name "*.png" 2>/dev/null || echo "  (none found)"
