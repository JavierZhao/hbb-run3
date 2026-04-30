#!/usr/bin/env python3
"""
Replot Gen Higgs pT trigger efficiency from saved .coffea outputs,
starting the pT axis at 300 GeV instead of 0.

Usage (from condor_sf/):
    python replot_gen_H_pt.py --year 2022 --mode ggF
    python replot_gen_H_pt.py --year 2023 --mode both
"""
import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import hist
from coffea import util
from trig_eff import plot_1d_trigger_soup_cms

# Map year labels to their eff output subdirectories
YEAR_TO_DIR = {
    '2022':     'eff_22',
    '2022EE':   'eff_22ee',
    '2023':     'eff_23',
    '2023BPix': 'eff_23bpix',
}

# Only the gen_H_pt variable, axis starting at 300
trig_vars_replot = {
    'gen_H_pt': {
        'label': "Gen Higgs pT [GeV]",
        'axis': hist.axis.Regular(bins=30, start=300, stop=1200,
                                  name="gen_H_pt", label="Gen Higgs pT [GeV]"),
    }
}


def run(year, mode):
    eff_dir = YEAR_TO_DIR[year]
    coffea_path = os.path.join(eff_dir, "output", f"trig_soup_eff_{year}_{mode}.coffea")

    if not os.path.isfile(coffea_path):
        print(f"WARNING: {coffea_path} not found — skipping {year} {mode}")
        return

    print(f"Loading {coffea_path}")
    output = util.load(coffea_path)

    save_dir = os.path.join("figures_gen_H_pt_300", year, mode)
    print(f"Plotting gen_H_pt (300–1200 GeV) → {save_dir}/")

    plot_1d_trigger_soup_cms(
        output,
        trig_vars_replot,
        save_dir=save_dir,
        year=int(year[:4]),   # mplhep wants an int; use the 4-digit year prefix
        data=False,
        show_img=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replot gen_H_pt efficiency starting at 300 GeV")
    parser.add_argument('--year', required=True, choices=YEAR_TO_DIR.keys(),
                        help="Data-taking period")
    parser.add_argument('--mode', required=True, choices=['ggF', 'VBF', 'both'],
                        help="Production mode (or 'both')")
    args = parser.parse_args()

    modes = ['ggF', 'VBF'] if args.mode == 'both' else [args.mode]
    for m in modes:
        run(args.year, m)
