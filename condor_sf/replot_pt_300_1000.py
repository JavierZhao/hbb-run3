#!/usr/bin/env python3
"""
Re-plot the 2D trigger efficiency for 2022 ggF with jet pT range 300-1000 GeV.
Uses the already-saved .coffea output — no need to re-run the processor.

Usage (from condor_sf/):
    python replot_pt_300_1000.py
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import hist
from coffea import util
from trig_eff_2d import trig_vars_2d, plot_2d_trigger_efficiency

# Override the pT axis: 300-1000 GeV, 17 bins (~41 GeV each, matches original binning)
trig_vars_2d['pt_vs_msd']['axis_y'] = hist.axis.Regular(
    bins=17, start=300, stop=1000,
    name="pt", label="Leading Jet $p_{T}$ [GeV]"
)

# Load the already-saved output
output_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "trig_eff_2d_2022", "output", "trig_eff_2d_2022_ggF.coffea"
)
print(f"Loading {output_path}")
output = util.load(output_path)

# Re-plot
save_dir = os.path.join("figures_2d_reranged", "2022", "ggF")
plot_2d_trigger_efficiency(
    output,
    trig_vars_2d,
    save_dir=save_dir,
    year=2022,
    data=False,
    show_img=False
)
print(f"Plots saved to {save_dir}/")
