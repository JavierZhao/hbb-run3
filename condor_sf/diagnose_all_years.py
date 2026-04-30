"""
Diagnose the low data efficiency by examining saved coffea outputs.
Computes efficiency in different mSD windows to identify root cause.
"""
import numpy as np
from coffea import util
import os

years = {
    '2022':     'year_2022/output/sf_2d_2022_Inclusive',
    '2022EE':   'year_2022EE/output/sf_2d_2022EE_Inclusive',
    '2023':     'year_2023/output/sf_2d_2023_Inclusive',
    '2023BPix': 'year_2023BPix/output/sf_2d_2023BPix_Inclusive',
}

var = 'pt_vs_msd'

# mSD windows to compare
msd_windows = [
    ('All mSD',       0,   300),
    ('mSD < 40',      0,    40),
    ('mSD 40-80',    40,    80),
    ('mSD 80-150',   80,   150),   # Higgs mass window
    ('mSD > 150',   150,   300),
]

pt_cut = 300  # minimum pT

print(f"{'Year':10s} {'Window':15s} {'Data':>10s} {'MC':>10s} {'SF':>8s}  (pT > {pt_cut} GeV)")
print('-' * 65)

for year, base_path in years.items():
    data_path = base_path + '_data.coffea'
    mc_path   = base_path + '_mc.coffea'

    if not os.path.exists(data_path):
        print(f"{year}: output not found at {data_path}")
        continue

    d = util.load(data_path)
    m = util.load(mc_path)

    base_x_d = np.array(d['Baseline'][var]['x'])  # mSD
    base_y_d = np.array(d['Baseline'][var]['y'])  # pT
    pass_x_d = np.array(d['Numerator'][var]['pass_x'])
    pass_y_d = np.array(d['Numerator'][var]['pass_y'])

    base_x_m = np.array(m['Baseline'][var]['x'])
    base_y_m = np.array(m['Baseline'][var]['y'])
    pass_x_m = np.array(m['Numerator'][var]['pass_x'])
    pass_y_m = np.array(m['Numerator'][var]['pass_y'])

    for label, msd_lo, msd_hi in msd_windows:
        cut_base_d = (base_y_d >= pt_cut) & (base_x_d >= msd_lo) & (base_x_d < msd_hi)
        cut_pass_d = (pass_y_d >= pt_cut) & (pass_x_d >= msd_lo) & (pass_x_d < msd_hi)
        cut_base_m = (base_y_m >= pt_cut) & (base_x_m >= msd_lo) & (base_x_m < msd_hi)
        cut_pass_m = (pass_y_m >= pt_cut) & (pass_x_m >= msd_lo) & (pass_x_m < msd_hi)

        n_base_d = np.sum(cut_base_d);  n_pass_d = np.sum(cut_pass_d)
        n_base_m = np.sum(cut_base_m);  n_pass_m = np.sum(cut_pass_m)

        eff_d = n_pass_d / n_base_d if n_base_d > 0 else 0
        eff_m = n_pass_m / n_base_m if n_base_m > 0 else 0
        sf    = eff_d / eff_m if eff_m > 0 else 0

        print(f"{year:10s} {label:15s} {eff_d:9.1%}  {eff_m:9.1%}  {sf:7.3f}")

    print()
