"""
Diagnose the low 2023BPix data efficiency by examining the saved coffea outputs.
"""
import numpy as np
from coffea import util

output_data = util.load('year_2023BPix/output/sf_2d_2023BPix_Inclusive_data.coffea')  # pass region
output_mc   = util.load('year_2023BPix/output/sf_2d_2023BPix_Inclusive_mc.coffea')

print('=== DATA OUTPUT ===')
print('Keys:', list(output_data.keys()))

var = 'pt_vs_msd'

baseline_x_data = np.array(output_data['Baseline'][var]['x'])  # mSD
baseline_y_data = np.array(output_data['Baseline'][var]['y'])  # pT
pass_x_data     = np.array(output_data['Numerator'][var]['pass_x'])
pass_y_data     = np.array(output_data['Numerator'][var]['pass_y'])

baseline_x_mc = np.array(output_mc['Baseline'][var]['x'])
baseline_y_mc = np.array(output_mc['Baseline'][var]['y'])
pass_x_mc     = np.array(output_mc['Numerator'][var]['pass_x'])
pass_y_mc     = np.array(output_mc['Numerator'][var]['pass_y'])

print(f'\nData  - Baseline: {len(baseline_x_data):8d} events')
print(f'Data  - Numerator:{len(pass_x_data):8d} events  => eff = {len(pass_x_data)/len(baseline_x_data):.2%}')
print(f'MC    - Baseline: {len(baseline_x_mc):8d} events')
print(f'MC    - Numerator:{len(pass_x_mc):8d} events  => eff = {len(pass_x_mc)/len(baseline_x_mc):.2%}')

print('\n=== DATA pT distribution (baseline) ===')
for lo, hi in [(300,400),(400,500),(500,600),(600,700),(700,1000)]:
    cut = (baseline_y_data >= lo) & (baseline_y_data < hi)
    n   = np.sum(cut)
    cut_pass = (pass_y_data >= lo) & (pass_y_data < hi)
    np_ = np.sum(cut_pass)
    eff = np_ / n if n > 0 else 0
    print(f'  pT [{lo:4d},{hi:4d}): baseline={n:6d}, pass={np_:6d}, eff={eff:.2%}')

print('\n=== MC pT distribution (baseline) ===')
for lo, hi in [(300,400),(400,500),(500,600),(600,700),(700,1000)]:
    cut = (baseline_y_mc >= lo) & (baseline_y_mc < hi)
    n   = np.sum(cut)
    cut_pass = (pass_y_mc >= lo) & (pass_y_mc < hi)
    np_ = np.sum(cut_pass)
    eff = np_ / n if n > 0 else 0
    print(f'  pT [{lo:4d},{hi:4d}): baseline={n:6d}, pass={np_:6d}, eff={eff:.2%}')

print('\n=== PER-TRIGGER BREAKDOWN (DATA) ===')
for trig in output_data['PerTrigger']:
    px = np.array(output_data['PerTrigger'][trig][var]['pass_x'])
    py = np.array(output_data['PerTrigger'][trig][var]['pass_y'])
    eff = len(px) / len(baseline_x_data) if len(baseline_x_data) > 0 else 0
    print(f'  {trig[:55]:55s}: {len(px):6d} events ({eff:.2%})')

print('\n=== PER-TRIGGER BREAKDOWN (MC) ===')
for trig in output_mc['PerTrigger']:
    px = np.array(output_mc['PerTrigger'][trig][var]['pass_x'])
    py = np.array(output_mc['PerTrigger'][trig][var]['pass_y'])
    eff = len(px) / len(baseline_x_mc) if len(baseline_x_mc) > 0 else 0
    print(f'  {trig[:55]:55s}: {len(px):6d} events ({eff:.2%})')
