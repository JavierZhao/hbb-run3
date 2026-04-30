"""
Replot Inclusive pass scale factors for all years with:
  - One giant bin for pt > 600 GeV
  - Title showing efficiency in pt 300-700, msd 40-300
  - Produces 3-panel (Data Eff, MC Eff, SF) and 2-panel (SF, uncertainty) figures
"""
import warnings
warnings.filterwarnings('ignore')
import os
import sys
import types
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import cloudpickle

plt.style.use(hep.style.ROOT)

# --- Inject legacy coffea 0.7 accumulator classes so cloudpickle can deserialise ---
class _dict_accumulator(dict):
    def add(self, other): self.update(other)
    def identity(self): return type(self)()

class _list_accumulator(list):
    def add(self, other): self.extend(other)
    def identity(self): return type(self)()

class _value_accumulator:
    def __init__(self, val=0): self.value = val
    def add(self, other): self.value += other.value
    def identity(self): return type(self)(0)

class _set_accumulator(set):
    def add(self, other): self.update(other)
    def identity(self): return type(self)()

class _column_accumulator:
    def __init__(self, val=None): self.value = val
    def add(self, other): pass
    def identity(self): return type(self)()

_legacy_classes = {
    'dict_accumulator':   _dict_accumulator,
    'list_accumulator':   _list_accumulator,
    'value_accumulator':  _value_accumulator,
    'set_accumulator':    _set_accumulator,
    'column_accumulator': _column_accumulator,
}

# Create stub modules for both possible import paths coffea used
import coffea.processor as _cp
for _mod_name in ['coffea.processor.accumulator', 'coffea.processor._accumulator']:
    _mod = types.ModuleType(_mod_name)
    for _n, _cls in _legacy_classes.items():
        setattr(_mod, _n, _cls)
    sys.modules[_mod_name] = _mod

# Also inject directly onto coffea.processor
for _n, _cls in _legacy_classes.items():
    if not hasattr(_cp, _n):
        setattr(_cp, _n, _cls)

def _load_coffea(path):
    import lz4.frame
    with lz4.frame.open(path, 'rb') as f:
        return cloudpickle.load(f)

NOBACKUP = "/uscms/home/zzhao2/nobackup/hbb-run3/condor_sf"

YEARS = ["2022", "2022EE", "2023", "2023BPix"]

# msd: 20 GeV bins from 0 to 300 (same as original)
msd_bins = np.arange(0, 301, 20)
# pt: 50 GeV bins from 300 to 600, then ONE giant bin 600+
pt_bins_hist    = np.array([300, 350, 400, 450, 500, 550, 600, 2000])
pt_bins_display = np.array([300, 350, 400, 450, 500, 550, 600, 1000])

TITLE_REGION = r"$300<p_T<700$, $40<m_{SD}<300$ GeV"
VAR = "pt_vs_msd"


def process_year(year):
    base = os.path.join(NOBACKUP, f"year_{year}", "output")
    data_path = os.path.join(base, f"sf_2d_{year}_Inclusive_pass_data.coffea")
    mc_path   = os.path.join(base, f"sf_2d_{year}_Inclusive_pass_mc.coffea")
    if not os.path.exists(data_path) or not os.path.exists(mc_path):
        print(f"  Skipping {year}: coffea files not found")
        return

    out = _load_coffea(data_path)
    mc  = _load_coffea(mc_path)

    base_msd_d = np.array(out['Baseline'][VAR]['x'])
    base_pt_d  = np.array(out['Baseline'][VAR]['y'])
    base_msd_m = np.array(mc['Baseline'][VAR]['x'])
    base_pt_m  = np.array(mc['Baseline'][VAR]['y'])
    pass_msd_d = np.array(out['Numerator'][VAR]['pass_x'])
    pass_pt_d  = np.array(out['Numerator'][VAR]['pass_y'])
    pass_msd_m = np.array(mc['Numerator'][VAR]['pass_x'])
    pass_pt_m  = np.array(mc['Numerator'][VAR]['pass_y'])

    h_base_d, xedges, _ = np.histogram2d(base_msd_d, base_pt_d, bins=[msd_bins, pt_bins_hist])
    h_base_m, _,      _ = np.histogram2d(base_msd_m, base_pt_m, bins=[msd_bins, pt_bins_hist])
    h_pass_d, _,      _ = np.histogram2d(pass_msd_d, pass_pt_d, bins=[msd_bins, pt_bins_hist])
    h_pass_m, _,      _ = np.histogram2d(pass_msd_m, pass_pt_m, bins=[msd_bins, pt_bins_hist])

    yedges = pt_bins_display.astype(float)

    with np.errstate(divide='ignore', invalid='ignore'):
        eff_d = np.where(h_base_d > 0, h_pass_d / h_base_d, np.nan)
        eff_m = np.where(h_base_m > 0, h_pass_m / h_base_m, np.nan)

    valid = (h_base_d > 0) & (h_base_m > 0) & (eff_m > 0)
    sf    = np.where(valid, eff_d / eff_m, np.nan)

    with np.errstate(divide='ignore', invalid='ignore'):
        sf_unc = np.where(
            valid & (h_pass_d > 0) & (h_pass_m > 0),
            np.abs(sf) * np.sqrt(1.0 / h_pass_d + 1.0 / h_pass_m),
            np.nan
        )

    # Title efficiencies in signal region
    rm = lambda pt, msd: (pt >= 300) & (pt < 700) & (msd >= 40) & (msd < 300)
    def _reff(px, py, bx, by):
        n = np.sum(rm(py, px)); d = np.sum(rm(by, bx))
        return n / d if d > 0 else np.nan
    teff_d = _reff(pass_msd_d, pass_pt_d, base_msd_d, base_pt_d)
    teff_m = _reff(pass_msd_m, pass_pt_m, base_msd_m, base_pt_m)
    tsf    = teff_d / teff_m if teff_m > 0 else np.nan

    save_dir = os.path.join(NOBACKUP, f"year_{year}", "figures_sf_2d", year, "Inclusive", "pass")
    os.makedirs(save_dir, exist_ok=True)

    def _pc(ax, Z, cmap, vmin, vmax, title, cms_data, cbar_label=None):
        im = ax.pcolormesh(xedges, yedges, Z.T, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
        cb = plt.colorbar(im, ax=ax, pad=0.02)
        if cbar_label:
            cb.set_label(cbar_label, fontsize=11)
        cb.ax.tick_params(labelsize=10)
        ax.set_xlabel(r"Leading Jet $m_{SD}$ [GeV]", fontsize=11)
        ax.set_ylabel(r"Leading Jet $p_{T}$ [GeV]", fontsize=11)
        ax.set_title(title, fontsize=11, pad=4)
        ax.axhline(600, color='black', linestyle='--', linewidth=1.0, alpha=0.6)
        ax.text(xedges[-1] - 5, 620, r'$p_T>600$ bin', ha='right', va='bottom',
                fontsize=8, color='black', alpha=0.8)
        if cms_data:
            hep.cms.label(ax=ax, data=True, year=year, com="13.6", fontsize=10)
        else:
            hep.cms.label(ax=ax, data=False, year=year, com="13.6", label="Simulation", fontsize=10)
        return im

    # Figure 1: Data Eff | MC Eff | SF
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(24, 7))
    _pc(a1, eff_d, 'RdYlGn', 0, 1, f"Data Efficiency\n{TITLE_REGION}: {teff_d:.2%}", True)
    _pc(a2, eff_m, 'RdYlGn', 0, 1, f"MC Efficiency\n{TITLE_REGION}: {teff_m:.2%}",   False)
    _pc(a3, sf,   'RdBu_r', 0.8, 1.2, f"Scale Factor\n{TITLE_REGION}: {tsf:.3f}",   True,
        cbar_label='Scale Factor (Data/MC)')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "scale_factor_pt_vs_msd_2d.png"), dpi=200, bbox_inches='tight')
    print(f"  Saved: scale_factor_pt_vs_msd_2d.png")
    plt.close(fig)

    # Figure 2: SF | Uncertainty (reference style, viridis)
    fig2, (b1, b2) = plt.subplots(1, 2, figsize=(16, 7))
    _pc(b1, sf,     'viridis', 0, 1.2, f"Scale Factor\n{TITLE_REGION}: {tsf:.3f}", True,
        cbar_label='Trigger efficiency scale factor')
    _pc(b2, sf_unc, 'viridis', 0, 0.1, "Scale Factor Statistical Uncertainty",     True,
        cbar_label='Trigger efficiency scale factor uncertainty')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "scale_factor_pt_vs_msd_2d_unc.png"), dpi=200, bbox_inches='tight')
    print(f"  Saved: scale_factor_pt_vs_msd_2d_unc.png")
    plt.close(fig2)


for year in YEARS:
    print(f"\nProcessing {year}...")
    process_year(year)
