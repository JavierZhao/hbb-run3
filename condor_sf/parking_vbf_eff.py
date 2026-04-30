"""
ParkingVBF Trigger Efficiency Study
=====================================
Measures the efficiency of the VBF parking trigger:
    HLT_VBF_DiPFJet125_45_Mjj720_Detajj3p0

Method: muon tag-and-probe on ParkingVBF NanoAOD data
  - Reference (denominator): events passing single-muon triggers
  - Signal (numerator):      reference AND VBF trigger fires

Variables:
  - mjj vs delta_eta (primary 2D map)
  - leading VBF jet pT vs subleading VBF jet pT
  - leading VBF jet pT (1D)

Output:
  figures_parking_vbf/{year}/vbf_trig_eff_mjj_vs_deta.png
  output/{year}/vbf_parking_trig_eff_{year}.csv

Usage:
    ./shell coffeateam/coffea-dask:latest-py3.9 python3 parking_vbf_eff.py --year 2023
    ./shell coffeateam/coffea-dask:latest-py3.9 python3 parking_vbf_eff.py --year 2023BPix
    ./shell coffeateam/coffea-dask:latest-py3.9 python3 parking_vbf_eff.py --year 2023 --test

Note: ParkingVBF NanoAOD only exists for Run2023 (not 2022/2022EE).
"""

import warnings
warnings.filterwarnings('ignore')
import os
import json
import csv

import awkward as ak
import numpy as np
import uproot

import hist
from coffea import processor, util
from coffea.nanoevents import NanoAODSchema
from coffea.processor import dict_accumulator, list_accumulator

import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
plt.style.use(hep.style.ROOT)

# ─── Trigger configuration ────────────────────────────────────────────────────

# The signal trigger we want to measure efficiency for
VBF_SIGNAL_TRIGGER = 'VBF_DiPFJet125_45_Mjj720_Detajj3p0'

# Muon reference triggers (unbiased w.r.t. jet measurement)
MUON_TRIGGERS = {
    '2023':     ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023BPix': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
}

# ParkingVBF dataset keys to use from the infiles JSON
# These are the highest-statistics versions of Run2023C and Run2023D
PARKING_VBF_KEYS = {
    '2023':     ['ParkingVBF_Run2023C-v3', 'ParkingVBF_Run2023C-v4'],
    '2023BPix': ['ParkingVBF_Run2023D-v1', 'ParkingVBF_Run2023D-v2'],
}

# ─── 2D efficiency variable definitions ───────────────────────────────────────

def compute_vbf_dijet(events):
    """
    Compute VBF dijet kinematics: the two highest-pT AK4 jets
    outside any leading FatJet (dR > 0.8).

    Returns: (mjj [GeV], delta_eta, jet1_pt, jet2_pt)
    """
    fatjets = events.FatJet
    candidatejet = ak.firsts(
        fatjets[(fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight][:, :2]
    )

    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 5.0) & jets.isTight][:, :4]

    ak4_outside = jets[ak.fill_none(jets.delta_r(candidatejet), 999) > 0.8]

    jet1 = ak.firsts(ak4_outside[:, 0:1])
    jet2 = ak.firsts(ak4_outside[:, 1:2])

    mjj  = ak.fill_none((jet1 + jet2).mass,  np.nan)
    deta = ak.fill_none(abs(jet1.eta - jet2.eta), np.nan)
    pt1  = ak.fill_none(jet1.pt, np.nan)
    pt2  = ak.fill_none(jet2.pt, np.nan)

    return mjj, deta, pt1, pt2


TRIG_VARS = {
    'mjj_vs_deta': {
        'label_x': r"Dijet $m_{jj}$ [GeV]",
        'label_y': r"Dijet $|\Delta\eta|$",
        'axis_x': hist.axis.Regular(bins=20, start=0,    stop=2000, name="mjj",  label=r"$m_{jj}$ [GeV]"),
        'axis_y': hist.axis.Regular(bins=16, start=0,    stop=8,    name="deta", label=r"$|\Delta\eta|$"),
        'proc_x': lambda events: compute_vbf_dijet(events)[0],
        'proc_y': lambda events: compute_vbf_dijet(events)[1],
    },
    'pt1_vs_pt2': {
        'label_x': r"Leading VBF jet $p_{T}$ [GeV]",
        'label_y': r"Subleading VBF jet $p_{T}$ [GeV]",
        'axis_x': hist.axis.Regular(bins=20, start=30,  stop=500,  name="pt1",  label=r"Leading VBF jet $p_T$ [GeV]"),
        'axis_y': hist.axis.Regular(bins=20, start=30,  stop=300,  name="pt2",  label=r"Subleading VBF jet $p_T$ [GeV]"),
        'proc_x': lambda events: compute_vbf_dijet(events)[2],
        'proc_y': lambda events: compute_vbf_dijet(events)[3],
    },
}


# ─── Baseline selection ────────────────────────────────────────────────────────

def create_baseline_mask(events):
    """
    Select events suitable for VBF trigger efficiency measurement.

    Requirements:
      - Fired a single-muon reference trigger
      - ≥1 loose muon (pT > 25 GeV, |η| < 2.4, isolated, dR > 0.8 from leading FatJet)
      - ≥1 FatJet (pT > 200 GeV, |η| < 2.5, tight)
      - ≥2 AK4 jets outside leading FatJet with VBF-like topology (loose)
    """
    n = len(events)
    mask = np.ones(n, dtype=bool)

    # Require ≥1 FatJet (loose, just to anchor the muon dR cut)
    fatjets = events.FatJet
    has_fatjet = ak.sum(
        (fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight, axis=1
    ) >= 1
    mask &= ak.to_numpy(has_fatjet)

    # Require ≥1 loose muon far from leading FatJet
    muons = events.Muon
    loose_mu = muons[
        (muons.pt > 25) & (abs(muons.eta) < 2.4) & muons.looseId & (muons.pfRelIso04_all < 0.25)
    ]
    has_loose_mu = ak.num(loose_mu) >= 1

    leading_fj = ak.firsts(ak.pad_none(events.FatJet, 1, clip=True))
    dphi = np.abs(((loose_mu.phi - leading_fj.phi + np.pi) % (2 * np.pi)) - np.pi)
    deta_mu = loose_mu.eta - leading_fj.eta
    dr = np.sqrt(dphi**2 + deta_mu**2)
    far_enough = ak.fill_none(ak.any(dr > 0.8, axis=1), True)

    mask &= ak.to_numpy(has_loose_mu & far_enough)

    # Require ≥2 AK4 jets outside FatJet (needed to compute mjj/deta)
    candidatejet = ak.firsts(
        fatjets[(fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight][:, :2]
    )
    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 5.0) & jets.isTight][:, :4]
    ak4_outside = jets[ak.fill_none(jets.delta_r(candidatejet), 999) > 0.8]
    has_two_vbf_jets = ak.num(ak4_outside) >= 2
    mask &= ak.to_numpy(has_two_vbf_jets)

    return mask


# ─── Coffea Processor ─────────────────────────────────────────────────────────

class ParkingVBFEffProcessor(processor.ProcessorABC):
    """
    Measure VBF trigger efficiency on ParkingVBF data using muon tag-and-probe.
    """

    def __init__(self, year, trig_vars=None):
        self.year = str(year)
        self.trig_vars = trig_vars or TRIG_VARS
        self.muon_triggers = MUON_TRIGGERS.get(self.year, [])

        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                vname: dict_accumulator({'x': list_accumulator(), 'y': list_accumulator()})
                for vname in self.trig_vars
            }),
            'VBFPass': dict_accumulator({
                vname: dict_accumulator({'x': list_accumulator(), 'y': list_accumulator()})
                for vname in self.trig_vars
            }),
        })

    def process(self, events):
        out = self.output

        # ── Step 1: require muon reference trigger ──────────────────────────
        fired_ref = ak.zeros_like(events.event, dtype=bool)
        for trig in self.muon_triggers:
            field = trig[4:] if trig.startswith('HLT_') else trig
            if field in events.HLT.fields:
                fired_ref = fired_ref | events.HLT[field]
        events = events[fired_ref]

        if len(events) == 0:
            return out

        # ── Step 2: baseline selection ──────────────────────────────────────
        baseline = create_baseline_mask(events)

        # ── Step 3: compute kinematic variables ─────────────────────────────
        variables = {}
        for vname, vinfo in self.trig_vars.items():
            vx = vinfo['proc_x'](events)
            vy = vinfo['proc_y'](events)
            variables[vname] = {'x': vx, 'y': vy}

        valid_mask = np.ones(len(events), dtype=bool)
        for vdata in variables.values():
            valid_mask &= ~np.isnan(ak.to_numpy(vdata['x']))
            valid_mask &= ~np.isnan(ak.to_numpy(vdata['y']))

        sel_baseline = baseline & valid_mask

        # ── Step 4: check VBF signal trigger ────────────────────────────────
        vbf_fired = ak.zeros_like(events.event, dtype=bool)
        if VBF_SIGNAL_TRIGGER in events.HLT.fields:
            vbf_fired = events.HLT[VBF_SIGNAL_TRIGGER]
        else:
            print(f"WARNING: {VBF_SIGNAL_TRIGGER} not found in HLT fields!")

        sel_pass = baseline & vbf_fired & valid_mask

        # ── Step 5: fill histograms ──────────────────────────────────────────
        for vname in self.trig_vars:
            vx = variables[vname]['x']
            vy = variables[vname]['y']

            out['Baseline'][vname]['x'].extend(ak.to_numpy(vx[sel_baseline]).tolist())
            out['Baseline'][vname]['y'].extend(ak.to_numpy(vy[sel_baseline]).tolist())
            out['VBFPass'][vname]['x'].extend(ak.to_numpy(vx[sel_pass]).tolist())
            out['VBFPass'][vname]['y'].extend(ak.to_numpy(vy[sel_pass]).tolist())

        return out

    def postprocess(self, accumulator):
        for region in ('Baseline', 'VBFPass'):
            for vname in accumulator[region]:
                for coord in ('x', 'y'):
                    accumulator[region][vname][coord] = np.array(
                        accumulator[region][vname][coord]
                    )
        return accumulator


# ─── Plotting ─────────────────────────────────────────────────────────────────

def plot_2d_efficiency(output, trig_vars, save_dir, year):
    os.makedirs(save_dir, exist_ok=True)

    for vname, vinfo in trig_vars.items():
        bins_x = vinfo['axis_x'].edges
        bins_y = vinfo['axis_y'].edges

        base_x = np.array(output['Baseline'][vname]['x'])
        base_y = np.array(output['Baseline'][vname]['y'])
        pass_x = np.array(output['VBFPass'][vname]['x'])
        pass_y = np.array(output['VBFPass'][vname]['y'])

        h_base, xedges, yedges = np.histogram2d(base_x, base_y, bins=[bins_x, bins_y])
        h_pass, _, _           = np.histogram2d(pass_x, pass_y, bins=[bins_x, bins_y])

        with np.errstate(divide='ignore', invalid='ignore'):
            eff = np.where(h_base > 0, h_pass / h_base, 0.0)

        total_base = int(np.sum(h_base))
        total_pass = int(np.sum(h_pass))
        total_eff  = total_pass / total_base if total_base > 0 else 0.0

        fig, ax = plt.subplots(figsize=(9, 7))
        im = ax.pcolormesh(xedges, yedges, eff.T, cmap='RdYlGn', vmin=0, vmax=1, shading='flat')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Trigger Efficiency', fontsize=12)
        ax.set_xlabel(vinfo['label_x'], fontsize=12)
        ax.set_ylabel(vinfo['label_y'], fontsize=12)
        ax.set_title(
            f"VBF Trigger Efficiency\n"
            f"HLT_{VBF_SIGNAL_TRIGGER}\n"
            f"Overall: {total_eff:.2%}  ({total_pass}/{total_base} events)",
            fontsize=11
        )
        hep.cms.label(ax=ax, data=True, year=year, com="13.6", fontsize=10)

        # Overlay efficiency values in each bin
        for ix in range(len(bins_x) - 1):
            for iy in range(len(bins_y) - 1):
                if h_base[ix, iy] > 0:
                    ax.text(
                        0.5 * (xedges[ix] + xedges[ix+1]),
                        0.5 * (yedges[iy] + yedges[iy+1]),
                        f"{eff[ix, iy]:.2f}",
                        ha='center', va='center', fontsize=6, color='black'
                    )

        plt.tight_layout()
        fname = os.path.join(save_dir, f"vbf_trig_eff_{vname}.png")
        plt.savefig(fname, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {fname}")


def generate_efficiency_csv(output, trig_vars, save_dir, year,
                            mjj_min=720, deta_min=3.0):
    """
    Report overall VBF trigger efficiency in the signal region
    (mjj > mjj_min, |Δη| > deta_min).
    """
    os.makedirs(save_dir, exist_ok=True)

    vname = 'mjj_vs_deta'
    if vname not in trig_vars:
        return

    base_x = np.array(output['Baseline'][vname]['x'])  # mjj
    base_y = np.array(output['Baseline'][vname]['y'])  # deta
    pass_x = np.array(output['VBFPass'][vname]['x'])
    pass_y = np.array(output['VBFPass'][vname]['y'])

    # Signal region cut
    sig_cut_base = (base_x > mjj_min) & (base_y > deta_min)
    sig_cut_pass = (pass_x > mjj_min) & (pass_y > deta_min)

    rows = [['Region', 'Trigger', 'Efficiency', 'Pass', 'Total']]

    for label, bc, pc in [
        ('Inclusive',      np.ones(len(base_x), dtype=bool), np.ones(len(pass_x), dtype=bool)),
        (f'mjj>{mjj_min},deta>{deta_min}', sig_cut_base, sig_cut_pass),
    ]:
        n_base = int(np.sum(bc))
        n_pass = int(np.sum(pc))
        eff = n_pass / n_base if n_base > 0 else 0.0
        rows.append([label, f"HLT_{VBF_SIGNAL_TRIGGER}", f"{eff:.4f}", str(n_pass), str(n_base)])
        print(f"  {label:35s} | Eff: {eff:.2%} ({n_pass}/{n_base})")

    csv_path = os.path.join(save_dir, f"vbf_parking_trig_eff_{year}.csv")
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f"Saved CSV: {csv_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="ParkingVBF Trigger Efficiency Study")
    parser.add_argument('--year', required=False,
                        choices=['2023', '2023BPix'], default='2023',
                        help="Year/era: 2023 (Run2023C) or 2023BPix (Run2023D)")
    parser.add_argument('--test', action='store_true',
                        help="Test mode: use only first 5 files per dataset")
    parser.add_argument('--streams', type=str, default='0,1,2,3,4,5,6,7,8,9',
                        help="Comma-separated ParkingVBF stream numbers to use (default: 0-9)")
    args = parser.parse_args()

    year = args.year
    test_mode = args.test
    stream_nums = [int(s) for s in args.streams.split(',')]

    os.makedirs("figures_parking_vbf", exist_ok=True)
    os.makedirs("output", exist_ok=True)

    print(f"\n{'='*60}")
    print(f"ParkingVBF Trigger Efficiency Study — {year}")
    print(f"Signal trigger: HLT_{VBF_SIGNAL_TRIGGER}")
    print(f"Streams: ParkingVBF{stream_nums}")
    if test_mode:
        print("TEST MODE: 5 files per dataset")
    print(f"{'='*60}\n")

    # ── Load infiles ───────────────────────────────────────────────────────
    infile_path = f"infiles/{year}/{year}_ParkingVBFData.json"
    if not os.path.exists(infile_path):
        print(f"ERROR: {infile_path} not found.")
        print("Run make_parking_vbf_infiles.py first (needs dasgoclient + grid proxy).")
        import sys; sys.exit(1)

    with open(infile_path) as f:
        parking_data = json.load(f)

    # Collect files across requested streams and versions
    fileset = {}
    stream_prefixes = [f"ParkingVBF{i}_" for i in stream_nums]

    for key, files in parking_data.items():
        if not any(key.startswith(pfx) for pfx in stream_prefixes):
            continue
        if not files:
            continue
        if test_mode:
            files = files[:5]
        # Merge all matching keys into one fileset entry
        if 'ParkingVBF' not in fileset:
            fileset['ParkingVBF'] = []
        fileset['ParkingVBF'].extend(files)

    if not fileset:
        print("ERROR: No files loaded. Check infiles JSON and --streams argument.")
        import sys; sys.exit(1)

    print(f"Loaded {len(fileset['ParkingVBF'])} total files")

    # ── Run processor ──────────────────────────────────────────────────────
    iterative_run = processor.Runner(
        executor=processor.FuturesExecutor(compression=None, workers=2),
        schema=NanoAODSchema,
        skipbadfiles=True,
        savemetrics=True,
    )

    print("\nRunning processor...")
    out, metrics = iterative_run(
        fileset,
        treename="Events",
        processor_instance=ParkingVBFEffProcessor(year, TRIG_VARS),
    )

    # ── Save and plot ──────────────────────────────────────────────────────
    coffea_path = os.path.join("output", f"parking_vbf_eff_{year}.coffea")
    util.save(out, coffea_path)
    print(f"Saved output: {coffea_path}")

    FIGURES_DIR = os.path.join("figures_parking_vbf", year)
    plot_2d_efficiency(out, TRIG_VARS, save_dir=FIGURES_DIR, year=year)

    CSV_DIR = os.path.join("output", year)
    os.makedirs(CSV_DIR, exist_ok=True)
    generate_efficiency_csv(out, TRIG_VARS, save_dir=CSV_DIR, year=year)

    print(f"\nDone. Figures: {FIGURES_DIR}/")
