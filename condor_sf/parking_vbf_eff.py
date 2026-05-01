"""
VBF Parking Trigger Efficiency / Scale-Factor Study
====================================================
Measures the efficiency of HLT_VBF_DiPFJet125_45_Mjj720_Detajj3p0 (and friends)
in two complementary ways, selected via --method:

  1. mu_tag_parkingvbf  (legacy, data-only)
       Sample:      ParkingVBF NanoAOD
       Reference:   single-muon triggers (HLT_Mu50 etc.)
       Numerator:   reference & VBF parking trigger
       Output:      data efficiency map only

  2. zmumu_vbf          (new, AN-style with MC for SF)
       Sample:      Muon PD (data) + DYJetsToLL (MC, key "DYJets" in <year>_Zjets.json)
       Reference:   HLT_IsoMu24 (unbiased w.r.t. the hadronic activity)
       Selection:   Tomas's DY+0L CR (VBF_Studies_Tomas_NewFW_090525.pdf, slide 5),
                    adapted as  Z(mumu) + 2 VBF-like AK4 jets:
                       - exactly 2 tight isolated muons, pT > 30 / 15, |eta| < 2.1, OS
                       - 60 < m(mumu) < 120
                       - electron / extra-muon vetoes
                       - b-tag veto (medium DeepFlavB)
                       - >= 2 AK4 jets pT > 30, |eta| < 4.7, opposite hemispheres,
                         |dEta(j1,j2)| > 3.8, m(jj) > 500
       Numerator:   reference & VBF parking trigger
       Output:      data eff, MC eff, and SF = data/MC maps with stat. uncertainty.

Variables (both methods):
  - mjj vs |Delta eta|   (primary 2D)
  - leading VBF jet pT vs subleading VBF jet pT

Outputs:
  figures_parking_vbf/<year>/<method>/...
  output/<year>/<method>/...

Usage:
  python3 parking_vbf_eff.py --year 2023 --method zmumu_vbf
  python3 parking_vbf_eff.py --year 2023 --method mu_tag_parkingvbf

Note: ParkingVBF NanoAOD only exists for Run2023, but zmumu_vbf works for any
year that has Muon PD data + DYJetsToLL MC + the VBF parking trigger in HLT.
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

VBF_SIGNAL_TRIGGER = 'VBF_DiPFJet125_45_Mjj720_Detajj3p0'

# Reference triggers per method (HLT_ prefix optional in the dict, stripped at use)
MU_TAG_REF_TRIGGERS = {
    '2022':     ["HLT_Mu50", "HLT_IsoMu24"],
    '2022EE':   ["HLT_Mu50", "HLT_IsoMu24"],
    '2023':     ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023BPix': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
}

# For zmumu_vbf we want a single, unbiased single-muon trigger. IsoMu24 is the
# Run-3 equivalent of Tomas's IsoMu27 from his Run-2 study.
ZMUMU_REF_TRIGGER = "HLT_IsoMu24"

# Data-fileset keys per period (legacy method uses ParkingVBF, new uses Muon PD)
PARKING_VBF_KEYS = {
    '2023':     ['ParkingVBF_Run2023C-v3', 'ParkingVBF_Run2023C-v4'],
    '2023BPix': ['ParkingVBF_Run2023D-v1', 'ParkingVBF_Run2023D-v2'],
}

MUON_DATA_KEYS = {
    '2022':     ["Muon_Run2022C"],
    '2022EE':   ["Muon_Run2022E"],
    '2023':     [
        "Muon0_Run2023C-v1", "Muon0_Run2023C-v2", "Muon0_Run2023C-v3", "Muon0_Run2023C-v4",
        "Muon1_Run2023C-v1", "Muon1_Run2023C-v2", "Muon1_Run2023C-v3", "Muon1_Run2023C-v4",
    ],
    '2023BPix': ["Muon0_Run2023D-v1", "Muon0_Run2023D-v2"],
}

# DYJetsToLL MC key in <year>_Zjets.json
DY_MC_KEY = "DYJets"

# ─── 2D variable definitions ──────────────────────────────────────────────────

def compute_vbf_dijet_outside_fatjet(events):
    """
    mu_tag_parkingvbf method: leading 2 AK4 jets outside any leading FatJet (dR > 0.8).
    """
    fatjets = events.FatJet
    candidatejet = ak.firsts(
        fatjets[(fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight][:, :2]
    )

    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 5.0) & jets.isTight][:, :4]
    ak4_outside = jets[ak.fill_none(jets.delta_r(candidatejet), 999) > 0.8]

    j1 = ak.firsts(ak4_outside[:, 0:1])
    j2 = ak.firsts(ak4_outside[:, 1:2])

    mjj  = ak.fill_none((j1 + j2).mass, np.nan)
    deta = ak.fill_none(abs(j1.eta - j2.eta), np.nan)
    pt1  = ak.fill_none(j1.pt, np.nan)
    pt2  = ak.fill_none(j2.pt, np.nan)
    return mjj, deta, pt1, pt2


def compute_zmumu_vbf_dijet(events):
    """
    zmumu_vbf method: leading 2 AK4 jets in |eta| < 4.7 (no FatJet veto needed).
    """
    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 4.7) & jets.isTight]

    j1 = ak.firsts(jets[:, 0:1])
    j2 = ak.firsts(jets[:, 1:2])

    mjj  = ak.fill_none((j1 + j2).mass, np.nan)
    deta = ak.fill_none(abs(j1.eta - j2.eta), np.nan)
    pt1  = ak.fill_none(j1.pt, np.nan)
    pt2  = ak.fill_none(j2.pt, np.nan)
    return mjj, deta, pt1, pt2


def trig_vars_for_method(method):
    """Pick the right VBF-pair function depending on the selection method."""
    proc = compute_zmumu_vbf_dijet if method == 'zmumu_vbf' else compute_vbf_dijet_outside_fatjet
    return {
        'mjj_vs_deta': {
            'label_x': r"Dijet $m_{jj}$ [GeV]",
            'label_y': r"Dijet $|\Delta\eta|$",
            'axis_x': hist.axis.Regular(bins=20, start=0,   stop=2000, name="mjj",  label=r"$m_{jj}$ [GeV]"),
            'axis_y': hist.axis.Regular(bins=16, start=0,   stop=8,    name="deta", label=r"$|\Delta\eta|$"),
            'proc_x': lambda events, _p=proc: _p(events)[0],
            'proc_y': lambda events, _p=proc: _p(events)[1],
        },
        'pt1_vs_pt2': {
            'label_x': r"Leading VBF jet $p_{T}$ [GeV]",
            'label_y': r"Subleading VBF jet $p_{T}$ [GeV]",
            'axis_x': hist.axis.Regular(bins=20, start=30, stop=500, name="pt1", label=r"Leading VBF jet $p_T$ [GeV]"),
            'axis_y': hist.axis.Regular(bins=20, start=30, stop=300, name="pt2", label=r"Subleading VBF jet $p_T$ [GeV]"),
            'proc_x': lambda events, _p=proc: _p(events)[2],
            'proc_y': lambda events, _p=proc: _p(events)[3],
        },
    }


# ─── Selections ───────────────────────────────────────────────────────────────

def create_mu_tag_baseline_mask(events):
    """Legacy ParkingVBF muon-tag selection (kept for back-compat)."""
    n = len(events)
    mask = np.ones(n, dtype=bool)

    fatjets = events.FatJet
    has_fatjet = ak.sum(
        (fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight, axis=1
    ) >= 1
    mask &= ak.to_numpy(has_fatjet)

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

    candidatejet = ak.firsts(
        fatjets[(fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight][:, :2]
    )
    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 5.0) & jets.isTight][:, :4]
    ak4_outside = jets[ak.fill_none(jets.delta_r(candidatejet), 999) > 0.8]
    has_two_vbf_jets = ak.num(ak4_outside) >= 2
    mask &= ak.to_numpy(has_two_vbf_jets)

    return mask


def create_zmumu_vbf_mask(events):
    """
    Z(mumu) + 2 VBF-like jets selection adapted from Tomas's DY+0L CR (slide 5):
      - exactly 2 tight isolated muons, leading pT > 30 / sub-leading > 15, |eta| < 2.1, opposite sign
      - 60 < m(mumu) < 120 GeV
      - electron veto (loose, pT > 5)
      - b-tag veto (medium DeepFlavB on AK4 with pT > 30, |eta| < 2.4)
      - >= 2 AK4 jets (PUPPI tight, pT > 30, |eta| < 4.7) with
            opposite hemispheres, |dEta(j1,j2)| > 3.8, m(jj) > 500 GeV.

    Note: tau veto is omitted (NanoAOD tau collection isn't reliable for this purpose
    in muon-streamed data). Photon veto is also omitted.
    """
    n = len(events)
    mask = np.ones(n, dtype=bool)

    # Tight isolated muons
    muons = events.Muon
    tight_mu = muons[
        (muons.pt > 15.0)
        & (abs(muons.eta) < 2.1)
        & muons.tightId
        & (muons.pfRelIso04_all < 0.15)
    ]
    mask &= ak.to_numpy(ak.num(tight_mu) == 2)

    mu1 = ak.pad_none(tight_mu, 2, clip=True)[:, 0]
    mu2 = ak.pad_none(tight_mu, 2, clip=True)[:, 1]

    mu1_pt_ok  = ak.fill_none(mu1.pt > 30.0, False)
    opp_sign   = ak.fill_none(mu1.charge * mu2.charge < 0, False)
    mll        = ak.fill_none((mu1 + mu2).mass, 0.0)
    z_window   = (mll > 60.0) & (mll < 120.0)
    mask &= ak.to_numpy(mu1_pt_ok & opp_sign & z_window)

    # Electron veto
    electrons = events.Electron
    veto_e = electrons[
        (electrons.pt > 5.0)
        & (abs(electrons.eta) < 2.5)
    ]
    mask &= ak.to_numpy(ak.num(veto_e) == 0)

    # b-tag veto (medium DeepFlavB ~ 0.30; conservative, dataset-dependent WP)
    jets_all = events.Jet
    btag_field = "btagDeepFlavB" if "btagDeepFlavB" in jets_all.fields else None
    if btag_field is not None:
        btag_cand = jets_all[
            (jets_all.pt > 30.0)
            & (abs(jets_all.eta) < 2.4)
            & (jets_all[btag_field] > 0.30)
        ]
        mask &= ak.to_numpy(ak.num(btag_cand) == 0)

    # >= 2 AK4 jets in |eta| < 4.7, tight ID
    vbf_jets = jets_all[
        (jets_all.pt > 30.0) & (abs(jets_all.eta) < 4.7) & jets_all.isTight
    ]
    mask &= ak.to_numpy(ak.num(vbf_jets) >= 2)

    j1 = ak.pad_none(vbf_jets, 2, clip=True)[:, 0]
    j2 = ak.pad_none(vbf_jets, 2, clip=True)[:, 1]

    opp_hemi = ak.fill_none(j1.eta * j2.eta < 0, False)
    deta     = ak.fill_none(abs(j1.eta - j2.eta), 0.0)
    mjj      = ak.fill_none((j1 + j2).mass, 0.0)
    vbf_topo = opp_hemi & (deta > 3.8) & (mjj > 500.0)
    mask &= ak.to_numpy(vbf_topo)

    return mask


# ─── Coffea Processor ─────────────────────────────────────────────────────────

class VBFTrigEffProcessor(processor.ProcessorABC):
    """
    Generic VBF-trigger efficiency processor that supports both methods.
    """

    def __init__(self, year, method, trig_vars):
        self.year = str(year)
        self.method = method
        self.trig_vars = trig_vars

        if method == 'zmumu_vbf':
            self._ref_triggers = [ZMUMU_REF_TRIGGER]
            self._selection_fn = create_zmumu_vbf_mask
        else:  # mu_tag_parkingvbf
            self._ref_triggers = MU_TAG_REF_TRIGGERS.get(self.year, [])
            self._selection_fn = create_mu_tag_baseline_mask

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

        # 1. require reference trigger
        fired_ref = ak.zeros_like(events.event, dtype=bool)
        for trig in self._ref_triggers:
            field = trig[4:] if trig.startswith('HLT_') else trig
            if field in events.HLT.fields:
                fired_ref = fired_ref | events.HLT[field]
        events = events[fired_ref]
        if len(events) == 0:
            return out

        # 2. baseline selection (method-specific)
        baseline = self._selection_fn(events)

        # 3. compute kinematic variables
        variables = {}
        for vname, vinfo in self.trig_vars.items():
            variables[vname] = {
                'x': vinfo['proc_x'](events),
                'y': vinfo['proc_y'](events),
            }
        valid_mask = np.ones(len(events), dtype=bool)
        for vdata in variables.values():
            valid_mask &= ~np.isnan(ak.to_numpy(vdata['x']))
            valid_mask &= ~np.isnan(ak.to_numpy(vdata['y']))

        sel_baseline = baseline & valid_mask

        # 4. VBF signal trigger
        vbf_fired = ak.zeros_like(events.event, dtype=bool)
        if VBF_SIGNAL_TRIGGER in events.HLT.fields:
            vbf_fired = events.HLT[VBF_SIGNAL_TRIGGER]
        sel_pass = baseline & vbf_fired & valid_mask

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

def _binomial_eff_err(k, n):
    k = np.asarray(k, dtype=float); n = np.asarray(n, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.where(n > 0, k / n, np.nan)
        err = np.where(n > 0, np.sqrt(np.clip(eff * (1.0 - eff), 0.0, None) / n), np.nan)
    return eff, err


def plot_2d_efficiency(output, trig_vars, save_dir, year, label_suffix=""):
    """Single-panel 2D efficiency map (used for both data-only and MC-only)."""
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

        total_eff = (np.sum(h_pass) / np.sum(h_base)) if np.sum(h_base) > 0 else 0.0

        fig, ax = plt.subplots(figsize=(9, 7))
        im = ax.pcolormesh(xedges, yedges, eff.T, cmap='RdYlGn', vmin=0, vmax=1, shading='flat')
        plt.colorbar(im, ax=ax, label='Trigger Efficiency')
        ax.set_xlabel(vinfo['label_x']); ax.set_ylabel(vinfo['label_y'])
        ax.set_title(
            f"VBF Trigger Efficiency {label_suffix}\n"
            f"HLT_{VBF_SIGNAL_TRIGGER}\n"
            f"Overall: {total_eff:.2%}"
        )
        hep.cms.label(ax=ax, data=(label_suffix.lower().startswith('data') or label_suffix == ''),
                      year=year, com="13.6", fontsize=10)
        plt.tight_layout()
        fname = os.path.join(save_dir, f"vbf_trig_eff_{vname}{('_'+label_suffix) if label_suffix else ''}.png")
        plt.savefig(fname, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {fname}")


def plot_2d_data_mc_sf(output_data, output_mc, trig_vars, save_dir, year):
    """3-panel data-eff / MC-eff / SF (data/MC) map for the zmumu_vbf method."""
    os.makedirs(save_dir, exist_ok=True)
    for vname, vinfo in trig_vars.items():
        bins_x = vinfo['axis_x'].edges
        bins_y = vinfo['axis_y'].edges

        def _hist2d(out, sub):
            x = np.array(out[sub][vname]['x'])
            y = np.array(out[sub][vname]['y'])
            h, xe, ye = np.histogram2d(x, y, bins=[bins_x, bins_y])
            return h, xe, ye

        h_base_d, xe, ye = _hist2d(output_data, 'Baseline')
        h_pass_d, _, _   = _hist2d(output_data, 'VBFPass')
        h_base_m, _, _   = _hist2d(output_mc,   'Baseline')
        h_pass_m, _, _   = _hist2d(output_mc,   'VBFPass')

        eff_d, err_d = _binomial_eff_err(h_pass_d, h_base_d)
        eff_m, err_m = _binomial_eff_err(h_pass_m, h_base_m)

        with np.errstate(divide='ignore', invalid='ignore'):
            sf = np.where((eff_m > 0) & np.isfinite(eff_d), eff_d / eff_m, np.nan)
            rel = np.where(h_pass_d > 0, (1 - eff_d) / h_pass_d, np.nan) \
                + np.where(h_pass_m > 0, (1 - eff_m) / h_pass_m, np.nan)
            sf_err = np.where(np.isfinite(sf) & np.isfinite(rel),
                              np.abs(sf) * np.sqrt(np.clip(rel, 0.0, None)), np.nan)

        tot_d = (np.sum(h_pass_d) / np.sum(h_base_d)) if np.sum(h_base_d) > 0 else 0.0
        tot_m = (np.sum(h_pass_m) / np.sum(h_base_m)) if np.sum(h_base_m) > 0 else 0.0
        tot_sf = (tot_d / tot_m) if tot_m > 0 else 1.0

        fig, axes = plt.subplots(1, 3, figsize=(24, 7))
        for ax, mat, cmap, vmin, vmax, label, title, is_data in [
            (axes[0], eff_d, 'RdYlGn', 0, 1,   'Efficiency',          f"Data Eff\nOverall: {tot_d:.2%}", True),
            (axes[1], eff_m, 'RdYlGn', 0, 1,   'Efficiency',          f"MC Eff\nOverall: {tot_m:.2%}", False),
            (axes[2], sf,    'RdBu_r', 0.5, 1.5, 'Scale Factor (D/MC)', f"Scale Factor\nOverall: {tot_sf:.3f}", True),
        ]:
            im = ax.pcolormesh(xe, ye, mat.T, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
            plt.colorbar(im, ax=ax, label=label)
            ax.set_xlabel(vinfo['label_x']); ax.set_ylabel(vinfo['label_y'])
            ax.set_title(title)
            hep.cms.label(ax=ax, data=is_data, year=year, com="13.6",
                          label=("Simulation" if not is_data else None), fontsize=9)
        plt.tight_layout()
        fname = os.path.join(save_dir, f"vbf_trig_sf_{vname}.png")
        plt.savefig(fname, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {fname}")

        # Companion 2-panel: SF + SF stat. uncertainty
        fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        im1 = ax1.pcolormesh(xe, ye, sf.T,     cmap='viridis', vmin=0, vmax=1.5, shading='flat')
        im2 = ax2.pcolormesh(xe, ye, sf_err.T, cmap='viridis', vmin=0, vmax=0.3, shading='flat')
        plt.colorbar(im1, ax=ax1, label='Scale Factor')
        plt.colorbar(im2, ax=ax2, label='SF stat. uncertainty')
        for ax, t in [(ax1, "SF"), (ax2, "SF stat. uncertainty")]:
            ax.set_xlabel(vinfo['label_x']); ax.set_ylabel(vinfo['label_y']); ax.set_title(t)
            hep.cms.label(ax=ax, data=True, year=year, com="13.6", fontsize=9)
        plt.tight_layout()
        fname2 = os.path.join(save_dir, f"vbf_trig_sf_{vname}_unc.png")
        plt.savefig(fname2, dpi=200, bbox_inches='tight')
        plt.close(fig2)
        print(f"Saved: {fname2}")


def generate_efficiency_csv(output, trig_vars, save_dir, year, suffix="",
                            mjj_min=720, deta_min=3.0):
    """Inclusive + signal-region (mjj > mjj_min, dEta > deta_min) efficiency table."""
    os.makedirs(save_dir, exist_ok=True)
    vname = 'mjj_vs_deta'
    if vname not in trig_vars:
        return

    base_x = np.array(output['Baseline'][vname]['x'])
    base_y = np.array(output['Baseline'][vname]['y'])
    pass_x = np.array(output['VBFPass'][vname]['x'])
    pass_y = np.array(output['VBFPass'][vname]['y'])

    sig_b = (base_x > mjj_min) & (base_y > deta_min)
    sig_p = (pass_x > mjj_min) & (pass_y > deta_min)

    rows = [['Region', 'Trigger', 'Efficiency', 'Eff_Err', 'Pass', 'Total']]
    for label, bc, pc in [
        ('Inclusive',                       np.ones(len(base_x), dtype=bool), np.ones(len(pass_x), dtype=bool)),
        (f'mjj>{mjj_min},deta>{deta_min}',  sig_b, sig_p),
    ]:
        nb = int(np.sum(bc)); npass = int(np.sum(pc))
        eff_arr, err_arr = _binomial_eff_err(np.array([npass]), np.array([nb]))
        e = float(eff_arr[0]) if np.isfinite(eff_arr[0]) else 0.0
        u = float(err_arr[0]) if np.isfinite(err_arr[0]) else 0.0
        rows.append([label, f"HLT_{VBF_SIGNAL_TRIGGER}", f"{e:.4f}", f"{u:.4f}", str(npass), str(nb)])
        print(f"  {label:35s} | Eff: {e:.2%} ± {u:.2%} ({npass}/{nb})")

    csv_path = os.path.join(save_dir, f"vbf_parking_trig_eff_{year}{('_'+suffix) if suffix else ''}.csv")
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f"Saved CSV: {csv_path}")


def generate_sf_csv(output_data, output_mc, trig_vars, save_dir, year,
                    mjj_min=720, deta_min=3.0):
    """SF table for the zmumu_vbf method (inclusive + signal region)."""
    os.makedirs(save_dir, exist_ok=True)
    vname = 'mjj_vs_deta'
    if vname not in trig_vars:
        return

    def _eff(out, mjj_cut, deta_cut):
        bx = np.array(out['Baseline'][vname]['x']); by = np.array(out['Baseline'][vname]['y'])
        px = np.array(out['VBFPass'][vname]['x']);  py = np.array(out['VBFPass'][vname]['y'])
        bcut = (bx > mjj_cut) & (by > deta_cut) if mjj_cut else np.ones(len(bx), dtype=bool)
        pcut = (px > mjj_cut) & (py > deta_cut) if mjj_cut else np.ones(len(px), dtype=bool)
        nb = int(np.sum(bcut)); npass = int(np.sum(pcut))
        e_arr, u_arr = _binomial_eff_err(np.array([npass]), np.array([nb]))
        return (float(e_arr[0]) if np.isfinite(e_arr[0]) else 0.0,
                float(u_arr[0]) if np.isfinite(u_arr[0]) else 0.0,
                npass, nb)

    rows = [['Region', 'Eff_data', 'Err_data', 'Eff_mc', 'Err_mc', 'SF', 'SF_err']]
    for label, mc, dc in [('Inclusive', None, None),
                          (f'mjj>{mjj_min},deta>{deta_min}', mjj_min, deta_min)]:
        ed, ud, kpd, kbd = _eff(output_data, mc, dc)
        em, um, kpm, kbm = _eff(output_mc,   mc, dc)
        sf = ed / em if em > 0 else 1.0
        rel = (((1 - ed) / kpd) if kpd > 0 else np.nan) + (((1 - em) / kpm) if kpm > 0 else np.nan)
        sf_err = abs(sf) * np.sqrt(max(rel, 0.0)) if np.isfinite(rel) else 0.0
        rows.append([label, f"{ed:.4f}", f"{ud:.4f}", f"{em:.4f}", f"{um:.4f}",
                     f"{sf:.4f}", f"{sf_err:.4f}"])
        print(f"  {label:35s} | data {ed:.2%}±{ud:.2%}  mc {em:.2%}±{um:.2%}  SF {sf:.3f}±{sf_err:.3f}")

    csv_path = os.path.join(save_dir, f"vbf_parking_trig_sf_{year}.csv")
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f"Saved CSV: {csv_path}")


# ─── Fileset construction ─────────────────────────────────────────────────────

def _truncate(samples, n):
    return samples[:n] if (n and len(samples) > n) else samples


def build_filesets(year, method, test_n_files, mc_test_n_files):
    """
    Returns dict { 'data': [...], 'mc': [...]  }  with file lists ready for coffea.
    'mc' is empty for the mu_tag_parkingvbf method.
    """
    filesets = {'data': [], 'mc': []}

    if method == 'mu_tag_parkingvbf':
        infile_path = f"infiles/{year}/{year}_ParkingVBFData.json"
        if not os.path.exists(infile_path):
            raise FileNotFoundError(
                f"{infile_path} missing. Run make_parking_vbf_infiles.py first."
            )
        with open(infile_path) as f:
            d = json.load(f)
        for k in PARKING_VBF_KEYS.get(year, []):
            files = d.get(k, [])
            filesets['data'].extend(_truncate(files, test_n_files))
        return filesets

    # zmumu_vbf
    muon_path = f"infiles/{year}/{year}_MuonData.json"
    zjets_path = f"infiles/{year}/{year}_Zjets.json"
    if not os.path.exists(muon_path):
        raise FileNotFoundError(f"{muon_path} missing.")
    if not os.path.exists(zjets_path):
        raise FileNotFoundError(f"{zjets_path} missing.")

    with open(muon_path) as f:
        mu_data = json.load(f)
    with open(zjets_path) as f:
        zjets = json.load(f)

    for k in MUON_DATA_KEYS.get(year, []):
        files = mu_data.get(k, [])
        filesets['data'].extend(_truncate(files, test_n_files))

    dy_files = zjets.get(DY_MC_KEY, [])
    filesets['mc'].extend(_truncate(dy_files, mc_test_n_files))

    return filesets


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="VBF Parking Trigger Efficiency / SF")
    parser.add_argument('--year', default='2023',
                        choices=['2022', '2022EE', '2023', '2023BPix'],
                        help="Period to process. Note: ParkingVBF NanoAOD only "
                             "exists for 2023 / 2023BPix.")
    parser.add_argument('--method', default='zmumu_vbf',
                        choices=['mu_tag_parkingvbf', 'zmumu_vbf'],
                        help="Selection method. 'zmumu_vbf' = Z(mumu)+2 VBF jets "
                             "(produces SF maps); 'mu_tag_parkingvbf' = legacy data-only.")
    parser.add_argument('--test', action='store_true',
                        help="Test mode: cap files per data/MC source.")
    parser.add_argument('--test-n-data', type=int, default=10,
                        help="Files-per-data-key cap in --test mode (default 10).")
    parser.add_argument('--test-n-mc', type=int, default=20,
                        help="Files cap on the DY MC in --test mode (default 20).")
    parser.add_argument('--workers', type=int, default=2,
                        help="Number of FuturesExecutor workers.")
    parser.add_argument('--outdir', default='.',
                        help="Base output directory (figures & coffea files go under here).")
    args = parser.parse_args()

    year = args.year
    method = args.method
    test_n_d = args.test_n_data if args.test else None
    test_n_m = args.test_n_mc   if args.test else None

    base_outdir   = os.path.abspath(args.outdir)
    figures_dir   = os.path.join(base_outdir, "figures_parking_vbf", year, method)
    output_dir    = os.path.join(base_outdir, "output", year, method)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(output_dir,  exist_ok=True)

    print(f"\n{'='*60}")
    print(f"VBF Trigger Efficiency Study — year={year}, method={method}")
    print(f"Signal trigger: HLT_{VBF_SIGNAL_TRIGGER}")
    if args.test:
        print(f"TEST MODE: data<= {test_n_d} files/key, mc<= {test_n_m} files")
    print(f"{'='*60}\n")

    trig_vars = trig_vars_for_method(method)
    filesets = build_filesets(year, method, test_n_d, test_n_m)

    print(f"data files: {len(filesets['data'])}")
    print(f"mc   files: {len(filesets['mc'])}")
    if not filesets['data']:
        print("ERROR: no data files loaded.")
        import sys; sys.exit(1)

    iterative_run = processor.Runner(
        executor=processor.FuturesExecutor(compression=None, workers=args.workers),
        schema=NanoAODSchema,
        skipbadfiles=True,
        savemetrics=True,
    )

    # ── Run on data ───────────────────────────────────────────────────────────
    print("\nRunning processor on DATA ...")
    out_data, _ = iterative_run(
        {'Data': filesets['data']},
        treename="Events",
        processor_instance=VBFTrigEffProcessor(year, method, trig_vars),
    )
    util.save(out_data, os.path.join(output_dir, f"parking_vbf_eff_{year}_data.coffea"))

    plot_2d_efficiency(out_data, trig_vars, save_dir=figures_dir, year=year, label_suffix='data')
    generate_efficiency_csv(out_data, trig_vars, save_dir=output_dir, year=year, suffix='data')

    # ── Run on MC for the SF (zmumu_vbf method only) ──────────────────────────
    if method == 'zmumu_vbf' and filesets['mc']:
        print("\nRunning processor on MC (DYJetsToLL) ...")
        out_mc, _ = iterative_run(
            {'MC': filesets['mc']},
            treename="Events",
            processor_instance=VBFTrigEffProcessor(year, method, trig_vars),
        )
        util.save(out_mc, os.path.join(output_dir, f"parking_vbf_eff_{year}_mc.coffea"))

        plot_2d_efficiency(out_mc, trig_vars, save_dir=figures_dir, year=year, label_suffix='mc')
        generate_efficiency_csv(out_mc, trig_vars, save_dir=output_dir, year=year, suffix='mc')

        plot_2d_data_mc_sf(out_data, out_mc, trig_vars, save_dir=figures_dir, year=year)
        generate_sf_csv(out_data, out_mc, trig_vars, save_dir=output_dir, year=year)
    elif method == 'zmumu_vbf':
        print("WARNING: --method zmumu_vbf requested but DY MC fileset is empty. SF maps skipped.")

    print(f"\nDone. Figures: {figures_dir}/")
    print(f"      Output:  {output_dir}/")
