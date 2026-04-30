import warnings
warnings.filterwarnings('ignore')
import os
import awkward as ak
import uproot
import hist
import numpy as np
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.processor import dict_accumulator, list_accumulator
from coffea import util

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep as hep
plt.style.use(hep.style.ROOT)

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium'}
pylab.rcParams.update(params)

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 5
import itertools
import json


# txbb+txcc working point for pass/fail region split
TXBB_WP = 0.82

# Run 2 AN kinematic region for trigger efficiency/SF reporting
PT_MIN_AN = 450.0
PT_MAX_AN = 1000.0
MSD_MIN_AN = 40.0
MSD_MAX_AN = 300.0

# POG correctionlib paths (used for v14+ JetID correction)
POG_CORRECTION_PATH = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration"
POG_JSONS = {
    "jetid": ("JME", "jetid.json.gz"),
}
POG_YEARS = {
    "2022": "2022_Summer22",
    "2022EE": "2022_Summer22EE",
    "2023": "2023_Summer23",
    "2023BPix": "2023_Summer23BPix",
}

_TXCC_FALLBACK_WARNED = False
_JETID_FALLBACK_WARNED = False


def ak_clip(arr, min_value, max_value):
    """Clip awkward array values to the valid range expected by correctionlib payloads."""
    return ak.where(arr < min_value, min_value, ak.where(arr > max_value, max_value, arr))


def get_pog_json(obj, year):
    """Return path to a POG json file for a given object/year."""
    group, filename = POG_JSONS[obj]
    year_key = POG_YEARS[year]
    return f"{POG_CORRECTION_PATH}/POG/{group}/{year_key}/{filename}"


def get_txbb_txcc_score(jets):
    """
    Preferred discriminator for pass/fail split: txbb + txcc.
    Fallback to particleNet_XbbVsQCD when txbb/txcc branches are unavailable.
    """
    global _TXCC_FALLBACK_WARNED
    fields = set(jets.fields)
    if "particleNet_Xbb" in fields and "particleNet_Xcc" in fields:
        xbb = ak.fill_none(jets.particleNet_Xbb, 0.0)
        xcc = ak.fill_none(jets.particleNet_Xcc, 0.0)
        return xbb + xcc

    if "particleNet_XbbVsQCD" in fields:
        if not _TXCC_FALLBACK_WARNED:
            print("Warning: particleNet_Xbb/Xcc not found. Falling back to particleNet_XbbVsQCD.")
            _TXCC_FALLBACK_WARNED = True
        return ak.fill_none(jets.particleNet_XbbVsQCD, 0.0)

    if not _TXCC_FALLBACK_WARNED:
        print("Warning: no txbb/txcc discriminator fields found. Using score=0.")
        _TXCC_FALLBACK_WARNED = True
    return ak.zeros_like(jets.pt)


def corrected_jetid_tight_mask(jets, jet_type, year, use_jetid_correction=True):
    """
    For v14+ periods (2023/2023BPix), compute Tight JetID from correctionlib.
    Falls back to NanoAOD isTight if correction payload/inputs are unavailable.
    """
    global _JETID_FALLBACK_WARNED

    if not use_jetid_correction or year not in {"2023", "2023BPix"}:
        return jets.isTight

    if jet_type == "AK8":
        corr_name = "AK8PUPPI_Tight"
    elif jet_type == "AK4":
        corr_name = "AK4PUPPI_Tight"
    else:
        return jets.isTight

    required_fields = {
        "eta", "chHEF", "neHEF", "chEmEF", "neEmEF", "muEF",
        "chMultiplicity", "neMultiplicity"
    }
    if not required_fields.issubset(set(jets.fields)):
        if not _JETID_FALLBACK_WARNED:
            print("Warning: missing JetID input branches. Falling back to NanoAOD isTight.")
            _JETID_FALLBACK_WARNED = True
        return jets.isTight

    try:
        import correctionlib  # local import to keep script usable when unavailable

        evaluator = correctionlib.CorrectionSet.from_file(get_pog_json("jetid", year))
        flat_jets = ak.flatten(jets)
        njets = ak.num(jets)

        eta = ak_clip(flat_jets.eta, -4.7, 4.7)
        chHEF = ak_clip(flat_jets.chHEF, 0.0, 1.0)
        neHEF = ak_clip(flat_jets.neHEF, 0.0, 1.0)
        chEmEF = ak_clip(flat_jets.chEmEF, 0.0, 1.0)
        neEmEF = ak_clip(flat_jets.neEmEF, 0.0, 1.0)
        muEF = ak_clip(flat_jets.muEF, 0.0, 1.0)
        chMult = ak_clip(flat_jets.chMultiplicity, 0.0, 100.0)
        neMult = ak_clip(flat_jets.neMultiplicity, 0.0, 100.0)
        nConst = chMult + neMult

        mask_flat = evaluator[corr_name].evaluate(
            eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, nConst
        )
        return ak.values_astype(ak.unflatten(mask_flat, njets), bool)
    except Exception as exc:
        if not _JETID_FALLBACK_WARNED:
            print(f"Warning: JetID correction unavailable ({exc}). Falling back to NanoAOD isTight.")
            _JETID_FALLBACK_WARNED = True
        return jets.isTight


def round_robin_sample_lists(list_of_lists, n_total):
    """
    Round-robin pick up to n_total items across multiple lists.
    Preserves per-list order and avoids bias toward earlier lists.
    """
    if n_total <= 0:
        return []
    sampled = []
    idx = 0
    # Continue until we collect enough or no list has remaining items.
    while len(sampled) < n_total:
        made_progress = False
        for lst in list_of_lists:
            if idx < len(lst):
                sampled.append(lst[idx])
                made_progress = True
                if len(sampled) >= n_total:
                    break
        if not made_progress:
            break
        idx += 1
    return sampled


# 2D trigger variables configuration
trig_vars_2d = {
    'pt_vs_msd': {
        'label_x': "Leading Jet $m_{SD}$ [GeV]",
        'label_y': "Leading Jet $p_{T}$ [GeV]",
        'axis_x': hist.axis.Regular(bins=15, start=0, stop=300, name="msd", label="Leading Jet $m_{SD}$ [GeV]"),
        'axis_y': hist.axis.Regular(bins=20, start=200, stop=1000, name="pt", label="Leading Jet $p_{T}$ [GeV]"),
        'proc_x': lambda events: ak.fill_none(ak.pad_none(events.FatJet.msoftdrop, 1, clip=True)[:, 0], np.nan),
        'proc_y': lambda events: ak.fill_none(ak.pad_none(events.FatJet.pt, 1, clip=True)[:, 0], np.nan)
    }
}

trigger_dict_periods = {
    '2022': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65'
    ],
    '2022EE': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65'
    ],
    '2023': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        # VBF_DiPFJet125_45_Mjj720_Detajj3p0 uses a different calibration strategy.
        # Excluded from trigger SF OR and per-trigger table in this workflow.
    ],
    '2023BPix': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        # VBF_DiPFJet125_45_Mjj720_Detajj3p0 uses a different calibration strategy.
        # It is also in ParkingVBF stream only and not stored in Muon0 PD for Run2023D.
    ]
}

# Muon reference triggers for SF calculation
muon_triggers = {
    '2022': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2022EE': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023BPix': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
}

def create_baseline_selection_mask_new(events, tag, year, txbb_region='inclusive', use_jetid_correction=True):
    """
    txbb_region: 'pass'      -> candidatejet (txbb+txcc) >= TXBB_WP (0.82)
                 'fail'      -> candidatejet (txbb+txcc) <  TXBB_WP (with loose preselection > 0.4)
    """
    n_events = len(events)
    baseline_mask = np.ones(n_events, dtype=bool)

    # Step 1: Common jet selection
    jets_all = events.Jet[events.Jet.pt > 15]
    has_jets = ak.num(jets_all) >= 1
    baseline_mask &= ak.to_numpy(has_jets)

    # Step 2: General selection on FatJets (loose preselection, txbb+txcc > 0.4)
    fatjets = events.FatJet
    fatjet_score = get_txbb_txcc_score(fatjets)
    fatjet_mask = (fatjets.pt > 250) & (abs(fatjets.eta) < 2.5) & (fatjet_score > 0.4)
    has_fatjet = ak.sum(fatjet_mask, axis=1) >= 1
    baseline_mask &= ak.to_numpy(has_fatjet)

    fatjet_tight = corrected_jetid_tight_mask(
        fatjets, jet_type="AK8", year=year, use_jetid_correction=use_jetid_correction
    )
    candidatejet = fatjets[
            (fatjets.pt > 200)
            & (abs(fatjets.eta) < 2.5)
            & fatjet_tight
        ]
    candidatejet = candidatejet[:, :2]
    candidate_scores = get_txbb_txcc_score(candidatejet)
    candidatejet = ak.firsts(candidatejet[ak.argmax(candidate_scores, axis=1, keepdims=True)])
    candidate_score = ak.fill_none(get_txbb_txcc_score(candidatejet), -1)

    # Step 2b: txbb pass/fail split on the candidatejet
    if txbb_region == 'pass':
        txbb_cut = candidate_score >= TXBB_WP
        baseline_mask &= ak.to_numpy(txbb_cut)
    elif txbb_region == 'fail':
        txbb_cut = candidate_score < TXBB_WP
        baseline_mask &= ak.to_numpy(txbb_cut)

    # Step 3: Optional VBF topology selection
    # NOTE: For trigger SF measurement, use tag='Inclusive' (no VBF cuts)
    # VBF cuts should NOT be used for trigger SF - triggers are production-mode agnostic
    if tag == "VBF":
        jets = events.Jet
        jet_tight = corrected_jetid_tight_mask(
            jets, jet_type="AK4", year=year, use_jetid_correction=use_jetid_correction
        )
        jets = jets[
            (jets.pt > 30.)
            & (abs(jets.eta) < 5.0)
            & jet_tight
        ]
        jets = jets[:, :4]
        has_two_jets = ak.num(jets) >= 2

        dR = jets.delta_r(candidatejet)
        ak4_outside_ak8 = jets[dR > 0.8]

        jet1 = ak4_outside_ak8[:, 0:1]
        jet2 = ak4_outside_ak8[:, 1:2]

        deta = abs(ak.firsts(jet1).eta - ak.firsts(jet2).eta)
        mjj = ( ak.firsts(jet1) + ak.firsts(jet2) ).mass

        has_valid_pair = ((deta > 3.5) & (mjj > 1000))
        vbf_mask = has_two_jets & has_valid_pair
        baseline_mask &= ak.to_numpy(vbf_mask)

    # Step 4: Muon veto - require at least one loose muon far from leading FatJet
    muons = events.Muon
    loose_muon_selection = (
        (muons.pt > 25) & (abs(muons.eta) < 2.4) & (muons.looseId) & (muons.pfRelIso04_all < 0.25)
    )
    loose_muons = muons[loose_muon_selection]
    has_loose_muon = ak.num(loose_muons) >= 1

    # Veto jets near the muon: compute ΔR(mu, leading FatJet) and require > 0.8
    leading_fj = ak.firsts(ak.pad_none(events.FatJet, 1, clip=True))
    dphi = np.abs(((loose_muons.phi - leading_fj.phi + np.pi) % (2 * np.pi)) - np.pi)
    deta = loose_muons.eta - leading_fj.eta
    dr = np.sqrt(dphi**2 + deta**2)
    far_enough = ak.fill_none(ak.any(dr > 0.8, axis=1), True)

    muon_veto = has_loose_muon & far_enough
    baseline_mask &= ak.to_numpy(muon_veto)

    assert len(baseline_mask) == n_events, "Baseline mask length does not match number of events!"
    return baseline_mask


class ScaleFactor2DProcessor(processor.ProcessorABC):
    """
    Calculate 2D trigger scale factors (leading jet pT vs mass) by comparing data and MC efficiencies.
    """
    def __init__(self, year, trig_vars_2d, baseline_key='Inclusive', txbb_region='inclusive', use_jetid_correction=True):
        self.year = str(year)
        self.trig_vars_2d = trig_vars_2d
        self.baseline_key = baseline_key
        self.txbb_region = txbb_region
        self.use_jetid_correction = use_jetid_correction

        # Get trigger list and muon triggers for this year
        self.triggers = trigger_dict_periods.get(self.year, [])
        self.muon_triggers = muon_triggers.get(self.year, [])

        # Initialize output using accumulators
        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                var_name: dict_accumulator({
                    'x': list_accumulator(),
                    'y': list_accumulator()
                })
                for var_name in trig_vars_2d
            }),
            'Numerator': dict_accumulator({
                var_name: dict_accumulator({
                    'pass_x': list_accumulator(),
                    'pass_y': list_accumulator()
                })
                for var_name in trig_vars_2d
            }),
            'PerTrigger': dict_accumulator({
                trigger: dict_accumulator({
                    var_name: dict_accumulator({
                        'pass_x': list_accumulator(),
                        'pass_y': list_accumulator()
                    })
                    for var_name in trig_vars_2d
                })
                for trigger in self.triggers
            })
        })

    def process(self, events):
        output = self.output
        dataset = events.metadata['dataset']
        is_mc = hasattr(events, "genWeight")

        # First, apply reference muon triggers
        reference_trigger_paths = self.muon_triggers
        fired_reference_triggers = ak.zeros_like(events.event, dtype=bool)
        for trig in reference_trigger_paths:
            tfield = trig[4:] if trig.startswith('HLT_') else trig
            if tfield in events.HLT.fields:
                fired_reference_triggers = fired_reference_triggers | events.HLT[tfield]

        events = events[fired_reference_triggers]

        if len(events) == 0:
            return output

        # Apply baseline selection
        baseline = create_baseline_selection_mask_new(
            events,
            tag=self.baseline_key,
            year=self.year,
            txbb_region=self.txbb_region,
            use_jetid_correction=self.use_jetid_correction
        )

        # Compute 2D variables
        variables = {}
        for var_name, var_info in self.trig_vars_2d.items():
            var_x = var_info['proc_x'](events)
            var_y = var_info['proc_y'](events)
            variables[var_name] = {'x': var_x, 'y': var_y}

        # Create a valid mask where all variables are not NaN
        valid_vars_mask = np.ones(len(events), dtype=bool)
        for var_data in variables.values():
            valid_vars_mask &= ~np.isnan(ak.to_numpy(var_data['x']))
            valid_vars_mask &= ~np.isnan(ak.to_numpy(var_data['y']))

        # Combined selection for baseline (denominator - no HLT trigger requirement)
        selection_baseline = baseline & valid_vars_mask

        # Fill baseline histograms (denominator)
        for var_name in self.trig_vars_2d:
            var_x = variables[var_name]['x']
            var_y = variables[var_name]['y']

            var_x_baseline = var_x[selection_baseline]
            var_y_baseline = var_y[selection_baseline]

            output['Baseline'][var_name]['x'].extend(ak.to_numpy(var_x_baseline).tolist())
            output['Baseline'][var_name]['y'].extend(ak.to_numpy(var_y_baseline).tolist())

        # Build numerator from the OR of all HLT triggers configured for this year
        if len(self.triggers) > 0:
            combined_trigger_mask = ak.zeros_like(events.event, dtype=bool)
            for trigger in self.triggers:
                tfield = trigger[4:] if trigger.startswith('HLT_') else trigger
                if tfield in events.HLT.fields:
                    combined_trigger_mask = combined_trigger_mask | events.HLT[tfield]
        else:
            combined_trigger_mask = ak.zeros_like(baseline, dtype=bool)

        selection_pass = baseline & combined_trigger_mask & valid_vars_mask

        # Fill numerator histograms (combined trigger OR)
        for var_name in self.trig_vars_2d:
            var_x = variables[var_name]['x']
            var_y = variables[var_name]['y']

            var_x_pass = var_x[selection_pass]
            var_y_pass = var_y[selection_pass]

            output['Numerator'][var_name]['pass_x'].extend(ak.to_numpy(var_x_pass).tolist())
            output['Numerator'][var_name]['pass_y'].extend(ak.to_numpy(var_y_pass).tolist())

        # Fill per-trigger histograms
        for trigger in self.triggers:
            tfield = trigger[4:] if trigger.startswith('HLT_') else trigger
            if tfield in events.HLT.fields:
                trigger_mask = events.HLT[tfield]
                selection_pass_single = baseline & trigger_mask & valid_vars_mask

                for var_name in self.trig_vars_2d:
                    var_x = variables[var_name]['x']
                    var_y = variables[var_name]['y']

                    var_x_pass = var_x[selection_pass_single]
                    var_y_pass = var_y[selection_pass_single]

                    output['PerTrigger'][trigger][var_name]['pass_x'].extend(ak.to_numpy(var_x_pass).tolist())
                    output['PerTrigger'][trigger][var_name]['pass_y'].extend(ak.to_numpy(var_y_pass).tolist())

        return output

    def postprocess(self, accumulator):
        # Convert lists to arrays
        for key in accumulator:
            if key == 'Baseline':
                for var_name in accumulator[key]:
                    accumulator[key][var_name]['x'] = np.array(accumulator[key][var_name]['x'])
                    accumulator[key][var_name]['y'] = np.array(accumulator[key][var_name]['y'])
            elif key == 'Numerator':
                for var_name in accumulator[key]:
                    accumulator[key][var_name]['pass_x'] = np.array(accumulator[key][var_name]['pass_x'])
                    accumulator[key][var_name]['pass_y'] = np.array(accumulator[key][var_name]['pass_y'])
            elif key == 'PerTrigger':
                for trigger in accumulator[key]:
                    for var_name in accumulator[key][trigger]:
                        accumulator[key][trigger][var_name]['pass_x'] = np.array(accumulator[key][trigger][var_name]['pass_x'])
                        accumulator[key][trigger][var_name]['pass_y'] = np.array(accumulator[key][trigger][var_name]['pass_y'])
        return accumulator


def plot_2d_scale_factors(output_data, output_mc, trig_vars_2d, save_dir=None, year=2022, show_img=True):
    """
    Plot 2D scale factors by comparing data and MC trigger efficiencies.

    Parameters:
    - output_data: dict - Output from ScaleFactor2DProcessor for data
    - output_mc: dict - Output from ScaleFactor2DProcessor for MC
    - trig_vars_2d: dict - Dictionary of 2D trigger variables
    - save_dir: str, optional - Directory to save the figures
    - year: int, optional - Year label for the plots
    - show_img: bool, optional - Whether to display the images
    """
    if save_dir is None:
        save_dir = "./figures_sf_2d"
    os.makedirs(save_dir, exist_ok=True)

    plt.rcParams.update({"font.size": 12})

    for var_name, var_info in trig_vars_2d.items():
        # Get bin edges
        bins_x = var_info['axis_x'].edges
        bins_y = var_info['axis_y'].edges

        # Get baseline (total) data for both data and MC
        baseline_x_data = np.array(output_data['Baseline'][var_name]['x'])
        baseline_y_data = np.array(output_data['Baseline'][var_name]['y'])
        baseline_x_mc = np.array(output_mc['Baseline'][var_name]['x'])
        baseline_y_mc = np.array(output_mc['Baseline'][var_name]['y'])

        # Get numerator data for both data and MC
        pass_x_data = np.array(output_data['Numerator'][var_name]['pass_x'])
        pass_y_data = np.array(output_data['Numerator'][var_name]['pass_y'])
        pass_x_mc = np.array(output_mc['Numerator'][var_name]['pass_x'])
        pass_y_mc = np.array(output_mc['Numerator'][var_name]['pass_y'])

        # Create 2D histograms for baseline (denominator)
        hist_baseline_data, xedges, yedges = np.histogram2d(
            baseline_x_data, baseline_y_data,
            bins=[bins_x, bins_y]
        )
        hist_baseline_mc, _, _ = np.histogram2d(
            baseline_x_mc, baseline_y_mc,
            bins=[bins_x, bins_y]
        )

        # Create 2D histograms for passing events (numerator)
        hist_pass_data, _, _ = np.histogram2d(
            pass_x_data, pass_y_data,
            bins=[bins_x, bins_y]
        )
        hist_pass_mc, _, _ = np.histogram2d(
            pass_x_mc, pass_y_mc,
            bins=[bins_x, bins_y]
        )

        # Calculate efficiencies
        with np.errstate(divide='ignore', invalid='ignore'):
            eff_data = np.divide(
                hist_pass_data, hist_baseline_data,
                out=np.zeros_like(hist_pass_data, dtype=float),
                where=hist_baseline_data > 0
            )
            eff_mc = np.divide(
                hist_pass_mc, hist_baseline_mc,
                out=np.zeros_like(hist_pass_mc, dtype=float),
                where=hist_baseline_mc > 0
            )

        # Calculate scale factor
        valid_sf = (hist_baseline_data > 0) & (hist_baseline_mc > 0) & (eff_mc > 0)
        sf = np.divide(
            eff_data, eff_mc,
            out=np.ones_like(eff_data, dtype=float),
            where=valid_sf
        )

        # Mask invalid SF bins
        sf[~valid_sf] = np.nan

        # Calculate total efficiency and SF
        total_pass_data = np.sum(hist_pass_data)
        total_baseline_data = np.sum(hist_baseline_data)
        total_eff_data = total_pass_data / total_baseline_data if total_baseline_data > 0 else 0.0

        total_pass_mc = np.sum(hist_pass_mc)
        total_baseline_mc = np.sum(hist_baseline_mc)
        total_eff_mc = total_pass_mc / total_baseline_mc if total_baseline_mc > 0 else 0.0

        total_sf = total_eff_data / total_eff_mc if total_eff_mc > 0 else 1.0

        # Create figure with 3 subplots (Data Eff, MC Eff, SF)
        fig, axes = plt.subplots(1, 3, figsize=(24, 7))

        # Plot 1: Data Efficiency
        ax1 = axes[0]
        im1 = ax1.pcolormesh(
            xedges, yedges, eff_data.T,
            cmap='RdYlGn',
            vmin=0, vmax=1,
            shading='flat'
        )
        cbar1 = plt.colorbar(im1, ax=ax1)
        cbar1.set_label('Efficiency', fontsize=12)
        ax1.set_xlabel(var_info['label_x'], fontsize=12)
        ax1.set_ylabel(var_info['label_y'], fontsize=12)
        ax1.set_title(f"Data Efficiency\nOverall: {total_eff_data:.2%}", fontsize=12)
        hep.cms.label(ax=ax1, data=True, year=year, com="13.6", fontsize=10)

        # Plot 2: MC Efficiency
        ax2 = axes[1]
        im2 = ax2.pcolormesh(
            xedges, yedges, eff_mc.T,
            cmap='RdYlGn',
            vmin=0, vmax=1,
            shading='flat'
        )
        cbar2 = plt.colorbar(im2, ax=ax2)
        cbar2.set_label('Efficiency', fontsize=12)
        ax2.set_xlabel(var_info['label_x'], fontsize=12)
        ax2.set_ylabel(var_info['label_y'], fontsize=12)
        ax2.set_title(f"MC Efficiency\nOverall: {total_eff_mc:.2%}", fontsize=12)
        hep.cms.label(ax=ax2, data=False, year=year, com="13.6", label="Simulation", fontsize=10)

        # Plot 3: Scale Factor
        ax3 = axes[2]
        im3 = ax3.pcolormesh(
            xedges, yedges, sf.T,
            cmap='RdBu_r',
            vmin=0.8, vmax=1.2,
            shading='flat'
        )
        cbar3 = plt.colorbar(im3, ax=ax3)
        cbar3.set_label('Scale Factor (Data/MC)', fontsize=12)
        ax3.set_xlabel(var_info['label_x'], fontsize=12)
        ax3.set_ylabel(var_info['label_y'], fontsize=12)
        ax3.set_title(f"Scale Factor\nOverall: {total_sf:.3f}", fontsize=12)
        hep.cms.label(ax=ax3, data=True, year=year, com="13.6", fontsize=10)

        plt.tight_layout()

        # Save figure
        save_path = os.path.join(save_dir, f"scale_factor_{var_name}_2d.png")
        print(f"Saving 2D SF plot to: {save_path}")
        plt.savefig(save_path, dpi=200, bbox_inches='tight')

        if show_img:
            plt.show()
        else:
            print(f"Saved {save_path}")

        plt.close(fig)


def plot_1d_efficiency_projections(output_data, output_mc, trig_vars_2d, triggers,
                                   save_dir=None, year='2022', show_img=False,
                                   pt_min_plateau=PT_MIN_AN, msd_min_plot=MSD_MIN_AN):
    """
    Plot 1D trigger efficiency projections matching AN Figs 3-4, 6-7, 10-11 style.

    Produces two-panel figures (data left, MC right):
      - efficiency_vs_pt_*.png  : efficiency vs leading jet pT (full range, all mSD values)
      - efficiency_vs_msd_*.png : efficiency vs leading jet mSD (full range, pT >= pt_min_plateau)

    The kinematic range used for *calculating* the plateau efficiency values in tables is
    450 <= pT < 1000 GeV, 40 <= mSD < 300 GeV (enforced in generate_per_trigger_efficiency_csv).
    Here we plot the full range so the turn-on is visible, matching the AN figure style.
    """
    if save_dir is None:
        save_dir = "./figures_sf_2d"
    os.makedirs(save_dir, exist_ok=True)

    var_name = 'pt_vs_msd'
    if var_name not in trig_vars_2d:
        print(f"Warning: {var_name} not in trig_vars_2d, skipping 1D projections.")
        return

    # Full-range bins for plotting (to show turn-on)
    pt_bins  = np.linspace(200, 1000, 33)   # ~25 GeV per bin
    msd_bins = np.linspace(0, 300, 31)      # 10 GeV per bin

    # Build ordered list: individual triggers first, then OR ("Soup")
    trigger_names  = list(triggers)
    trigger_labels = [t.replace('HLT_', '') for t in trigger_names]
    all_keys   = trigger_names + ['OR']
    all_labels = trigger_labels + ['Soup']

    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = [prop_cycle[i % len(prop_cycle)] for i in range(len(all_keys))]

    def _eff_arrays(output, key, axis, bins, plateau_mask_fn=None):
        """Return (bin_centers, efficiency) arrays for one trigger key and one axis."""
        if key == 'OR':
            pass_x = np.array(output['Numerator'][var_name]['pass_x'])
            pass_y = np.array(output['Numerator'][var_name]['pass_y'])
        else:
            if 'PerTrigger' not in output or key not in output['PerTrigger']:
                return None, None
            pass_x = np.array(output['PerTrigger'][key][var_name]['pass_x'])
            pass_y = np.array(output['PerTrigger'][key][var_name]['pass_y'])

        base_x = np.array(output['Baseline'][var_name]['x'])
        base_y = np.array(output['Baseline'][var_name]['y'])

        if plateau_mask_fn is not None:
            pm_base = plateau_mask_fn(base_x, base_y)
            pm_pass = plateau_mask_fn(pass_x, pass_y)
            base_x, base_y = base_x[pm_base], base_y[pm_base]
            pass_x, pass_y = pass_x[pm_pass], pass_y[pm_pass]

        if axis == 'pt':
            n_base, _ = np.histogram(base_y, bins=bins)
            n_pass, _ = np.histogram(pass_y, bins=bins)
        else:  # msd
            n_base, _ = np.histogram(base_x, bins=bins)
            n_pass, _ = np.histogram(pass_x, bins=bins)

        with np.errstate(divide='ignore', invalid='ignore'):
            eff = np.where(n_base > 0, n_pass / n_base.astype(float), np.nan)

        centers = 0.5 * (bins[:-1] + bins[1:])
        return centers, eff

    def _draw_panel(ax, output, is_data, axis, bins, plateau_mask_fn=None):
        for i, (key, label) in enumerate(zip(all_keys, all_labels)):
            centers, eff = _eff_arrays(output, key, axis, bins, plateau_mask_fn)
            if centers is None:
                continue
            ls  = '-'
            mk  = '+' if label == 'Soup' else None
            mks = 8   if label == 'Soup' else None
            ax.plot(centers, eff, label=label, color=colors[i],
                    linestyle=ls, marker=mk, markersize=mks, linewidth=2)
        ax.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax.set_ylim(0, 1.4)
        ax.legend(fontsize=8, loc='lower right', ncol=1)
        hep.cms.label(ax=ax, data=is_data, year=year, com="13.6", fontsize=10)

    # ------------------------------------------------------------------ #
    #  Figure 1: Efficiency vs leading jet pT  (no mSD cut, full range)   #
    # ------------------------------------------------------------------ #
    fig, (ax_d, ax_mc) = plt.subplots(1, 2, figsize=(18, 7))
    _draw_panel(ax_d,  output_data, is_data=True,  axis='pt', bins=pt_bins)
    _draw_panel(ax_mc, output_mc,   is_data=False, axis='pt', bins=pt_bins)
    for ax, title in [(ax_d, "Triggers ($\\mu$ ref, all $m_{SD}$)"),
                      (ax_mc, "Triggers (QCD MC, all $m_{SD}$)")]:
        ax.set_xlabel("Leading jet $p_T$ [GeV]", fontsize=13)
        ax.set_ylabel("Trigger Efficiency", fontsize=13)
        ax.set_xlim(pt_bins[0], pt_bins[-1])
        ax.set_title(title, fontsize=11)
    plt.tight_layout()
    save_path = os.path.join(save_dir, "efficiency_vs_pt.png")
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    if show_img:
        plt.show()
    else:
        print(f"Saved {save_path}")
    plt.close(fig)

    # ------------------------------------------------------------------ #
    #  Figure 2: Efficiency vs mSD  (pT >= pt_min_plateau)               #
    # ------------------------------------------------------------------ #
    plateau_fn = lambda px, py: py >= pt_min_plateau
    fig, (ax_d, ax_mc) = plt.subplots(1, 2, figsize=(18, 7))
    _draw_panel(ax_d,  output_data, is_data=True,  axis='msd', bins=msd_bins, plateau_mask_fn=plateau_fn)
    _draw_panel(ax_mc, output_mc,   is_data=False, axis='msd', bins=msd_bins, plateau_mask_fn=plateau_fn)
    for ax, title in [
        (ax_d,  f"Triggers ($\\mu$ ref, jet $p_T \\geq {pt_min_plateau:.0f}$)"),
        (ax_mc, f"Triggers (QCD MC, jet $p_T \\geq {pt_min_plateau:.0f}$)")
    ]:
        ax.set_xlabel("Leading jet $m_{{SD}}$ [GeV]", fontsize=13)
        ax.set_ylabel("Trigger Efficiency", fontsize=13)
        ax.set_xlim(msd_bins[0], msd_bins[-1])
        ax.set_title(title, fontsize=11)
    plt.tight_layout()
    save_path = os.path.join(save_dir, "efficiency_vs_msd.png")
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    if show_img:
        plt.show()
    else:
        print(f"Saved {save_path}")
    plt.close(fig)


def generate_per_trigger_efficiency_csv(output_data, output_mc, trig_vars_2d, triggers, save_dir=None, year='2022', pt_min=450, pt_max=1000, msd_min=0, msd_max=300):
    """
    Generate CSV files with per-trigger efficiencies for data and MC.
    Applies kinematic cuts (pT >= pt_min GeV) for efficiency calculation as per AN methodology.

    Parameters:
    - output_data: dict - Output from ScaleFactor2DProcessor for data
    - output_mc: dict - Output from ScaleFactor2DProcessor for MC
    - trig_vars_2d: dict - Dictionary of 2D trigger variables
    - triggers: list - List of trigger names for this period
    - save_dir: str, optional - Directory to save CSV files
    - year: str - Year/period label
    - pt_min: float - Minimum pT for efficiency calculation (default 450 GeV)
    - pt_max: float - Maximum pT for efficiency calculation (default 1000 GeV)
    - msd_min: float - Minimum mSD for efficiency calculation (default 0 GeV)
    - msd_max: float - Maximum mSD for efficiency calculation (default 300 GeV)
    """
    if save_dir is None:
        save_dir = "./output"
    os.makedirs(save_dir, exist_ok=True)

    import csv

    # Process pt_vs_msd variable
    var_name = 'pt_vs_msd'
    if var_name not in trig_vars_2d:
        print(f"Warning: {var_name} not found in trig_vars_2d")
        return

    # Get baseline data for data and MC
    baseline_x_data = np.array(output_data['Baseline'][var_name]['x'])  # mSD
    baseline_y_data = np.array(output_data['Baseline'][var_name]['y'])  # pT
    baseline_x_mc = np.array(output_mc['Baseline'][var_name]['x'])
    baseline_y_mc = np.array(output_mc['Baseline'][var_name]['y'])

    # Apply kinematic cuts for efficiency calculation (pT >= pt_min, mSD range)
    baseline_cut_data = (baseline_y_data >= pt_min) & (baseline_y_data < pt_max) & \
                        (baseline_x_data >= msd_min) & (baseline_x_data < msd_max)
    baseline_cut_mc = (baseline_y_mc >= pt_min) & (baseline_y_mc < pt_max) & \
                      (baseline_x_mc >= msd_min) & (baseline_x_mc < msd_max)

    n_baseline_data = np.sum(baseline_cut_data)
    n_baseline_mc = np.sum(baseline_cut_mc)

    print(f"\n{'='*60}")
    print(f"Per-Trigger Efficiency Calculation for {year}")
    print(f"Kinematic region: {pt_min} ≤ pT < {pt_max} GeV, {msd_min} ≤ mSD < {msd_max} GeV")
    print(f"Baseline events - Data: {n_baseline_data}, MC: {n_baseline_mc}")
    print(f"{'='*60}\n")

    # Prepare CSV data
    csv_data = []
    csv_data.append(['Trigger', 'Data_Efficiency', 'Data_Pass', 'Data_Total',
                     'MC_Efficiency', 'MC_Pass', 'MC_Total', 'Scale_Factor'])

    # Calculate efficiency for each trigger
    for trigger in triggers:
        # Get per-trigger passing events
        if trigger not in output_data['PerTrigger']:
            print(f"Warning: Trigger {trigger} not found in data output")
            continue

        pass_x_data = np.array(output_data['PerTrigger'][trigger][var_name]['pass_x'])
        pass_y_data = np.array(output_data['PerTrigger'][trigger][var_name]['pass_y'])
        pass_x_mc = np.array(output_mc['PerTrigger'][trigger][var_name]['pass_x'])
        pass_y_mc = np.array(output_mc['PerTrigger'][trigger][var_name]['pass_y'])

        # Apply kinematic cuts
        pass_cut_data = (pass_y_data >= pt_min) & (pass_y_data < pt_max) & \
                        (pass_x_data >= msd_min) & (pass_x_data < msd_max)
        pass_cut_mc = (pass_y_mc >= pt_min) & (pass_y_mc < pt_max) & \
                      (pass_x_mc >= msd_min) & (pass_x_mc < msd_max)

        n_pass_data = np.sum(pass_cut_data)
        n_pass_mc = np.sum(pass_cut_mc)

        # Calculate efficiencies
        eff_data = n_pass_data / n_baseline_data if n_baseline_data > 0 else 0.0
        eff_mc = n_pass_mc / n_baseline_mc if n_baseline_mc > 0 else 0.0
        sf = eff_data / eff_mc if eff_mc > 0 else 1.0

        # Simplify trigger name for display
        trigger_display = trigger.replace('HLT_', '')

        csv_data.append([
            trigger_display,
            f"{eff_data:.4f}",
            f"{n_pass_data}",
            f"{n_baseline_data}",
            f"{eff_mc:.4f}",
            f"{n_pass_mc}",
            f"{n_baseline_mc}",
            f"{sf:.4f}"
        ])

        print(f"{trigger_display:60s} | Data: {eff_data:.2%} ({n_pass_data}/{n_baseline_data}) | MC: {eff_mc:.2%} ({n_pass_mc}/{n_baseline_mc}) | SF: {sf:.4f}")

    # Calculate combined (OR) efficiency
    pass_x_data = np.array(output_data['Numerator'][var_name]['pass_x'])
    pass_y_data = np.array(output_data['Numerator'][var_name]['pass_y'])
    pass_x_mc = np.array(output_mc['Numerator'][var_name]['pass_x'])
    pass_y_mc = np.array(output_mc['Numerator'][var_name]['pass_y'])

    pass_cut_data = (pass_y_data >= pt_min) & (pass_y_data < pt_max) & \
                    (pass_x_data >= msd_min) & (pass_x_data < msd_max)
    pass_cut_mc = (pass_y_mc >= pt_min) & (pass_y_mc < pt_max) & \
                  (pass_x_mc >= msd_min) & (pass_x_mc < msd_max)

    n_pass_data = np.sum(pass_cut_data)
    n_pass_mc = np.sum(pass_cut_mc)

    eff_data = n_pass_data / n_baseline_data if n_baseline_data > 0 else 0.0
    eff_mc = n_pass_mc / n_baseline_mc if n_baseline_mc > 0 else 0.0
    sf = eff_data / eff_mc if eff_mc > 0 else 1.0

    csv_data.append([
        'Combined (OR)',
        f"{eff_data:.4f}",
        f"{n_pass_data}",
        f"{n_baseline_data}",
        f"{eff_mc:.4f}",
        f"{n_pass_mc}",
        f"{n_baseline_mc}",
        f"{sf:.4f}"
    ])

    print(f"{'Combined (OR)':60s} | Data: {eff_data:.2%} ({n_pass_data}/{n_baseline_data}) | MC: {eff_mc:.2%} ({n_pass_mc}/{n_baseline_mc}) | SF: {sf:.4f}")

    # Write CSV file
    csv_filename = os.path.join(save_dir, f"per_trigger_efficiency_{year}.csv")
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_data)

    print(f"\nSaved per-trigger efficiency table to: {csv_filename}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="2D Scale Factor Analysis")
    parser.add_argument('--year', required=False, choices=list(trigger_dict_periods.keys()), default='2022',
                        help="Year of the dataset (e.g., 2022, 2023, etc.).")
    parser.add_argument('--test', action='store_true',
                        help="Test mode: only process first 10 ROOT files.")
    parser.add_argument('--quick-test', action='store_true',
                        help="Quick debug mode: process only a couple of files/chunks to finish in <10 min.")
    parser.add_argument('--quick-files-per-dataset', type=int, default=2,
                        help="Number of files per dataset in --quick-test mode (default: 2).")
    parser.add_argument('--quick-workers', type=int, default=1,
                        help="Number of workers in --quick-test mode (default: 1).")
    parser.add_argument('--quick-maxchunks', type=int, default=2,
                        help="Maximum chunks per dataset in --quick-test mode (default: 2).")
    parser.add_argument('--chunksize', type=int, default=100000,
                        help="Events per chunk for coffea Runner (default: 100000).")
    parser.add_argument('--txbb-region', choices=['pass', 'fail', 'inclusive', 'all'], default='all',
                        help="txbb region: 'pass' (>=0.82), 'fail' (<0.82), 'inclusive' (no cut), 'all' (run pass+fail).")
    parser.add_argument('--pt-min', type=float, default=PT_MIN_AN,
                        help=f"Minimum pT [GeV] for per-trigger efficiency CSV calculation (default: {PT_MIN_AN:g}).")
    parser.add_argument('--pt-max', type=float, default=PT_MAX_AN,
                        help=f"Maximum pT [GeV] for per-trigger efficiency CSV calculation (default: {PT_MAX_AN:g}).")
    parser.add_argument('--msd-min', type=float, default=MSD_MIN_AN,
                        help=f"Minimum mSD [GeV] for per-trigger efficiency CSV calculation (default: {MSD_MIN_AN:g}).")
    parser.add_argument('--msd-max', type=float, default=MSD_MAX_AN,
                        help=f"Maximum mSD [GeV] for per-trigger efficiency CSV calculation (default: {MSD_MAX_AN:g}).")
    parser.add_argument('--disable-jetid-correction', action='store_true',
                        help="Disable correctionlib-based JetID recomputation and use NanoAOD isTight directly.")
    parser.add_argument('--outdir', default='.',
                        help="Base output directory. Results are saved under <outdir>/figures_sf_2d and <outdir>/output.")
    args = parser.parse_args()

    year = args.year
    test_mode = args.test
    quick_test_mode = args.quick_test
    pt_min_csv = args.pt_min
    pt_max_csv = args.pt_max
    msd_min_csv = args.msd_min
    msd_max_csv = args.msd_max
    use_jetid_correction = not args.disable_jetid_correction
    base_outdir = os.path.abspath(args.outdir)
    figures_base_dir = os.path.join(base_outdir, "figures_sf_2d")
    output_base_dir = os.path.join(base_outdir, "output")

    # Create output directories early (required for condor file transfer)
    os.makedirs(figures_base_dir, exist_ok=True)
    os.makedirs(output_base_dir, exist_ok=True)

    if test_mode:
        print(f'Running in TEST mode: only processing first 10 files')
    if quick_test_mode:
        print(
            f"Running in QUICK TEST mode: "
            f"{args.quick_files_per_dataset} files/dataset, "
            f"maxchunks={args.quick_maxchunks}, workers={args.quick_workers}"
        )
    print(f'Processing period {year}')
    print(f"JetID correction: {'enabled' if use_jetid_correction else 'disabled'}")
    print(
        "Efficiency/SF kinematic region: "
        f"{pt_min_csv:g} <= pT < {pt_max_csv:g} GeV, "
        f"{msd_min_csv:g} <= mSD < {msd_max_csv:g} GeV"
    )
    print(f"Output base directory: {base_outdir}")

    # File dictionaries
    file_dict_periods = {
        '2022': {'ttbar': "TTtoLNu2Q", 'MuonData': "Muon_Run2022C"},
        '2022EE': {'ttbar': "TTtoLNu2Q", 'MuonData': "Muon_Run2022E"},
        '2023': {'ttbar': "TTtoLNu2Q", 'MuonData': [
            "Muon0_Run2023C-v1", "Muon0_Run2023C-v2", "Muon0_Run2023C-v3", "Muon0_Run2023C-v4",
            "Muon1_Run2023C-v1", "Muon1_Run2023C-v2", "Muon1_Run2023C-v3", "Muon1_Run2023C-v4"
        ]},
        '2023BPix': {'ttbar': "TTtoLNu2Q", 'MuonData': ["Muon0_Run2023D-v1", "Muon0_Run2023D-v2"]}
    }

    prod_modes_map = file_dict_periods[year]

    # Trigger SFs are measured inclusively (not separated by production mode)
    prod_modes = ['Inclusive']

    # txbb pass/fail regions
    if args.txbb_region == 'all':
        txbb_regions = ['pass', 'fail']
    else:
        txbb_regions = [args.txbb_region]

    for prod_mode in prod_modes:
        for txbb_region in txbb_regions:
            print(f"\n{'='*60}")
            print(f"Measuring trigger SF: {prod_mode}, txbb={txbb_region} (WP={TXBB_WP})")
            print(f"{'='*60}\n")

        # Load fileset once per period (shared across txbb regions)
        fileset = {}

        for dataset_type, file_key in prod_modes_map.items():
            json_path = os.path.join(f'infiles/{year}', f"{year}_{dataset_type}.json")
            print(f"Reading: {json_path} (key: {file_key})")

            if not os.path.exists(json_path):
                print(f"  Warning: {json_path} not found. Skipping {dataset_type}.")
                continue

            with open(json_path, 'r') as file:
                data = json.load(file)

            # file_key can be a string or list of strings (to combine multiple run versions)
            if isinstance(file_key, list):
                samples_by_key = []
                for k in file_key:
                    s = data.get(k, [])
                    samples_by_key.append(s)
                    print(f"  + {k}: {len(s)} files")
                samples = [x for sub in samples_by_key for x in sub]
            else:
                samples = data.get(file_key, [])

            if quick_test_mode:
                if isinstance(file_key, list):
                    samples = round_robin_sample_lists(samples_by_key, args.quick_files_per_dataset)
                    print(
                        "QUICK TEST MODE: Round-robin sampled "
                        f"{len(samples)} files across {len(file_key)} keys"
                    )
                else:
                    samples = samples[:args.quick_files_per_dataset]
                    print(f"QUICK TEST MODE: Limited to {len(samples)} files")
            elif test_mode and isinstance(file_key, list):
                samples = round_robin_sample_lists(samples_by_key, 10)
                print(
                    "TEST MODE: Round-robin sampled "
                    f"{len(samples)} files across {len(file_key)} keys"
                )
            elif test_mode:
                samples = samples[:10]
                print(f"TEST MODE: Limited to {len(samples)} files")

            dataset_key = f"{dataset_type}_{year}"
            fileset[dataset_key] = samples
            print(f"Loaded {len(samples)} files for dataset '{dataset_key}'")

        if len(fileset) == 0:
            print(f"Error: No data loaded for production mode {prod_mode}")
            continue

        print("Constructed fileset:")
        print(json.dumps(fileset, indent=2))

        workers = args.quick_workers if quick_test_mode else 2
        runner_maxchunks = args.quick_maxchunks if quick_test_mode else None
        iterative_run = processor.Runner(
            executor=processor.FuturesExecutor(compression=None, workers=workers),
            schema=NanoAODSchema,
            skipbadfiles=True,
            savemetrics=True,
            chunksize=args.chunksize,
            maxchunks=runner_maxchunks,
        )

        # --- Inner loop: run processor separately for each txbb region ---
        for txbb_region in txbb_regions:
            print(f"\n{'-'*60}")
            print(f"  txbb region: {txbb_region} (WP = {TXBB_WP})")
            print(f"{'-'*60}")

            outputs = {}

            for dataset_key in fileset.keys():
                if len(fileset[dataset_key]) == 0:
                    print(f"\nSkipping {dataset_key} - no files found")
                    continue

                print(f"\nRunning 2D SF processor for {dataset_key} [{txbb_region}]...")
                temp_fileset = {dataset_key: fileset[dataset_key]}

                out = iterative_run(
                    temp_fileset,
                    treename="Events",
                    processor_instance=ScaleFactor2DProcessor(
                        year, trig_vars_2d,
                        baseline_key=prod_mode,
                        txbb_region=txbb_region,
                        use_jetid_correction=use_jetid_correction
                    ),
                )
                outputs[dataset_key] = out[0]

            data_key = [k for k in outputs.keys() if 'Muon' in k][0] if any('Muon' in k for k in outputs.keys()) else None
            mc_key   = [k for k in outputs.keys() if 'ttbar' in k][0] if any('ttbar' in k for k in outputs.keys()) else None

            if data_key is None or mc_key is None:
                print("Error: Missing data or MC output. Skipping plots.")
                continue

            output_data = outputs[data_key]
            output_mc   = outputs[mc_key]

            # Output paths include txbb_region
            FIGURES_DIR = os.path.join(figures_base_dir, year, prod_mode, txbb_region)
            os.makedirs(FIGURES_DIR, exist_ok=True)

            output_path_data = os.path.join(output_base_dir, f"sf_2d_{year}_{prod_mode}_{txbb_region}_data.coffea")
            output_path_mc   = os.path.join(output_base_dir, f"sf_2d_{year}_{prod_mode}_{txbb_region}_mc.coffea")
            os.makedirs(os.path.dirname(output_path_data), exist_ok=True)

            util.save(output_data, output_path_data)
            util.save(output_mc,   output_path_mc)
            print(f"Saved data output to {output_path_data}")
            print(f"Saved MC output   to {output_path_mc}")

            plot_2d_scale_factors(
                output_data, output_mc, trig_vars_2d,
                save_dir=FIGURES_DIR, year=year, show_img=False
            )

            plot_1d_efficiency_projections(
                output_data, output_mc, trig_vars_2d,
                triggers=trigger_dict_periods[year],
                save_dir=FIGURES_DIR,
                year=year,
                show_img=False,
                pt_min_plateau=pt_min_csv,
                msd_min_plot=msd_min_csv,
            )
            print(f"\nAll plots saved to {FIGURES_DIR}")

            csv_output_dir = os.path.join(output_base_dir, year, prod_mode, txbb_region)
            os.makedirs(csv_output_dir, exist_ok=True)

            generate_per_trigger_efficiency_csv(
                output_data, output_mc, trig_vars_2d,
                triggers=trigger_dict_periods[year],
                save_dir=csv_output_dir,
                year=f"{year}_{prod_mode}_{txbb_region}",
                pt_min=pt_min_csv,
                pt_max=pt_max_csv,
                msd_min=msd_min_csv,
                msd_max=msd_max_csv
            )

    print(f"\n{'='*60}")
    print(f"Completed processing for all production modes")
    print(f"{'='*60}")
