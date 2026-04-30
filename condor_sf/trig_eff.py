import warnings
warnings.filterwarnings('ignore')
import os
import awkward as ak
import uproot
import hist
import numpy as np
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
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

#line thickness
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 5
import itertools
import json


# allow importing module from the sibling directory ../trig_eff
from pathlib import Path
import sys


def compute_mjj(events):
    """Dijet invariant mass of the two leading AK4 jets outside the candidate AK8 jet.
    Returns NaN for events without a valid pair (e.g. ggF)."""
    fatjets = events.FatJet
    candidatejet = fatjets[
        (fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjets.isTight
    ]
    candidatejet = candidatejet[:, :2]
    candidatejet = ak.firsts(
        candidatejet[ak.argmax(candidatejet.particleNet_XbbVsQCD, axis=1, keepdims=True)]
    )

    jets = events.Jet
    jets = jets[(jets.pt > 30.) & (abs(jets.eta) < 5.0) & jets.isTight]
    jets = jets[:, :4]

    dR = jets.delta_r(candidatejet)
    ak4_outside_ak8 = jets[dR > 0.8]

    jet1 = ak.firsts(ak4_outside_ak8[:, 0:1])
    jet2 = ak.firsts(ak4_outside_ak8[:, 1:2])

    return ak.fill_none((jet1 + jet2).mass, np.nan)


trig_vars = {
    'ht': {
        'label': "H_{T} [GeV]",
        'axis': hist.axis.Regular(bins=100, start=0, stop=2000, name="ht", label="H_{T} [GeV]"),
        'proc': lambda events: ak.sum(events.FatJet.pt, axis=1)
    },
    'pt': {
        'label': "Leading Jet $p_{T}$ [GeV]",
        'axis': hist.axis.Regular(bins=30, start=0, stop=1200, name="pt", label="Leading Jet $p_{T}$ [GeV]"),
        'proc': lambda events: ak.fill_none(ak.pad_none(events.FatJet.pt, 1, clip=True)[:, 0], np.nan)
    },
    'msd': {
        'label': "Leading Jet $m_{SD}$ [GeV]",
        'axis': hist.axis.Regular(bins=15, start=0, stop=300, name="msd", label="Leading Jet $m_{SD}$ [GeV]"),
        'proc': lambda events: ak.fill_none(ak.pad_none(events.FatJet.msoftdrop, 1, clip=True)[:, 0], np.nan)
    },
    'num_ak4': {
        'label': "Number of AK4 Jets",
        'axis': hist.axis.Integer(0, 20, name="num_ak4", label="Number of AK4 Jets"),
        'proc': lambda events: ak.num(events.Jet)
    },
    'gen_H_pt': {
        'label': "Gen Higgs pT [GeV]",
        'axis': hist.axis.Regular(bins=30, start=0, stop=1200, name="gen_H_pt", label="Gen Higgs pT [GeV]"),
        'proc': lambda events: events.HTXS.Higgs_pt
    },
    'particleNet_XbbVsQCD':{
        'label': "Leading Particle Net TXbb score",
        'axis': hist.axis.Regular(bins=30, start=0, stop=1, name="particleNet_XbbVsQCD", label="Leading Particle Net TXbb score"),
        'proc': lambda events: ak.fill_none(
            ak.pad_none(events.FatJet.particleNet_XbbVsQCD, 1, clip=True)[:, 0],
            np.nan
        )
    },
    'mjj': {
        'label': "$m_{jj}$ [GeV]",
        'axis': hist.axis.Regular(bins=20, start=1000, stop=3000, name="mjj", label="$m_{jj}$ [GeV]"),
        'proc': compute_mjj
    }
}

trigger_dict_periods = {
    '2022': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        # 'AK8PFJet425_SoftDropMass40',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65'
    ],
    '2022EE': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        # 'AK8PFJet425_SoftDropMass40',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65'
    ],
    '2023': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        # 'AK8PFJet425_SoftDropMass40',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'VBF_DiPFJet125_45_Mjj720_Detajj3p0'
    ],
    '2023BPix': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        # 'AK8PFJet425_SoftDropMass40',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'VBF_DiPFJet125_45_Mjj720_Detajj3p0'
    ]
}

# Which variables to plot for each individual trigger
trigger_var_map = {
    'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35': ['pt', 'msd'],
    'AK8PFJet250_SoftDropMass40_PNetBB0p06':             ['pt', 'msd'],
    'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65': ['ht', 'particleNet_XbbVsQCD'],
    'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70':      ['ht', 'particleNet_XbbVsQCD'],
    'VBF_DiPFJet125_45_Mjj720_Detajj3p0':                  ['mjj'],
}

def create_baseline_selection_mask_new(events, tag):
    n_events = len(events)
    # Start with all events selected.
    baseline_mask = np.ones(n_events, dtype=bool)
    
    # ------------------------------------------------------------------
    # Step 1: Common jet selection, e.g. require at least one AK4 jet with pt > 15 GeV.
    jets_all = events.Jet[events.Jet.pt > 15]
    has_jets = ak.num(jets_all) >= 1
    baseline_mask &= ak.to_numpy(has_jets)
    
    # ------------------------------------------------------------------
    # Step 2: General selection on FatJets.
    # Select events that have at least one FatJet satisfying:
    #   - pt > 250,
    #   - |eta| < 2.5, and
    #   - particleNet_XbbVsQCD > 0.4.
    fatjets = events.FatJet
    fatjet_mask = (fatjets.pt > 250) & (abs(fatjets.eta) < 2.5) & (fatjets.particleNet_XbbVsQCD > 0.4)
    has_fatjet = ak.sum(fatjet_mask, axis=1) >= 1
    baseline_mask &= ak.to_numpy(has_fatjet)

    candidatejet = fatjets[
            (fatjets.pt > 200)
            & (abs(fatjets.eta) < 2.5)
            & fatjets.isTight 
        ]
    candidatejet = candidatejet[:, :2]
    candidatejet = ak.firsts(candidatejet[ak.argmax(candidatejet.particleNet_XbbVsQCD, axis=1, keepdims=True)])

    
    # ------------------------------------------------------------------
    # Step 3: If tag is "VBF", apply additional VBF selection on top of general selection.
    if tag == "VBF":
        # Select additional jets (AK4) with:
        #   - pt > 15 GeV (common threshold) and
        #   - |eta| < 5.0.
        jets = events.Jet
        jets = jets[
            (jets.pt > 30.)
            & (abs(jets.eta) < 5.0)
            & jets.isTight
        ]

        # only consider first 4 jets to be consistent with old framework
        jets = jets[:, :4]
        
        # Require at least two such jets per event.
        has_two_jets = ak.num(jets) >= 2
        
        # VBF specific variables                                                      
        dR = jets.delta_r(candidatejet)
        ak4_outside_ak8 = jets[dR > 0.8]

        jet1 = ak4_outside_ak8[:, 0:1]
        jet2 = ak4_outside_ak8[:, 1:2]

        deta = abs(ak.firsts(jet1).eta - ak.firsts(jet2).eta)
        mjj = ( ak.firsts(jet1) + ak.firsts(jet2) ).mass

        has_valid_pair = ((deta > 3.5) & (mjj > 1000))
        
        # # Form all combinations of two jets in each event.
        # jet_pairs = ak.combinations(additional_jets, 2, fields=["jet1", "jet2"])
        # # Calculate the absolute difference in eta (rapidity separation) of the two jets.
        # delta_eta = abs(jet_pairs.jet1.eta - jet_pairs.jet2.eta)
        # # Calculate the dijet invariant mass for the jet pair.
        # mjj = (jet_pairs.jet1 + jet_pairs.jet2).mass
        # # VBF criteria: at least one pair must satisfy:
        # #   - |delta_eta| > 3.5, and
        # #   - mjj > 1000 GeV.
        # has_valid_pair = ak.any((deta > 3.5) & (mjj > 1000), axis=1)
        
        # Combine the VBF-specific requirements.
        vbf_mask = has_two_jets & has_valid_pair
        
        # Update the cumulative baseline mask.
        baseline_mask &= ak.to_numpy(vbf_mask)
    
    # Final sanity check: baseline_mask length must equal the number of original events.
    assert len(baseline_mask) == n_events, "Baseline mask length does not match number of events!"
    
    return baseline_mask

class ParkingSoupProcessor(processor.ProcessorABC):
    """
    Calculate histograms for two scenarios:
      1) Baseline (no triggers)
      2) Original Trigger Soup
    """
    def __init__(self, orig_trig_soup, trig_vars, baseline_key='VBF'):
        self.orig_trig_soup = orig_trig_soup
        self.trig_vars = trig_vars
        self.baseline_key = baseline_key

        # Accumulators: Baseline (denominator) + OR + one key per individual trigger
        trigger_keys = ['OR'] + orig_trig_soup

        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                var_name: list_accumulator() for var_name in trig_vars
            }),
            **{
                trig_key: dict_accumulator({
                    var_name: dict_accumulator({'pass': list_accumulator()})
                    for var_name in trig_vars
                })
                for trig_key in trigger_keys
            }
        })

    def process(self, events):
        output = self.output

        # 1) Baseline mask
        baseline_mask = create_baseline_selection_mask_new(events, tag=self.baseline_key)

        # 2) Variable arrays
        variables = {}
        for var_name, var_info in self.trig_vars.items():
            variables[var_name] = var_info['proc'](events)

        # 3) Trigger masks: OR + each individual trigger
        if len(self.orig_trig_soup) > 0:
            combined_or_mask = ak.zeros_like(events.HLT[self.orig_trig_soup[0]], dtype=bool)
            for trig in self.orig_trig_soup:
                combined_or_mask = combined_or_mask | events.HLT[trig]
        else:
            combined_or_mask = ak.zeros_like(baseline_mask, dtype=bool)

        indiv_masks = {trig: ak.to_numpy(events.HLT[trig]) for trig in self.orig_trig_soup}
        or_mask_np  = ak.to_numpy(combined_or_mask)

        # 4) Fill arrays per variable, each with its own NaN-validity mask.
        #    This keeps variables like mjj (NaN for ggF) from zeroing out
        #    the denominator of unrelated variables.
        for var_name in self.trig_vars:
            var_array = variables[var_name]
            valid = ~np.isnan(ak.to_numpy(var_array))

            sel_base = baseline_mask & valid
            sel_or   = sel_base & or_mask_np

            output['Baseline'][var_name].extend(ak.to_numpy(var_array[sel_base]).tolist())
            output['OR'][var_name]['pass'].extend(ak.to_numpy(var_array[sel_or]).tolist())

            for trig in self.orig_trig_soup:
                sel_trig = sel_base & indiv_masks[trig]
                output[trig][var_name]['pass'].extend(ak.to_numpy(var_array[sel_trig]).tolist())

        return output

    def postprocess(self, accumulator):
        return accumulator

def plot_1d_trigger_soup_cms(output, trig_vars, save_dir=None, year=2022, data=False, show_img=True):
    """
    Plot the trigger efficiencies for each trigger combination.

    Parameters:
    - output: dict
        The output dictionary from the TriggerSoupProcessor.
    - trig_vars: dict
        Dictionary of trigger variables with their processing functions and axis information.
    - save_dir: str, optional
        Directory to save the figures. If None, defaults to "./figures".
    - year: int, optional
        Year label for the plots.
    - data: bool, optiona
        plotting with data or MC
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import mplhep as hep

    if save_dir is None:
        save_dir = "./figures"
    os.makedirs(save_dir, exist_ok=True)

    # Update font size for plots
    plt.rcParams.update({"font.size": 14})

    # Define colors for plotting
    colors = plt.cm.tab10.colors

    # Extract the (single) combination name from the output, excluding 'Baseline'
    combination_names = [key for key in output.keys() if key != 'Baseline']

    for var_name in trig_vars.keys():
        # Prepare bin edges for the variable
        bin_edges = trig_vars[var_name]['axis'].edges

        # Create figure with two subplots
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(12, 14),
            gridspec_kw=dict(height_ratios=[4, 1], hspace=0.07),
            sharex=True
        )

        # Initialize histograms and efficiencies dictionaries
        hists = {}
        ratios = {}
        total_efficiencies = {}

        # Get 'total' values from 'Baseline'
        total_values = np.array(output['Baseline'][var_name])
        weights_total = np.ones_like(total_values)  # Modify if you have actual weights

        # 'All' histogram (denominator for efficiency)
        hists['All'], _ = np.histogram(total_values, bins=bin_edges, weights=weights_total)

        # Plot 'Baseline' histogram (no efficiency shown for baseline)
        hep.histplot(
            hists['All'],
            bins=bin_edges,
            yerr=False,
            label='Baseline',
            ax=ax,
            color='black',
            linestyle='dashed',
            linewidth=2
        )

        max_efficiency = 0.0

        # Only plot the single combination present in output (if any)
        if len(combination_names) == 0:
            print(f"No trigger combinations found in output for variable {var_name}; skipping.")
            plt.close(fig)
            continue

        combination = combination_names[0]
        # compute pass histogram for the single combination
        pass_values = np.array(output[combination][var_name]['pass'])
        weights_pass = np.ones_like(pass_values)
        hists[combination], _ = np.histogram(pass_values, bins=bin_edges, weights=weights_pass)

        # ratio (efficiency) = pass / all
        ratios[combination] = np.divide(
            hists[combination], hists['All'],
            out=np.zeros_like(hists[combination], dtype=float),
            where=hists['All'] != 0
        )

        # total efficiency for this combination
        total_events = np.sum(hists['All'])
        passed_events = np.sum(hists[combination])
        total_efficiency = passed_events / total_events if total_events > 0 else 0.0
        total_efficiencies[combination] = total_efficiency
        max_efficiency = total_efficiency

        # Label the combination as 'All_triggers' in the plot legend
        label_str = f"{combination} ({total_efficiency:.2%})"

        # Plot the single combination on top of baseline
        hep.histplot(
            hists[combination],
            bins=bin_edges,
            yerr=False,
            label=label_str,
            ax=ax,
            color='red',
            linestyle='-',
            linewidth=2,
            zorder=5
        )

        # Plot ratio for this single combination
        hep.histplot(
            (ratios[combination], bin_edges),
            yerr=False,
            ax=rax,
            histtype='step',
            color='red',
            linewidth=2,
        )

        # Set labels and titles
        title = f"{var_name} (Max Efficiency: {max_efficiency:.2%})"
        ax.set_title(title)

        ax.set_ylabel("Events [A.U.]")
        ax.legend()

        rax.legend()
        rax.grid(axis='y')
        rax.set_xlabel(trig_vars[var_name]['label'])
        rax.set_ylabel("Efficiency")
        # print("Baseline histogram counts:", hists['All'])
        # print("'All_triggers' histogram counts:", hists['All_triggers'])

        # Add CMS label (modify according to your style)
        hep.cms.label(ax=ax, data=data, year=year, com="13.6")

        # Save the figure
        plt.savefig(f"{save_dir}/{var_name}_trigger_soup.png", dpi=200, bbox_inches='tight')
        if show_img:
            plt.show()
        else:
            print(f"Saved {save_dir}/{var_name}_trigger_soup.png")
        plt.close(fig)

if __name__ == "__main__":
    import argparse

    # --- Parse command line arguments ---
    parser = argparse.ArgumentParser(description="Trigger Efficiency Analysis")
    parser.add_argument('--year', required=False, choices=list(trigger_dict_periods.keys()), default='2022',
                        help="Year of the dataset (e.g., 2022, 2023, etc.).")
    args = parser.parse_args()
    
    year = args.year
    print(f'Processing period {year}')
    all_triggers = trigger_dict_periods[year]
    file_dict = {'VBF': "VBF_Hto2B", 'ggF': "GluGlu_Hto2B"}
    outputs = {key_prod: [] for key_prod in file_dict.keys()}
    for prod_mode, file_name in file_dict.items():
        print("Currently processing production mode", prod_mode)
        fileset = {}
        json_path = os.path.join(f'infiles/{year}', f"{year}_{prod_mode}.json")
        print(f"Reading: {json_path} (key: {file_name})")
        with open(json_path, 'r') as file:
            data = json.load(file)
        samples = data.get(file_name, [])
        dataset_key = f"{prod_mode}_{year}"
        # Keep previous behavior of using a small subset for quick tests
        fileset[dataset_key] = samples

        print("Constructed fileset:")
        print(json.dumps(fileset, indent=2))
        iterative_run = processor.Runner(
                executor = processor.FuturesExecutor(compression=None, workers=2),
                schema=NanoAODSchema,
                skipbadfiles=True,  # Skip files that fail to open
                savemetrics=True,   # Save metrics to understand where failures occur
        )
        out = iterative_run(
            fileset,
            treename="Events",
            processor_instance=ParkingSoupProcessor(all_triggers, trig_vars, baseline_key=prod_mode),
        )
        output = out[0]
        outputs[prod_mode] = output
    
        FIGURES_DIR = os.path.join(os.getcwd(), "figures/", f"{year}/{prod_mode}/")
        os.makedirs(FIGURES_DIR, exist_ok=True)
        output_path = os.path.join(os.getcwd(), f"output/trig_soup_eff_{year}_{prod_mode}.coffea")
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        util.save(output, output_path)
        plot_1d_trigger_soup_cms(output, trig_vars, save_dir=FIGURES_DIR, year=year, data=False, show_img=False)

        # Per-trigger plots: each trigger plotted against only its assigned variables
        for trig in all_triggers:
            vars_for_trig = trigger_var_map.get(trig, [])
            if not vars_for_trig:
                continue
            mini_output    = {'Baseline': output['Baseline'], trig: output[trig]}
            mini_trig_vars = {v: trig_vars[v] for v in vars_for_trig}
            trig_fig_dir   = os.path.join(FIGURES_DIR, trig)
            plot_1d_trigger_soup_cms(mini_output, mini_trig_vars, save_dir=trig_fig_dir,
                                     year=year, data=False, show_img=False)
