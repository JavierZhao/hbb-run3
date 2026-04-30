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

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 5
import itertools
import json


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
        'VBF_DiPFJet125_45_Mjj720_Detajj3p0'
    ],
    '2023BPix': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'VBF_DiPFJet125_45_Mjj720_Detajj3p0'
    ]
}

def create_baseline_selection_mask_new(events, tag):
    n_events = len(events)
    baseline_mask = np.ones(n_events, dtype=bool)
    cutflow = {"All events": n_events}

    # Step 1: Common jet selection
    jets_all = events.Jet[events.Jet.pt > 15]
    has_jets = ak.num(jets_all) >= 1
    baseline_mask &= ak.to_numpy(has_jets)
    cutflow[">=1 AK4 jet (pT > 15 GeV)"] = int(np.sum(baseline_mask))

    # Step 2: General selection on FatJets
    fatjets = events.FatJet
    fatjet_mask = (fatjets.pt > 250) & (abs(fatjets.eta) < 2.5) & (fatjets.particleNet_XbbVsQCD > 0.4)
    has_fatjet = ak.sum(fatjet_mask, axis=1) >= 1
    baseline_mask &= ak.to_numpy(has_fatjet)
    cutflow[">=1 FatJet (pT > 250, |eta| < 2.5, Xbb > 0.4)"] = int(np.sum(baseline_mask))

    candidatejet = fatjets[
            (fatjets.pt > 200)
            & (abs(fatjets.eta) < 2.5)
            & fatjets.isTight
        ]
    candidatejet = candidatejet[:, :2]
    candidatejet = ak.firsts(candidatejet[ak.argmax(candidatejet.particleNet_XbbVsQCD, axis=1, keepdims=True)])

    # Step 3: If tag is "VBF", apply additional VBF selection
    if tag == "VBF":
        jets = events.Jet
        jets = jets[
            (jets.pt > 30.)
            & (abs(jets.eta) < 5.0)
            & jets.isTight
        ]
        jets = jets[:, :4]
        has_two_jets = ak.num(jets) >= 2
        baseline_mask &= ak.to_numpy(has_two_jets)
        cutflow[">=2 tight AK4 jets (pT > 30, |eta| < 5)"] = int(np.sum(baseline_mask))

        dR = jets.delta_r(candidatejet)
        ak4_outside_ak8 = jets[dR > 0.8]

        jet1 = ak4_outside_ak8[:, 0:1]
        jet2 = ak4_outside_ak8[:, 1:2]

        deta = abs(ak.firsts(jet1).eta - ak.firsts(jet2).eta)
        mjj = ( ak.firsts(jet1) + ak.firsts(jet2) ).mass

        has_valid_pair = ((deta > 3.5) & (mjj > 1000))
        baseline_mask &= ak.to_numpy(has_valid_pair)
        cutflow["VBF pair (dR > 0.8, deta > 3.5, mjj > 1000)"] = int(np.sum(baseline_mask))

    assert len(baseline_mask) == n_events, "Baseline mask length does not match number of events!"
    return baseline_mask, cutflow


class TriggerEfficiency2DProcessor(processor.ProcessorABC):
    """
    Calculate 2D trigger efficiencies (leading jet pT vs mass) for individual triggers and their OR.
    """
    def __init__(self, triggers, trig_vars_2d, baseline_key='VBF'):
        self.triggers = triggers
        self.trig_vars_2d = trig_vars_2d
        self.baseline_key = baseline_key

        # Add 'OR' to the list of triggers
        self.trigger_list = triggers + ['OR']

        # Define cutflow steps (must match order in create_baseline_selection_mask_new)
        self.cutflow_steps = [
            "All events",
            ">=1 AK4 jet (pT > 15 GeV)",
            ">=1 FatJet (pT > 250, |eta| < 2.5, Xbb > 0.4)",
        ]
        if baseline_key == "VBF":
            self.cutflow_steps += [
                ">=2 tight AK4 jets (pT > 30, |eta| < 5)",
                "VBF pair (dR > 0.8, deta > 3.5, mjj > 1000)",
            ]
        self.cutflow_steps.append("Valid mSD & pT")

        # Initialize output using accumulators
        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                var_name: dict_accumulator({
                    'x': list_accumulator(),
                    'y': list_accumulator()
                })
                for var_name in trig_vars_2d
            }),
            'cutflow': dict_accumulator({
                step: list_accumulator() for step in self.cutflow_steps
            }),
            **{
                trigger: dict_accumulator({
                    var_name: dict_accumulator({
                        'pass_x': list_accumulator(),
                        'pass_y': list_accumulator()
                    }) for var_name in trig_vars_2d
                }) for trigger in self.trigger_list
            }
        })

    def process(self, events):
        output = self.output

        baseline, cutflow = create_baseline_selection_mask_new(events, tag=self.baseline_key)

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

        # Combined selection for baseline (no trigger requirement)
        selection_baseline = baseline & valid_vars_mask
        cutflow["Valid mSD & pT"] = int(np.sum(selection_baseline))

        # Accumulate cutflow counts
        for step in self.cutflow_steps:
            output['cutflow'][step].extend([cutflow[step]])

        # Fill baseline histograms
        for var_name in self.trig_vars_2d:
            var_x = variables[var_name]['x']
            var_y = variables[var_name]['y']

            var_x_baseline = var_x[selection_baseline]
            var_y_baseline = var_y[selection_baseline]

            output['Baseline'][var_name]['x'].extend(ak.to_numpy(var_x_baseline).tolist())
            output['Baseline'][var_name]['y'].extend(ak.to_numpy(var_y_baseline).tolist())

        # Process individual triggers
        for trigger in self.triggers:
            trigger_mask = events.HLT[trigger]
            selection_pass = baseline & trigger_mask & valid_vars_mask

            for var_name in self.trig_vars_2d:
                var_x = variables[var_name]['x']
                var_y = variables[var_name]['y']

                var_x_pass = var_x[selection_pass]
                var_y_pass = var_y[selection_pass]

                output[trigger][var_name]['pass_x'].extend(ak.to_numpy(var_x_pass).tolist())
                output[trigger][var_name]['pass_y'].extend(ak.to_numpy(var_y_pass).tolist())

        # Process OR of all triggers
        if len(self.triggers) > 0:
            combined_trigger_mask = ak.zeros_like(events.HLT[self.triggers[0]], dtype=bool)
            for trigger in self.triggers:
                combined_trigger_mask = combined_trigger_mask | events.HLT[trigger]
        else:
            combined_trigger_mask = ak.zeros_like(baseline, dtype=bool)

        selection_pass_or = baseline & combined_trigger_mask & valid_vars_mask

        for var_name in self.trig_vars_2d:
            var_x = variables[var_name]['x']
            var_y = variables[var_name]['y']

            var_x_pass_or = var_x[selection_pass_or]
            var_y_pass_or = var_y[selection_pass_or]

            output['OR'][var_name]['pass_x'].extend(ak.to_numpy(var_x_pass_or).tolist())
            output['OR'][var_name]['pass_y'].extend(ak.to_numpy(var_y_pass_or).tolist())

        return output

    def postprocess(self, accumulator):
        # Sum cutflow counts across partitions
        if 'cutflow' in accumulator:
            for step in accumulator['cutflow']:
                accumulator['cutflow'][step] = sum(accumulator['cutflow'][step])

        # Convert lists to arrays
        for key in accumulator:
            if key == 'Baseline':
                for var_name in accumulator[key]:
                    accumulator[key][var_name]['x'] = np.array(accumulator[key][var_name]['x'])
                    accumulator[key][var_name]['y'] = np.array(accumulator[key][var_name]['y'])
            elif key == 'cutflow':
                continue
            else:
                for var_name in accumulator[key]:
                    accumulator[key][var_name]['pass_x'] = np.array(accumulator[key][var_name]['pass_x'])
                    accumulator[key][var_name]['pass_y'] = np.array(accumulator[key][var_name]['pass_y'])
        return accumulator


def plot_2d_trigger_efficiency(output, trig_vars_2d, save_dir=None, year=2022, data=False, show_img=True):
    """
    Plot 2D trigger efficiencies for each trigger and the OR.

    Parameters:
    - output: dict
        The output dictionary from the TriggerEfficiency2DProcessor.
    - trig_vars_2d: dict
        Dictionary of 2D trigger variables.
    - save_dir: str, optional
        Directory to save the figures.
    - year: int, optional
        Year label for the plots.
    - data: bool, optional
        Plotting with data or MC.
    - show_img: bool, optional
        Whether to display the images.
    """
    if save_dir is None:
        save_dir = "./figures_2d"
    os.makedirs(save_dir, exist_ok=True)

    plt.rcParams.update({"font.size": 12})

    # Get list of triggers (excluding 'Baseline')
    trigger_names = [key for key in output.keys() if key != 'Baseline']

    for var_name, var_info in trig_vars_2d.items():
        # Get bin edges
        bins_x = var_info['axis_x'].edges
        bins_y = var_info['axis_y'].edges

        # Get baseline (total) data
        baseline_x = np.array(output['Baseline'][var_name]['x'])
        baseline_y = np.array(output['Baseline'][var_name]['y'])

        # Create 2D histogram for baseline (denominator)
        hist_baseline, xedges, yedges = np.histogram2d(
            baseline_x, baseline_y,
            bins=[bins_x, bins_y]
        )

        # Loop through each trigger
        for trigger in trigger_names:
            # Get pass data for this trigger
            pass_x = np.array(output[trigger][var_name]['pass_x'])
            pass_y = np.array(output[trigger][var_name]['pass_y'])

            # Create 2D histogram for passing events (numerator)
            hist_pass, _, _ = np.histogram2d(
                pass_x, pass_y,
                bins=[bins_x, bins_y]
            )

            # Calculate efficiency
            with np.errstate(divide='ignore', invalid='ignore'):
                efficiency = np.divide(
                    hist_pass, hist_baseline,
                    out=np.zeros_like(hist_pass, dtype=float),
                    where=hist_baseline > 0
                )

            # Calculate total efficiency
            total_pass = np.sum(hist_pass)
            total_baseline = np.sum(hist_baseline)
            total_efficiency = total_pass / total_baseline if total_baseline > 0 else 0.0

            # Create figure
            fig, ax = plt.subplots(figsize=(10, 8))

            # Plot 2D efficiency
            im = ax.pcolormesh(
                xedges, yedges, efficiency.T,
                cmap='RdYlGn',
                vmin=0, vmax=1,
                shading='flat'
            )

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Trigger Efficiency', fontsize=14)

            # Set labels
            ax.set_xlabel(var_info['label_x'], fontsize=14)
            ax.set_ylabel(var_info['label_y'], fontsize=14)

            # Add title with overall efficiency
            trigger_label = trigger if trigger != 'OR' else 'OR of all triggers'
            ax.set_title(f"{trigger_label}\nOverall Efficiency: {total_efficiency:.2%}", fontsize=14)

            # Add CMS label
            hep.cms.label(ax=ax, data=data, year=year, com="13.6")

            # Save figure
            trigger_safe_name = trigger.replace('/', '_').replace(' ', '_')
            save_path = os.path.join(save_dir, f"{trigger_safe_name}_{var_name}_2d.png")
            plt.savefig(save_path, dpi=200, bbox_inches='tight')

            if show_img:
                plt.show()
            else:
                print(f"Saved {save_path}")

            plt.close(fig)


def plot_cutflow(cutflow, tag, save_dir="./figures_cutflow", show_img=True):
    """
    Plot a horizontal-bar cutflow diagram.

    Parameters:
    - cutflow: dict
        Ordered dict mapping step labels to event counts.
    - tag: str
        Production mode label (e.g. 'VBF', 'ggF'), used in the title and filename.
    - save_dir: str
        Directory to save the figure.
    - show_img: bool
        Whether to display the figure interactively.
    """
    steps  = list(cutflow.keys())
    counts = [cutflow[s] for s in steps]
    total  = counts[0] if counts else 1

    fig, ax = plt.subplots(figsize=(10, 0.5 * len(steps) + 2))

    y_pos = np.arange(len(steps))
    bars  = ax.barh(y_pos, counts, color='#4C72B0', edgecolor='black', height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(steps, fontsize=11)
    ax.invert_yaxis()
    ax.set_xlabel("Number of events", fontsize=13)
    ax.set_title(f"Cutflow  —  {tag}", fontsize=14)

    # Annotate each bar with count and fraction relative to "All events"
    x_max = counts[0] if counts else 1
    for bar, count in zip(bars, counts):
        frac = count / total if total > 0 else 0
        ax.text(
            bar.get_width() + x_max * 0.012,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}  ({frac:.1%})",
            va='center', ha='left', fontsize=10
        )

    ax.set_xlim(0, x_max * 1.35)
    plt.tight_layout()

    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, f"cutflow_{tag}.png")
    plt.savefig(save_path, dpi=200, bbox_inches='tight')

    if show_img:
        plt.show()
    else:
        print(f"Saved {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="2D Trigger Efficiency Analysis")
    parser.add_argument('--year', required=False, choices=list(trigger_dict_periods.keys()), default='2022',
                        help="Year of the dataset (e.g., 2022, 2023, etc.).")
    parser.add_argument('--test', action='store_true',
                        help="Test mode: only process first 10 ROOT files.")
    args = parser.parse_args()

    year = args.year
    test_mode = args.test

    if test_mode:
        print(f'Running in TEST mode: only processing first 10 files')
    print(f'Processing period {year}')

    all_triggers = trigger_dict_periods[year]

    # Process both VBF and ggF production modes
    prod_modes = ['VBF', 'ggF']

    for prod_mode in prod_modes:
        print(f"\n{'='*60}")
        print(f"Processing production mode: {prod_mode}")
        print(f"{'='*60}\n")

        # Load fileset
        fileset = {}
        json_path = os.path.join(f'infiles/{year}', f"{year}_{prod_mode}.json")
        print(f"Reading: {json_path}")

        with open(json_path, 'r') as file:
            data = json.load(file)

        # Get the first key in the JSON (assuming single dataset per file)
        file_name = list(data.keys())[0] if data else None
        if file_name:
            samples = data.get(file_name, [])

            # Apply test mode limit if enabled
            if test_mode:
                samples = samples[:10]
                print(f"TEST MODE: Limited to {len(samples)} files")

            dataset_key = f"{prod_mode}_{year}"
            fileset[dataset_key] = samples
            print(f"Loaded {len(samples)} files for dataset '{dataset_key}'")
        else:
            print(f"Error: No data found in {json_path}")
            continue  # Skip to next production mode

        print("Constructed fileset:")
        print(json.dumps(fileset, indent=2))

        # Run processor
        iterative_run = processor.Runner(
            executor=processor.FuturesExecutor(compression=None, workers=2),
            schema=NanoAODSchema,
            skipbadfiles=True,
            savemetrics=True,
        )

        print(f"\nRunning 2D trigger efficiency processor for {prod_mode}...")
        out = iterative_run(
            fileset,
            treename="Events",
            processor_instance=TriggerEfficiency2DProcessor(all_triggers, trig_vars_2d, baseline_key=prod_mode),
        )
        output = out[0]

        # Save output
        FIGURES_DIR = os.path.join(os.getcwd(), "figures_2d/", f"{year}/{prod_mode}/")
        os.makedirs(FIGURES_DIR, exist_ok=True)

        output_path = os.path.join(os.getcwd(), f"output/trig_eff_2d_{year}_{prod_mode}.coffea")
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        util.save(output, output_path)
        print(f"Saved output to {output_path}")

        # Generate plots
        plot_2d_trigger_efficiency(
            output,
            trig_vars_2d,
            save_dir=FIGURES_DIR,
            year=year,
            data=False,
            show_img=False
        )

        # Generate cutflow plot and print table
        CUTFLOW_DIR = os.path.join(os.getcwd(), "figures_cutflow/", f"{year}/{prod_mode}/")
        plot_cutflow(output['cutflow'], tag=prod_mode, save_dir=CUTFLOW_DIR, show_img=False)

        print(f"\n{prod_mode} cutflow:")
        total_events = output['cutflow']['All events']
        for step, count in output['cutflow'].items():
            frac = count / total_events if total_events > 0 else 0
            print(f"  {step:<55} {count:>10,} ({frac:.1%})")

        print(f"\nAll plots for {prod_mode} saved to {FIGURES_DIR}")

    print(f"\n{'='*60}")
    print(f"Completed processing for all production modes")
    print(f"{'='*60}")
