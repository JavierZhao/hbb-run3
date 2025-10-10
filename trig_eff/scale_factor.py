import awkward as ak
import numpy as np
import coffea.processor as processor
from coffea.nanoevents import NanoAODSchema
import hist
import uproot
import json
import os
import shutil

# Plotting libraries (ensure you have them installed: pip install matplotlib mplhep)
import matplotlib.pyplot as plt
import mplhep as hep
# Apply the CMS plotting style
plt.style.use(hep.style.CMS)


trig_vars_data = {
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
        'axis': hist.axis.Regular(bins=10, start=40, stop=200, name="msd", label="Leading Jet $m_{SD}$ [GeV]"),
        'proc': lambda events: ak.fill_none(ak.pad_none(events.FatJet.msoftdrop, 1, clip=True)[:, 0], np.nan)
    },
    'num_ak4': {
        'label': "Number of AK4 Jets",
        'axis': hist.axis.Integer(0, 20, name="num_ak4", label="Number of AK4 Jets"),
        'proc': lambda events: ak.num(events.Jet)
    },
    'particleNet_XbbVsQCD':{
        'label': "Leading Particle Net TXbb score",
        'axis': hist.axis.Regular(bins=30, start=0, stop=1, name="particleNet_XbbVsQCD", label="Leading Particle Net TXbb score"),
        'proc': lambda events: ak.fill_none(
            ak.pad_none(events.FatJet.particleNet_XbbVsQCD, 1, clip=True)[:, 0],
            np.nan
        )
    }
}

class MuonTriggerSFProcessor(processor.ProcessorABC):
    """
    A coffea processor to calculate muon trigger scale factors.
    """
    def __init__(self, year='2022', trigger_of_interest='HLT_Mu50'):
        self.year = str(year)
        self.trigger_of_interest = trigger_of_interest
        self._reference_triggers = {
            # '2022': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100"],
            '2022':["HLT_IsoMu24"],
            # 	HLT_Mu50_v* OR HLT_CascadeMu100_v* OR HLT_HighPtTkMu100_v*
        }
        # HLT_Mu50_v* OR HLT_CascadeMu100_v* OR HLT_HighPtTkMu100_v*
        if self.year not in self._reference_triggers:
            raise ValueError(f"Year '{self.year}' not supported. Please add reference triggers.")
        self._accumulator = processor.dict_accumulator({
            "counts": processor.value_accumulator(int, initial=0),
            # Per-variable histograms: dataset, selection, variable axis
            **{
                f"eff_{name}": hist.Hist(
                    hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True),
                    hist.axis.StrCategory(["numerator", "denominator"], name="selection", label="Selection"),
                    cfg["axis"],
                ) for name, cfg in trig_vars_data.items()
            },
            "selected_counts": processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        # Build a fresh output; avoid calling identity() on scikit-hep hist
        output = {
            'counts': 0,
        }
        for _name in trig_vars_data.keys():
            output[f'eff_{_name}'] = self._accumulator[f'eff_{_name}'].copy()
        # initialize per-file selection counts accumulator for this chunk
        output['selected_counts'] = processor.defaultdict_accumulator(int)
        dataset = events.metadata['dataset']
        is_mc = hasattr(events, "genWeight")
        output['counts'] = len(events)
        reference_trigger_paths = self._reference_triggers[self.year]
        fired_reference_triggers = ak.zeros_like(events.event, dtype=bool)
        for trig in reference_trigger_paths:
            if trig in events.HLT.fields:
                fired_reference_triggers = fired_reference_triggers | events.HLT[trig]
        # if not (is_mc and ("QCD" in dataset)):
        #     events = events[fired_reference_triggers]
        # Build denominator: for QCD MC, do NOT require a tight muon; for others, keep muon selection
        # is_qcd_mc = is_mc and ("QCD" in dataset)
        # if is_qcd_mc:
        #     denominator_events = events
        # else:
        muons = events.Muon
        tight_muon_selection = (
            (muons.pt > 25) & (abs(muons.eta) < 2.4) & (muons.tightId) & (muons.pfRelIso04_all < 0.15)
        )
        tight_muons = muons[tight_muon_selection]
        has_one_tight_muon = ak.num(tight_muons) >= 1
        denominator_events = events[has_one_tight_muon]
        # Build variables once
        var_values_den = {name: cfg['proc'](denominator_events) for name, cfg in trig_vars_data.items()}
        # Always fill denominator first for each variable
        weights_den = denominator_events.genWeight if is_mc else np.ones(len(denominator_events))
        for name, values in var_values_den.items():
            output[f'eff_{name}'].fill(dataset=dataset, selection="denominator", **{name: values}, weight=weights_den)
        # Track numerator count for logging
        numerator_count = 0
        # Conditionally fill numerator if the trigger exists and fired
        if self.trigger_of_interest in denominator_events.HLT.fields:
            fired_trigger_of_interest = denominator_events.HLT[self.trigger_of_interest]
            numerator_events = denominator_events[fired_trigger_of_interest]
            if len(numerator_events) > 0:
                var_values_num = {name: cfg['proc'](numerator_events) for name, cfg in trig_vars_data.items()}
                weights_num = numerator_events.genWeight if is_mc else np.ones(len(numerator_events))
                numerator_count = len(numerator_events)
                for name, values in var_values_num.items():
                    output[f'eff_{name}'].fill(dataset=dataset, selection="numerator", **{name: values}, weight=weights_num)
        else:
            print(f"Warning: Trigger '{self.trigger_of_interest}' not found in HLT paths for dataset {dataset}! Skipping numerator fill.")
        # Record per-file counts into the accumulator so we can print from the main process
        try:
            filename = events.metadata.get('filename', '')
        except Exception:
            filename = ''
        fname = os.path.basename(filename) if filename else 'unknown'
        key_den = f"{dataset}::{fname}::den"
        key_num = f"{dataset}::{fname}::num"
        output['selected_counts'][key_den] += int(len(denominator_events))
        output['selected_counts'][key_num] += int(numerator_count)
        return output

    def postprocess(self, accumulator):
        return accumulator

def plot_trigger_sfs(accumulator, year, trigger_name, data_dataset, mc_dataset, save_dir):
    """
    Plots trigger efficiencies and scale factors from the processor output.
    """
    print(f"\n--- Generating plots for trigger '{trigger_name}' ---")
    os.makedirs(save_dir, exist_ok=True)
    # Print per-file selection counts collected from workers
    if 'selected_counts' in accumulator:
        print("Per-file selected counts (den/num):")
        for key, val in accumulator['selected_counts'].items():
            print(f"  {key} -> {val}")

    for name, cfg in trig_vars_data.items():
        h = accumulator[f'eff_{name}']
        axis = h.axes[name]

        def safe_project(ds, sel):
            try:
                return h[ds, sel, :].values()
            except Exception:
                return np.zeros(len(axis.edges) - 1)

        data_num = safe_project(data_dataset, "numerator")
        data_den = safe_project(data_dataset, "denominator")
        mc_num = safe_project(mc_dataset, "numerator")
        mc_den = safe_project(mc_dataset, "denominator")

        eff_data = np.divide(data_num, data_den, out=np.zeros_like(data_num), where=data_den!=0)
        eff_mc = np.divide(mc_num, mc_den, out=np.zeros_like(mc_num), where=mc_den!=0)
        sf = np.divide(eff_data, eff_mc, out=np.ones_like(eff_data), where=eff_mc!=0)

        fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=0.05)
        yerr_flag = bool(np.any(data_den > 0))
        hep.histplot(eff_data, bins=axis.edges, yerr=yerr_flag, ax=ax, histtype='errorbar', color='k', label='Data')
        hep.histplot(eff_mc, bins=axis.edges, ax=ax, histtype='step', color='r', linewidth=2, label='MC (Sim.)')
        ax.set_ylabel("Efficiency")
        ax.set_ylim(0, 1.1)
        ax.legend()
        hep.cms.label(f"Preliminary", data=True, year=year, ax=ax, lumi=27.7)
        ax.set_title(f"{trigger_name}", loc='right', style='italic')
        hep.histplot(sf, bins=axis.edges, ax=rax, histtype='errorbar', color='k', yerr=False)
        rax.set_ylabel("Data / MC")
        rax.set_xlabel(axis.label)
        rax.set_ylim(0.8, 1.2)
        rax.grid(True)
        save_path = os.path.join(save_dir, f"{trigger_name}/{name}.png")
        print(f"  - Saving 1D plot to: {save_path}")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

file_dict_periods = {'2022': {'QCD': "QCD_HT100to200", 'MuonData': "Muon_Run2022C"},
                     '2022EE': {'QCD': "QCD_HT100to200", 'MuonData': "Muon_Run2022E"},
                     '2023': {'QCD': "QCD_HT100to200", 'MuonData': "Muon0_Run2023C-v4"},
                     '2023BPix': {'QCD': "QCD_HT100to200", 'MuonData': "Muon0_Run2023D-v1"}
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


if __name__ == '__main__':

    # --- 1. Setup for Analysis ---
    print("--- Step 1: Setting up analysis environment ---")
    YEAR = '2022'
    TRIGGER_OF_INTEREST = 'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35'
    FIGURES_DIR = os.path.join(os.getcwd(), "figures/SF/", YEAR)

    # --- 2. File Loading from local JSONs and Processor Execution ---
    print("--- Step 2: Running processor ---")

    # Expect JSONs at /srv/{YEAR}/{YEAR}_{prod_mode}.json with keys from file_dict_periods[YEAR]
    # e.g. file_dict_periods['2022'] = {'QCD': 'QCD_HT100to200', 'MuonData': 'Muon_Run2022C'}
    prod_modes_map = file_dict_periods[YEAR]

    fileset = {}
    for prod_mode, file_key in prod_modes_map.items():
        json_path = os.path.join(YEAR, f"{YEAR}_{prod_mode}.json")
        print(f"Reading: {json_path} (key: {file_key})")
        with open(json_path, 'r') as file:
            data = json.load(file)
        samples = data.get(file_key, [])
        dataset_key = f"{prod_mode}_{YEAR}"
        fileset[dataset_key] = samples[:5]

    print("Constructed fileset:")
    print(json.dumps(fileset, indent=2))

    iterative_run = processor.Runner(
            executor = processor.FuturesExecutor(compression=None, workers=2),
            schema=NanoAODSchema,
            skipbadfiles=True,
    )
    
    print("\nExecuting processor... (This may take a few moments with remote files)")
    output = iterative_run(
        fileset,
        treename="Events",
        processor_instance=MuonTriggerSFProcessor(year=YEAR, trigger_of_interest=TRIGGER_OF_INTEREST),
    )
    print("Processing complete.")

    # --- 3. Plotting Results ---
    plot_trigger_sfs(
        accumulator=output,
        year=YEAR,
        trigger_name=TRIGGER_OF_INTEREST,
        data_dataset=f'MuonData_{YEAR}',
        mc_dataset=f'QCD_{YEAR}',
        save_dir=FIGURES_DIR
    )
    print("\nAll steps completed.")