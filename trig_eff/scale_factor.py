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
            '2022': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
        }
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
            "cutflow": processor.defaultdict_accumulator(int),
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
        # initialize cutflow accumulator for this chunk (so it merges back)
        output['cutflow'] = processor.defaultdict_accumulator(int)
        dataset = events.metadata['dataset']
        is_mc = hasattr(events, "genWeight")
        output['counts'] = len(events)
        reference_trigger_paths = self._reference_triggers[self.year]
        fired_reference_triggers = ak.zeros_like(events.event, dtype=bool)
        for trig in reference_trigger_paths:
            tfield = trig[4:] if trig.startswith('HLT_') else trig
            if tfield in events.HLT.fields:
                fired_reference_triggers = fired_reference_triggers | events.HLT[tfield]
        events = events[fired_reference_triggers]
        # Cutflow: after reference triggers
        output['cutflow'][f"{dataset}::ref_trig"] += int(len(events))
        n_or = ak.count_nonzero(fired_reference_triggers)
        # Denominator selection:
        muons = events.Muon
        # Loose muon: low pT, loose ID, modest isolation
        loose_muon_selection = (
            (muons.pt > 25) & (abs(muons.eta) < 2.4) & (muons.looseId) & (muons.pfRelIso04_all < 0.25)
        )
        loose_muons = muons[loose_muon_selection]
        has_loose_muon = ak.num(loose_muons) >= 1
        # Cutflow: after loose muon
        output['cutflow'][f"{dataset}::loose_muon"] += int(ak.count_nonzero(has_loose_muon))
        # Veto jets near the muon: compute ΔR(mu, leading FatJet) and require > 0.8
        leading_fj = ak.firsts(ak.pad_none(events.FatJet, 1, clip=True))
        # If no FatJet present in an event, treat as passing the ΔR veto
        dphi = np.abs(((loose_muons.phi - leading_fj.phi + np.pi) % (2 * np.pi)) - np.pi)
        deta = loose_muons.eta - leading_fj.eta
        dr = np.sqrt(dphi**2 + deta**2)
        # Pass if ANY loose muon is farther than 0.8 from the leading FatJet
        far_enough = ak.fill_none(ak.any(dr > 0.8, axis=1), True)
        denom_mask = has_loose_muon & far_enough
        # Cutflow: after far_enough
        output['cutflow'][f"{dataset}::far_enough"] += int(ak.count_nonzero(denom_mask))
        # Apply loose muon selection
        # keep events filtered only by the OR of reference muon triggers
        denominator_events = events[denom_mask]
        # Build variables once
        var_values_den = {name: cfg['proc'](denominator_events) for name, cfg in trig_vars_data.items()}
        # Always fill denominator first for each variable
        weights_den = denominator_events.genWeight if is_mc else np.ones(len(denominator_events))
        # Try with uniform weights
        # weights_den = np.ones(len(denominator_events))
        for name, values in var_values_den.items():
            output[f'eff_{name}'].fill(dataset=dataset, selection="denominator", **{name: values}, weight=weights_den)
        # Track numerator count for logging
        numerator_count = 0
        # Conditionally fill numerator if the trigger exists and fired
        t_interest = self.trigger_of_interest[4:] if self.trigger_of_interest.startswith('HLT_') else self.trigger_of_interest
        if t_interest in denominator_events.HLT.fields:
            fired_trigger_of_interest = denominator_events.HLT[t_interest]
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
        key_or = f"{dataset}::{fname}::or"  
        output['selected_counts'][key_den] += int(len(denominator_events))
        output['selected_counts'][key_num] += int(numerator_count)
        output['selected_counts'][key_or] += int(n_or)
        return output

    def postprocess(self, accumulator):
        return accumulator

def plot_trigger_sfs(accumulator, year, trigger_name, data_dataset, mc_dataset, save_dir):
    """
    Plots trigger efficiencies and scale factors from the processor output.
    - Shows Data and MC efficiencies.
    - Ratio (Data/MC) masks undefined bins (denominators==0) and adapts y-limits to visible data.
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
                # Cast to float for safe division downstream
                return h[ds, sel, :].values().astype(float)
            except Exception:
                return np.zeros(len(axis.edges) - 1, dtype=float)

        # numerator/denominator counts per bin
        data_num = safe_project(data_dataset, "numerator")
        data_den = safe_project(data_dataset, "denominator")
        mc_num   = safe_project(mc_dataset,   "numerator")
        mc_den   = safe_project(mc_dataset,   "denominator")

        # Efficiencies: undefined bins -> NaN (so the line breaks there)
        eff_data = np.divide(
            data_num, data_den,
            out=np.full_like(data_num, np.nan, dtype=float),
            where=(data_den > 0)
        )
        eff_mc = np.divide(
            mc_num, mc_den,
            out=np.full_like(mc_num, np.nan, dtype=float),
            where=(mc_den > 0)
        )

        # Scale factor (ratio): only when both dens>0 and efficiencies are finite
        valid_ratio = (data_den > 0) & (mc_den > 0) & np.isfinite(eff_data) & np.isfinite(eff_mc)
        sf = np.divide(
            eff_data, eff_mc,
            out=np.full_like(eff_data, np.nan, dtype=float),
            where=valid_ratio
        )

        # --- plotting ---
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(10, 10),
            gridspec_kw={"height_ratios": (3, 1)},
            sharex=True
        )
        fig.subplots_adjust(hspace=0.05)
        fontsize = 14

        # Top: efficiencies
        hep.histplot(eff_data, bins=axis.edges, yerr=False, ax=ax, histtype='step', color='k', label='Data')
        hep.histplot(eff_mc,   bins=axis.edges, yerr=False, ax=ax, histtype='step', color='r', linewidth=2, label='MC (Sim.)')
        ax.set_ylabel("Efficiency", fontsize=fontsize)
        ax.set_ylim(0.0, 1.1)
        ax.legend(fontsize=9)
        # Title above CMS label as figure-level title
        fig.suptitle(f"{trigger_name}", x=0.98, ha='right', fontsize=fontsize, fontstyle='italic')
        hep.cms.label("Preliminary", data=True, year=year, ax=ax, lumi=27.7)
        fig.subplots_adjust(top=0.88)
        ax.tick_params(labelsize=9)

        # Bottom: Data/MC ratio with adaptive y-limits
        hep.histplot(sf, bins=axis.edges, ax=rax, histtype='step', color='k', yerr=False)
        rax.set_ylabel("Data / MC", fontsize=fontsize)
        rax.set_xlabel(axis.label, fontsize=fontsize)
        rax.grid(True)
        rax.tick_params(labelsize=9)

        finite_sf = np.isfinite(sf)
        if np.any(finite_sf):
            lo = float(np.nanmin(sf[finite_sf]))
            hi = float(np.nanmax(sf[finite_sf]))
            if hi == lo:
                pad = 0.15 * max(1.0, abs(hi))
                rax.set_ylim(lo - pad, hi + pad)
            else:
                pad = 0.15 * (hi - lo)
                # keep ratio sensible (non-negative) if possible
                rax.set_ylim(max(0.0, lo - pad), hi + pad)
        else:
            # fallback when no valid ratio bins
            rax.set_ylim(0.8, 1.2)

        # Save
        save_path = os.path.join(save_dir, f"{trigger_name}/{name}.png")
        print(f"  - Saving 1D plot to: {save_path}")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Also render a simple cutflow bar plot per dataset if present
    if 'cutflow' in accumulator:
        cf = accumulator['cutflow']
        # group by dataset
        by_dataset = {}
        for key, val in cf.items():
            ds, step = key.split('::', 1)
            by_dataset.setdefault(ds, {})[step] = val
        for ds, steps in by_dataset.items():
            labels = ['ref_trig', 'loose_muon', 'far_enough']
            values = [int(steps.get(k, 0)) for k in labels]
            fig, ax = plt.subplots(figsize=(7, 4))
            x = np.arange(len(labels))
            ax.bar(x, values, color='tab:purple')
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=20, ha='right')
            ax.set_ylabel('Events', fontsize=10)
            ax.set_title(f"Cutflow: {ds} ({year})", fontsize=10)
            ax.tick_params(labelsize=9)
            for xi, v in zip(x, values):
                ax.annotate(str(v), (xi, v), textcoords='offset points', xytext=(0, 3), ha='center', fontsize=7)
            fig.subplots_adjust(left=0.12, right=0.98, bottom=0.22, top=0.88)
            hep.cms.label("Work in Progress", data=('Muon' in ds or 'SingleMuon' in ds), year=year, ax=ax, fontsize=9)
            out_path = os.path.join(save_dir, f"cutflow_{ds}.png")
            print(f"  - Saving cutflow to: {out_path}")
            plt.savefig(out_path, dpi=150, bbox_inches='tight')
            plt.close(fig)


muon_triggers = ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"]
class TriggerCutflowProcessor(processor.ProcessorABC):
    """Count events passing each trigger in trigger_dict_periods[year] and their OR."""
    def __init__(self, year: str):
        self.year = str(year)
        self.triggers = muon_triggers
        self._accumulator = processor.dict_accumulator({
            "cutflow": processor.defaultdict_accumulator(int),  # keys: f"{dataset}::{trigger}"
            "total": processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata.get('dataset', 'unknown')
        output['total'][dataset] += int(len(events))

        # Count per-trigger with awkward's count_nonzero (robust for booleans)
        # Note: NanoAODSchema drops the 'HLT_' prefix in events.HLT fields
        or_mask = None
        for tname in self.triggers:
            tfield = tname[4:] if tname.startswith('HLT_') else tname
            if tfield not in events.HLT.fields:
                continue
            mask = events.HLT[tfield]
            n_fired = int(ak.count_nonzero(mask))
            output['cutflow'][f"{dataset}::{tname}"] += n_fired
            or_mask = mask if or_mask is None else (or_mask | mask)

        # OR of all present triggers
        if or_mask is not None:
            output['cutflow'][f"{dataset}::OR"] += int(ak.count_nonzero(or_mask))
        else:
            output['cutflow'][f"{dataset}::OR"] += 0

        return output

    def postprocess(self, accumulator):
        return accumulator


def plot_trigger_cutflow(accumulator, year: str, save_dir: str):
    """Plot bar charts of counts per trigger and their OR for each dataset."""
    os.makedirs(save_dir, exist_ok=True)
    cutflow = accumulator['cutflow']

    dataset_to_counts = {}
    for key, val in cutflow.items():
        dataset, trigger = key.split('::', 1)
        dataset_to_counts.setdefault(dataset, {})[trigger] = val

    for dataset, counts in dataset_to_counts.items():
        # Use the configured muon triggers for label order
        trig_list = muon_triggers
        labels = trig_list + ['OR']
        values = [counts.get(t, 0) for t in trig_list] + [counts.get('OR', 0)]

        # Bar graph with generous margins so labels are not cut off
        fig, ax = plt.subplots(figsize=(max(9, 1.2 * len(labels)), 6))
        x = np.arange(len(labels))
        ax.bar(x, values, color='tab:blue')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha='right')
        ax.set_ylabel('Events passing')
        ax.set_title(f"Cutflow counts for {dataset} ({year})", fontsize=12)
        # small margins around plot
        ax.margins(x=0.02, y=0.10)
        # annotate counts above bars
        for xi, v in zip(x, values):
            ax.annotate(str(v), xy=(xi, v), xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)
        # adjust layout so nothing is clipped
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.26, top=0.88)
        hep.cms.label("Work in Progress", data=('Muon' in dataset or 'SingleMuon' in dataset), year=year, ax=ax)
        out_path = os.path.join(save_dir, f"cutflow_{dataset}.png")
        print(f"  - Saving cutflow to: {out_path}")
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)

file_dict_periods = {'2022': {'ttbar': "TTto4Q", 'MuonData': "Muon_Run2022C"},
                     '2022EE': {'ttbar': "TTto4Q", 'MuonData': "Muon_Run2022E"},
                     '2023': {'ttbar': "TTto4Q", 'MuonData': "Muon0_Run2023C-v4"},
                     '2023BPix': {'ttbar': "TTto4Q", 'MuonData': "Muon0_Run2023D-v1"}
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
    # e.g. file_dict_periods['2022'] = {'ttbar': 'TTto4Q', 'MuonData': 'Muon_Run2022C'}
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

    # --- Quick inspection: load first data file and print available HLT paths ---
    try:
        data_ds_key = next(k for k in fileset.keys() if ('Muon' in k or 'SingleMuon' in k))
        first_file = fileset[data_ds_key][0]
        print(f"\nInspecting HLT branches from first data file: {first_file}")
        with uproot.open(first_file) as fin:
            events_tree = fin["Events"]
            # uproot 4 filter_name uses shell-style wildcards
            try:
                hlt_branches = list(events_tree.keys(filter_name="HLT_*"))
            except Exception:
                branch_names = list(events_tree.keys())
                hlt_branches = [b for b in branch_names if str(b).startswith("HLT_")]
            hlt_branches = sorted(str(b) for b in hlt_branches)
            print(f"Found {len(hlt_branches)} HLT paths:")
            for name in hlt_branches:
                print(f"  {name}")
    except StopIteration:
        print("No Muon data dataset found in fileset; skipping HLT inspection.")
    except Exception as exc:
        print(f"Failed to inspect HLT branches: {exc}")

    
    # --- Quick inspection: load first ttbar MC file and print available HLT paths ---
    try:
        mc_ds_key = next(k for k in fileset.keys() if ('ttbar' in k or 'TT' in k))
        first_mc_file = fileset[mc_ds_key][0]
        print(f"\nInspecting HLT branches from first MC file: {first_mc_file}")
        with uproot.open(first_mc_file) as fin:
            events_tree = fin["Events"]
            try:
                hlt_branches = list(events_tree.keys(filter_name="HLT_*"))
            except Exception:
                branch_names = list(events_tree.keys())
                hlt_branches = [b for b in branch_names if str(b).startswith("HLT_")]
            hlt_branches = sorted(str(b) for b in hlt_branches)
            print(f"Found {len(hlt_branches)} HLT paths (MC):")
            for name in hlt_branches:
                print(f"  {name}")
    except StopIteration:
        print("No ttbar MC dataset found in fileset; skipping HLT inspection.")
    except Exception as exc:
        print(f"Failed to inspect HLT branches (MC): {exc}")

    # runner = processor.Runner(executor=processor.FuturesExecutor(workers=2), schema=NanoAODSchema)
    # out = runner(fileset, "Events", TriggerCutflowProcessor(year=YEAR))
    # plot_trigger_cutflow(out, YEAR, os.path.join(FIGURES_DIR, "cutflow"))

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
        mc_dataset=f'ttbar_{YEAR}',
        save_dir=FIGURES_DIR
    )
    print("\nAll steps completed.")