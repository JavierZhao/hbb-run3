import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import awkward as ak
import uproot
import hist
import numpy as np
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema

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
import json

class TriggerEfficiencyProcessor(processor.ProcessorABC):
    def __init__(self, triggers):
        self.triggers = triggers
        pt_axis = hist.axis.Regular(bins=100, start=0, stop=3000, name="pt", label="Leading Jet $p_{T}$ [GeV]")
        msd_axis = hist.axis.Regular(bins=100, start=0, stop=400, name="msd", label="Leading Jet $m_{SD}$ [GeV]")

        self.output = processor.dict_accumulator({
            trigger: {
                'pt_total': hist.Hist(pt_axis),
                'pt_pass': hist.Hist(pt_axis),
                'msd_total': hist.Hist(msd_axis),
                'msd_pass': hist.Hist(msd_axis)
            } for trigger in triggers
        })

    def process(self, events):
        output = self.output.copy()
        
        # Ensure there are events with at least one fat jet
        has_fatjets = ak.num(events.FatJet) > 0
        leading_fatjet_pt = events.FatJet[has_fatjets].pt[:, 0]  # Leading jet pT
        leading_fatjet_msd = events.FatJet[has_fatjets].msoftdrop[:, 0] # Leading jet mSD 
        baseline = ak.any(events.FatJet[has_fatjets].pt > 30, axis=1)  

        # Loop through each trigger
        for trigger in self.triggers:
            trigger_mask = events.HLT[trigger]
            output[trigger]['pt_total'].fill(pt=leading_fatjet_pt[baseline])
            output[trigger]['pt_pass'].fill(pt=leading_fatjet_pt[baseline & trigger_mask[has_fatjets]])
            output[trigger]['msd_total'].fill(msd=leading_fatjet_msd[baseline])  
            output[trigger]['msd_pass'].fill(msd=leading_fatjet_msd[baseline & trigger_mask[has_fatjets]])

        return output

    def postprocess(self, accumulator):
        return accumulator

def plot_trigger_efficiencies(trig_name, prod_mode, output, triggers_set, max_triggers_per_plot=5):
    save_dir = f"figures/{prod_mode}/leading_pt_mSD"
    os.makedirs(save_dir, exist_ok=True)
    # Function to chunk the triggers_set
    def chunk_triggers(triggers, chunk_size):
        for i in range(0, len(triggers), chunk_size):
            yield triggers[i:i + chunk_size]

    colors = ['blue', 'green', 'red', 'purple', 'orange']  # Adjust or add more colors as needed

    # Chunk the triggers_set into parts of size max_triggers_per_plot
    trigger_chunks = list(chunk_triggers(triggers_set, max_triggers_per_plot))

    for chunk_index, trigger_chunk in enumerate(trigger_chunks):
        fig, axs = plt.subplots(1, 2, figsize=(20, 8))
        
        # First subplot for H_T
        for i, trigger in enumerate(trigger_chunk):
            ht_total = output[trigger]['pt_total']
            ht_pass = output[trigger]['pt_pass']
            efficiency = ht_pass / ht_total
            efficiency.plot(ax=axs[0], label=f'{trigger}', color=colors[i % len(colors)])  # Plot directly on the first subplot

        axs[0].set_xlabel(r'Leading jet $p_T$')
        axs[0].set_ylabel('Efficiency')

        # Second subplot for m_H
        for i, trigger in enumerate(trigger_chunk):
            mh_total = output[trigger]['msd_total']
            mh_pass = output[trigger]['msd_pass']
            efficiency = mh_pass / mh_total
            efficiency.plot(ax=axs[1], label=f'{trigger}', color=colors[i % len(colors)])  # Plot directly on the second subplot

        axs[1].set_xlabel(r'Leading jet $m_{SD}$')
        axs[1].set_ylabel('Efficiency')

        # Place a shared legend outside the plots on the right side
        handles, labels = axs[0].get_legend_handles_labels()  # Only need to extract handles and labels from one axis
        fig.legend(handles, labels, title='Triggers', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # Adjust layout to accommodate legends outside of plots
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust rect if needed to fit legend
        plt.title(prod_mode)
        plt.savefig(f"{save_dir}/{trig_name}_TriggerEfficiencies_Set{chunk_index + 1}.png", dpi=200, bbox_inches='tight')
        plt.show()


def main(args):
    file_dict = {'VBF': 'VBF_Hto2B', 'ggF': 'GluGlu_Hto2B'}
    with open(f'2022/2022_{args.file_name}.json', 'r') as file:
        data = json.load(file)
        samples = data[file_dict[args.file_name]]
        sample = samples[1]

    # -------------------------------
    #        Triggers
    # -------------------------------  
    triggers_cmantill = ['AK8DiPFJet250_250_MassSD50',
     'AK8DiPFJet260_260_MassSD30',
     'AK8PFJet230_SoftDropMass40',
     'AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet400_SoftDropMass40',
     'AK8PFJet425_SoftDropMass40',
     'PFHT1050',
     'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65']
    triggers_pnet = ['QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',
     'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65',
     'QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65',
     'AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35',
     'AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30',
     'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30',
     'AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30']
    triggers_VBF_CSV = [trig[4:-2] for trig in [
    "HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",
    "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",
    "HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",
    "HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v"]]
    triggers_VBF_Jet =[trig[4:-2] for trig in[       "HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v",
        "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v",
        "HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v",
        "HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v",              
    ]]
    triggers_dict = {'Cristina' : triggers_cmantill,
                    'ParticleNet' : triggers_pnet,
                    'VBF_CSV': triggers_VBF_CSV,
                    'VBF_Jet': triggers_VBF_Jet}
    
    for name, triggers_set in triggers_dict.items():
        print(name)
        iterative_run = processor.Runner(
        executor = processor.FuturesExecutor(compression=None, workers=2),
        schema=NanoAODSchema,
        skipbadfiles=True,  # Skip files that fail to open
        savemetrics=True,   # Save metrics to understand where failures occur
    )
        fileset = {"Dataset" : samples}
        out = iterative_run(
            fileset,
            treename="Events",
            processor_instance=TriggerEfficiencyProcessor(triggers_set),
        )
        output = out[0]
        plot_trigger_efficiencies(name, args.file_name, output, triggers_set)
        


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file-name",
        type=str,
        action="store",
        default="VBF",
        help="which json file to load (which production mode)",
    )
    
    args = parser.parse_args()
    main(args)