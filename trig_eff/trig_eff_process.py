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

def create_baseline_selection_mask(events, tag):
    n_events = len(events)
    
    # Initialize cumulative event mask (all events start as True)
    baseline_mask = np.ones(n_events, dtype=bool)
    
    # Step 1: Select jets with pT > 15 GeV
    jets = events.Jet[events.Jet.pt > 15]
    has_jets = ak.num(jets) >= 1  # Events with at least one jet
    has_jets = ak.to_numpy(has_jets)  # Convert to NumPy array
    
    # Update baseline mask
    baseline_mask &= has_jets
    
    # Apply mask to events and jets
    # events = events[baseline_mask]
    # jets = jets[baseline_mask]
    
    # Sort jets by descending pT
    jets = jets[ak.argsort(jets.pt, ascending=False)]

    n_events = len(events)

    if tag == "VBF":
        # Step 2: Ensure there are at least four jets
        has_four_jets = ak.num(jets) >= 4  # Events with at least four jets
        has_four_jets = ak.to_numpy(has_four_jets)  # Convert to NumPy array
        
        # Update baseline mask
        baseline_mask[baseline_mask] &= has_four_jets  # Update only where baseline_mask is True
        
        # Apply mask to events and jets
        events = events[has_four_jets]
        jets = jets[has_four_jets]
        
        # Step 3: Apply pT thresholds to the leading jets
        leading_jet_pt_mask = (
            (jets[:, 0].pt > 105) &
            (jets[:, 1].pt > 88) &
            (jets[:, 2].pt > 76) &
            (jets[:, 3].pt > 15)
        )
        leading_jet_pt_mask = ak.to_numpy(leading_jet_pt_mask)  # Convert to NumPy array
        
        # Update baseline mask
        baseline_mask_indices = np.where(baseline_mask)[0]
        baseline_mask[baseline_mask_indices] &= leading_jet_pt_mask
        
        # Apply mask to events and jets
        events = events[leading_jet_pt_mask]
        jets = jets[leading_jet_pt_mask]
        
        # Step 4: Apply eta ranges
        eta_mask = (
            (abs(jets[:, 0].eta) < 5.0) &
            (abs(jets[:, 1].eta) < 5.0) &
            (abs(jets[:, 2].eta) < 5.0) &
            (abs(jets[:, 3].eta) < 5.0)
        )
        eta_mask = ak.to_numpy(eta_mask)  # Convert to NumPy array
        
        # Update baseline mask
        baseline_mask_indices = np.where(baseline_mask)[0]
        baseline_mask[baseline_mask_indices] &= eta_mask
        
        # Apply mask to events and jets
        events = events[eta_mask]
        jets = jets[eta_mask]
        
        # Step 5: Apply VBF criteria
        # Form all combinations of two jets
        jet_pairs = ak.combinations(jets, 2, fields=['jet1', 'jet2'])
        
        # Calculate delta_eta and invariant mass for each pair
        delta_eta = jet_pairs.jet1.eta - jet_pairs.jet2.eta
        mjj = (jet_pairs.jet1 + jet_pairs.jet2).mass
        
        # Apply VBF selection
        vbf_pair_mask = (abs(delta_eta) > 2.0) & (mjj > 500)
        vbf_event_mask = ak.any(vbf_pair_mask, axis=1)
        vbf_event_mask = ak.to_numpy(vbf_event_mask)  # Convert to NumPy array
        
        # Update baseline mask
        baseline_mask_indices = np.where(baseline_mask)[0]
        baseline_mask[baseline_mask_indices] &= vbf_event_mask
    else:
        # Initialize general mask
        general_mask = np.ones(n_events, dtype=bool)

        # Step 1: Process FatJets
        fatjets = events.FatJet
        # Create a mask for fatjets passing the criteria
        fatjet_mask = (fatjets.pt > 250) & (abs(fatjets.eta) < 2.5)
        # Count the number of fatjets passing the criteria per event
        n_fatjets = ak.sum(fatjet_mask, axis=1)
        # Events with at least one such fatjet
        has_fatjet = n_fatjets >= 1
        has_fatjet = ak.to_numpy(has_fatjet)

        # Update general mask
        general_mask &= has_fatjet

        # Step 2: Process AK4 jets
        jets = events.Jet
        # Create a mask for jets within eta range
        jet_eta_mask = abs(jets.eta) < 5.0
        # Apply mask to jets to get jets in eta range
        jets_in_eta_range = jets[jet_eta_mask]
        # Count the number of jets per event in eta range
        n_jets_in_eta_range = ak.num(jets_in_eta_range)
        # Events with at least 4 such jets
        has_four_jets = n_jets_in_eta_range >= 4
        has_four_jets = ak.to_numpy(has_four_jets)

        # Update general mask
        general_mask &= has_four_jets

        # Step 3: Combine baseline and general masks
        baseline_mask &= general_mask

        # Step 4: Form all combinations of two jets for events passing so far
        # Create a mask for the events passing so far
        selected_events_mask = baseline_mask

        # Select events and jets passing so far
        selected_events = events[selected_events_mask]
        selected_jets = jets_in_eta_range[selected_events_mask]

        # Form jet pairs
        jet_pairs = ak.combinations(selected_jets, 2, fields=['jet1', 'jet2'])

        # Calculate delta_eta and invariant mass for each pair
        delta_eta = jet_pairs.jet1.eta - jet_pairs.jet2.eta
        mjj = (jet_pairs.jet1 + jet_pairs.jet2).mass

        # Apply selection: any pair with delta_eta > 1.5 and mjj > 200
        pair_mask = (abs(delta_eta) > 1.5) & (mjj > 200)
        event_mask = ak.any(pair_mask, axis=1)
        event_mask = ak.to_numpy(event_mask)

        # Create a mask of length n_events, initialized to False
        vbf_event_mask = np.zeros(n_events, dtype=bool)
        # Set the events that passed the previous cuts and the VBF cuts to True
        vbf_event_mask[selected_events_mask] = event_mask

        # Update baseline mask
        baseline_mask &= vbf_event_mask
    assert len(baseline_mask) == len(events)
        
    
    # The baseline_mask now represents events that pass all cuts
    return baseline_mask

class TriggerEfficiencyProcessor(processor.ProcessorABC):
    def __init__(self, triggers, trig_vars, baseline_key = 'VBF'):
        self.triggers = triggers
        self.trig_vars = trig_vars
        self.baseline_key = baseline_key

        # Initialize output using accumulators
        self.output = dict_accumulator({
            trigger: dict_accumulator({
                var_name: dict_accumulator({
                    'total': list_accumulator(),
                    'pass': list_accumulator()
                }) for var_name in trig_vars
            }) for trigger in triggers
        })

    def process(self, events):
        output = self.output
        
        baseline = create_baseline_selection_mask(events, tag=self.baseline_key)

       # Compute variables using 'proc' functions from 'trig_vars'
        variables = {}
        for var_name, var_info in self.trig_vars.items():
            var_array = var_info['proc'](events)
            variables[var_name] = var_array
    
        # Create a valid mask where all variables are not NaN
        valid_vars_mask = np.ones(len(events), dtype=bool)
        for var_array in variables.values():
            valid_vars_mask &= ~np.isnan(ak.to_numpy(var_array))
    
        # Combined selections
        selection_total = baseline & valid_vars_mask
    
        # Loop through each trigger
        for trigger in self.triggers:
            trigger_mask = events.HLT[trigger]
    
            # Selection for passing events
            selection_pass = baseline & trigger_mask & valid_vars_mask
    
            # Collect per-event data for all variables
            for var_name in self.trig_vars:
                var_array = variables[var_name]
    
                var_total = var_array[selection_total]
                var_pass = var_array[selection_pass]
    
                # Append to accumulators
                output[trigger][var_name]['total'].extend(ak.to_numpy(var_total).tolist())
                output[trigger][var_name]['pass'].extend(ak.to_numpy(var_pass).tolist())
    
        return output

    def postprocess(self, accumulator):
        # Concatenate lists in the accumulator
        for trigger in accumulator:
            for var_name in accumulator[trigger]:
                accumulator[trigger][var_name]['total'] = np.array(accumulator[trigger][var_name]['total'])
                accumulator[trigger][var_name]['pass'] = np.array(accumulator[trigger][var_name]['pass'])
        return accumulator

class TriggerEfficiencyImprovementProcessor(processor.ProcessorABC):
    """
    Calculate the improvement in trigger efficiency after adding a logical OR with one of the VBF triggers, e.g.   QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
    """
    def __init__(self, triggers, trig_vars, vbf_trigger, baseline_key = 'VBF'):
        self.triggers = triggers
        self.trig_vars = trig_vars
        self.vbf_trigger = vbf_trigger
        self.baseline_key = baseline_key

        # Initialize output using accumulators
        self.output = dict_accumulator({
            trigger: dict_accumulator({
                var_name: dict_accumulator({
                    'total': list_accumulator(),
                    'pass': list_accumulator(),
                    'pass_or': list_accumulator() # events passing trigger or vbf_trigger
                }) for var_name in trig_vars
            }) for trigger in triggers
        })

    def process(self, events):
        output = self.output
        
        baseline = create_baseline_selection_mask(events, tag=self.baseline_key)

       # Compute variables using 'proc' functions from 'trig_vars'
        variables = {}
        for var_name, var_info in self.trig_vars.items():
            var_array = var_info['proc'](events)
            variables[var_name] = var_array
    
        # Create a valid mask where all variables are not NaN
        valid_vars_mask = np.ones(len(events), dtype=bool)
        for var_array in variables.values():
            valid_vars_mask &= ~np.isnan(ak.to_numpy(var_array))
    
        # Combined selections
        selection_total = baseline & valid_vars_mask
    
        # Loop through each trigger
        for trigger in self.triggers:
            trigger_mask = events.HLT[trigger]
            vbf_trigger_mask = events.HLT[self.vbf_trigger]
    
            # Selection for passing events
            selection_pass = baseline & trigger_mask & valid_vars_mask
            selection_pass_or = baseline & (trigger_mask | vbf_trigger_mask) & valid_vars_mask
    
            # Collect per-event data for all variables
            for var_name in self.trig_vars:
                var_array = variables[var_name]
    
                var_total = var_array[selection_total]
                var_pass = var_array[selection_pass]
                var_pass_or = var_array[selection_pass_or]
    
                # Append to accumulators
                output[trigger][var_name]['total'].extend(ak.to_numpy(var_total).tolist())
                output[trigger][var_name]['pass'].extend(ak.to_numpy(var_pass).tolist())
                output[trigger][var_name]['pass_or'].extend(ak.to_numpy(var_pass_or).tolist())
    
        return output

    def postprocess(self, accumulator):
        # Concatenate lists in the accumulator
        for trigger in accumulator:
            for var_name in accumulator[trigger]:
                accumulator[trigger][var_name]['total'] = np.array(accumulator[trigger][var_name]['total'])
                accumulator[trigger][var_name]['pass'] = np.array(accumulator[trigger][var_name]['pass'])
                accumulator[trigger][var_name]['pass_or'] = np.array(accumulator[trigger][var_name]['pass_or'])
        return accumulator

class VBFBaselineCutFlowProcessor(processor.ProcessorABC):
    def __init__(self):
        # Initialize the output accumulator
        self.output = dict_accumulator({
            'cutflow': dict_accumulator(),
            # Include other outputs if needed
        })

    def process(self, events):
        output = self.output
        cutflow = output['cutflow']

        # Total number of events before any cuts
        n_events = len(events)
        cutflow['All events'] = cutflow.get('All events', 0) + n_events

        ### Start of the selection steps ###

        # Initialize cumulative event mask
        cumulative_event_mask = np.ones(n_events, dtype=bool)

        # Step 1: Select jets with pT > 15 GeV
        jets = events.Jet[events.Jet.pt > 15]

        # Count events with at least one jet after pT > 15 GeV cut
        has_jets = ak.num(jets) >= 1
        n_events_jets_pt15 = ak.sum(has_jets)
        cutflow['pT > 15 GeV'] = cutflow.get('pT > 15 GeV', 0) + n_events_jets_pt15

        # Update cumulative event mask
        cumulative_event_mask = cumulative_event_mask & has_jets

        # Filter events and jets based on cumulative event mask
        events = events[cumulative_event_mask]
        jets = jets[cumulative_event_mask]

        # Sort jets by descending pT
        jets = jets[ak.argsort(jets.pt, ascending=False)]

        # Step 2: Ensure there are at least four jets
        has_four_jets = ak.num(jets) >= 4
        n_events_four_jets = ak.sum(has_four_jets)
        cutflow['At least 4 jets'] = cutflow.get('At least 4 jets', 0) + n_events_four_jets

        # Update cumulative event mask
        cumulative_event_mask = has_four_jets  # Already filtered events

        # Filter events and jets
        events = events[has_four_jets]
        jets = jets[has_four_jets]

        # Step 3: Apply pT thresholds to the leading jets
        leading_jet_pt_mask = (
            (jets[:, 0].pt > 105) &
            (jets[:, 1].pt > 88) &
            (jets[:, 2].pt > 76) &
            (jets[:, 3].pt > 15)
        )
        n_events_pt_thresholds = ak.sum(leading_jet_pt_mask)
        cutflow['Leading jet pT thresholds'] = cutflow.get('Leading jet pT thresholds', 0) + n_events_pt_thresholds

        # Update cumulative event mask
        cumulative_event_mask = leading_jet_pt_mask

        # Filter events and jets
        events = events[leading_jet_pt_mask]
        jets = jets[leading_jet_pt_mask]

        # Step 4: Apply eta ranges
        eta_mask = (
            (abs(jets[:, 0].eta) < 5.0) &
            (abs(jets[:, 1].eta) < 5.0) &
            (abs(jets[:, 2].eta) < 5.0) &
            (abs(jets[:, 3].eta) < 5.0)
        )
        n_events_eta = ak.sum(eta_mask)
        cutflow['Eta ranges'] = cutflow.get('Eta ranges', 0) + n_events_eta

        # Update cumulative event mask
        cumulative_event_mask = eta_mask

        # Filter events and jets
        events = events[eta_mask]
        jets = jets[eta_mask]

        # Step 5: Apply VBF criteria
        # Form all combinations of two jets
        jet_pairs = ak.combinations(jets, 2, fields=['jet1', 'jet2'])

        # Calculate delta_eta and invariant mass for each pair
        delta_eta = jet_pairs.jet1.eta - jet_pairs.jet2.eta
        mjj = (jet_pairs.jet1 + jet_pairs.jet2).mass

        # Apply VBF selection
        vbf_pair_mask = (abs(delta_eta) > 2.0) & (mjj > 500)
        vbf_event_mask = ak.any(vbf_pair_mask, axis=1)
        n_events_vbf = ak.sum(vbf_event_mask)
        cutflow['delta_eta and mjj'] = cutflow.get('delta_eta and mjj', 0) + n_events_vbf

        # Update cumulative event mask
        cumulative_event_mask = vbf_event_mask

        # Filter events and jets
        events = events[vbf_event_mask]
        jets = jets[vbf_event_mask]

        return output

    def postprocess(self, accumulator):
        # No postprocessing required in this case
        return accumulator

class GeneralBaselineCutFlowProcessor(processor.ProcessorABC):
    def __init__(self):
        # Initialize the output accumulator
        self.output = dict_accumulator({
            'cutflow': dict_accumulator(),
            # Include other outputs if needed
        })

    def process(self, events):
        output = self.output
        cutflow = output['cutflow']

        # Total number of events before any cuts
        n_events = len(events)
        cutflow['All events'] = cutflow.get('All events', 0) + n_events

        ### Start of the selection steps ###

        # Initialize cumulative event mask as an Awkward Array
        cumulative_event_mask = ak.from_numpy(np.ones(n_events, dtype=bool))

        # Step 1: Select FatJets with pt > 250 GeV and abs(eta) < 2.5
        fatjets_mask = (events.FatJet.pt > 250) & (abs(events.FatJet.eta) < 2.5)

        # Require at least one such fatjet per event
        has_fatjet = ak.sum(fatjets_mask, axis=-1) >= 1

        # Update cumulative event mask
        cumulative_event_mask = cumulative_event_mask & has_fatjet
        
        n_events_has_fatjet = ak.sum(cumulative_event_mask)
        cutflow['At least 1 FatJet with pt > 250 GeV and abs(eta) < 2.5'] = (
            cutflow.get('At least 1 FatJet with pt > 250 GeV and abs(eta) < 2.5', 0) + n_events_has_fatjet
        )

        

        # Step 2: Select AK4 jets with abs(eta) < 5
        jets_mask = abs(events.Jet.eta) < 5.0

        # Require at least 4 such AK4 jets per event
        has_four_jets = ak.sum(jets_mask, axis=-1) >= 4

        # Update cumulative event mask
        cumulative_event_mask = cumulative_event_mask & has_four_jets
        
        n_events_four_jets = ak.sum(cumulative_event_mask)
        cutflow['At least 4 AK4 jets with abs(eta) < 5'] = (
            cutflow.get('At least 4 AK4 jets with abs(eta) < 5', 0) + n_events_four_jets
        )

        # Step 3: Form all combinations of two jets (without filtering)
        jet_pairs = ak.combinations(events.Jet, 2, fields=['jet1', 'jet2'])

        # Create masks for jet pairs based on individual jet selections
        jet1_mask = abs(jet_pairs.jet1.eta) < 5.0
        jet2_mask = abs(jet_pairs.jet2.eta) < 5.0
        jet_pair_mask = jet1_mask & jet2_mask

        # Calculate delta_eta and invariant mass mjj for each pair
        delta_eta = jet_pairs.jet1.eta - jet_pairs.jet2.eta
        mjj = (jet_pairs.jet1 + jet_pairs.jet2).mass

        # Apply selection: any pair with delta_eta > 1.5 and mjj > 200
        selection_mask = jet_pair_mask & (abs(delta_eta) > 1.5) & (mjj > 200)
        event_has_good_pair = ak.any(selection_mask, axis=1)

        # Update cumulative event mask
        cumulative_event_mask = cumulative_event_mask & event_has_good_pair
        
        n_events_pairs = ak.sum(cumulative_event_mask)
        cutflow['Any jet pair with delta_eta > 1.5 and mjj > 200'] = (
            cutflow.get('Any jet pair with delta_eta > 1.5 and mjj > 200', 0) + n_events_pairs
        )

        

        # Note: Events and objects are not filtered at any step to ensure consistent shapes

        # Optionally, store the cumulative mask in the output if needed
        # output['cumulative_event_mask'] = cumulative_event_mask

        return output

    def postprocess(self, accumulator):
        # No postprocessing required in this case
        return accumulator

class TriggerSoupProcessor(processor.ProcessorABC):
    """
    Calculate the trigger efficiency of the soup without each trigger in the soup.
    """
    def __init__(self, trig_soup, trig_vars, baseline_key='VBF'):
        self.trig_soup = trig_soup
        self.trig_vars = trig_vars
        self.baseline_key = baseline_key

        # Generate trigger combinations
        self.trigger_combinations = []

        # All triggers (no triggers omitted)
        self.trigger_combinations.append(list(self.trig_soup))
# 
        # Now, for each trigger, create a combination without that trigger
        for idx, trig in enumerate(self.trig_soup):
            triggers_without_current = self.trig_soup[:idx] + self.trig_soup[idx+1:]
            self.trigger_combinations.append(triggers_without_current)

        # Create combination names
        self.combination_names = ['All_triggers'] + ['Without_{}'.format(trig) for trig in self.trig_soup]

        # Print trigger combinations for debugging
        print(f"Trigger combinations: {list(zip(self.combination_names, self.trigger_combinations))}")

        # Initialize output using accumulators
        # Add 'Baseline' key to store 'total' arrays
        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                var_name: list_accumulator()
                for var_name in trig_vars
            }),
            **{
                name: dict_accumulator({
                    var_name: dict_accumulator({
                        'pass': list_accumulator(),
                    }) for var_name in trig_vars
                }) for name in self.combination_names
            }
        })

    def process(self, events):
        output = self.output

        baseline = create_baseline_selection_mask(events, tag=self.baseline_key)

        # Compute variables using 'proc' functions from 'trig_vars'
        variables = {}
        for var_name, var_info in self.trig_vars.items():
            var_array = var_info['proc'](events)
            variables[var_name] = var_array

        # Create a valid mask where all variables are not NaN
        valid_vars_mask = np.ones(len(events), dtype=bool)
        for var_array in variables.values():
            valid_vars_mask &= ~np.isnan(ak.to_numpy(var_array))

        # Combined selection for total (baseline without triggers)
        selection_total = baseline & valid_vars_mask

        # Collect 'total' arrays under 'Baseline' key
        for var_name in self.trig_vars:
            var_array = variables[var_name]
            var_total = var_array[selection_total]
            output['Baseline'][var_name].extend(ak.to_numpy(var_total).tolist())

        # Loop through each trigger combination
        for name, triggers in zip(self.combination_names, self.trigger_combinations):
            # Create the trigger mask for this combination
            if triggers:
                trigger_masks = [events.HLT[trig] for trig in triggers]
                combined_trigger_mask = ak.zeros_like(events.HLT[triggers[0]], dtype=bool)
                for mask in trigger_masks:
                    combined_trigger_mask = combined_trigger_mask | mask
            else:
                # If triggers is empty, combined mask is all False
                combined_trigger_mask = ak.zeros_like(events.HLT[self.trig_soup[0]], dtype=bool)

            # Selection for passing events
            selection_pass = baseline & combined_trigger_mask & valid_vars_mask

            # Collect per-event data for all variables
            for var_name in self.trig_vars:
                var_array = variables[var_name]

                var_pass = var_array[selection_pass]

                # Append to accumulators
                output[name][var_name]['pass'].extend(ak.to_numpy(var_pass).tolist())

        # print(f"Baseline mask (total events: {len(baseline)}): {np.sum(baseline)} events pass baseline selection")
        # for trig in triggers:
        #     print(f"Trigger: {trig}, Events passing: {np.sum(events.HLT[trig])}")
        # print(f"Combined trigger mask (name: {name}, total events: {len(combined_trigger_mask)}): {np.sum(combined_trigger_mask)} events pass")

        return output

    def postprocess(self, accumulator):
        # Concatenate lists in the accumulator
        for key in accumulator:
            if key == 'Baseline':
                for var_name in accumulator[key]:
                    accumulator[key][var_name] = np.array(accumulator[key][var_name])
            else:
                for var_name in accumulator[key]:
                    accumulator[key][var_name]['pass'] = np.array(accumulator[key][var_name]['pass'])
        return accumulator
