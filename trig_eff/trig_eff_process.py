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

def create_baseline_selection_mask(events):
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
    events = events[baseline_mask]
    jets = jets[baseline_mask]
    
    # Sort jets by descending pT
    jets = jets[ak.argsort(jets.pt, ascending=False)]
    
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
    
    # The baseline_mask now represents events that pass all cuts
    return baseline_mask

class TriggerEfficiencyProcessor(processor.ProcessorABC):
    def __init__(self, triggers, trig_vars):
        self.triggers = triggers
        self.trig_vars = trig_vars

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
        
        baseline = create_baseline_selection_mask(events)

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
