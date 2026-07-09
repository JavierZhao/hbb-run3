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

# JetID + jet-veto-map helpers (for v15 / 2024 NanoAOD that often lack jetId).
from corrections import JET_VETOMAP_KEYS, correct_jetid, event_passes_jet_veto_map


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
    "2024": "2024_Summer24",
}

# Years where NanoAOD often lacks jetId; use POG jetid.json.gz via correctionlib.
_JETID_CORRECTION_YEARS = {"2023", "2023BPix", "2024"}

_JETID_FALLBACK_WARNED = False

# Branch-set detection cache: which tagger we're reading on the current worker.
# Set on first call, kept silent in subsequent calls so the log isn't spammy.
_TAGGER_SOURCE_CACHE = None  # one of: 'globalParT3', 'particleNet_v12', 'pnetVsQCD_legacy', 'none'


def _detect_tagger_source(fields):
    """Return a label for which fatjet tagger branches are available."""
    f = set(fields)
    if {'globalParT3_Xbb', 'globalParT3_Xcc', 'globalParT3_QCD'}.issubset(f):
        return 'globalParT3'
    if {'particleNet_XbbVsQCD', 'particleNet_XccVsQCD', 'particleNet_QCD'}.issubset(f):
        return 'particleNet_v12'
    if 'particleNet_XbbVsQCD' in f:
        return 'pnetVsQCD_legacy'  # only Xbb side; cc not separable
    return 'none'


def ak_clip(arr, min_value, max_value):
    """Clip awkward array values to the valid range expected by correctionlib payloads."""
    return ak.where(arr < min_value, min_value, ak.where(arr > max_value, max_value, arr))


def binomial_eff_err(k, n):
    """
    Per-bin statistical uncertainty on a binomial efficiency eff = k/n.
    Uses the Gaussian approximation: sigma(eff) = sqrt(eff*(1-eff)/n).

    Returns NaN where n == 0. This is the standard "sqrt(N) for binomial counts" form
    and is more accurate than sqrt(1/k) when eff is far from 0; for very low eff it
    reduces to the same numerical value.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.where(n > 0, k / n, np.nan)
        err = np.where(n > 0, np.sqrt(np.clip(eff * (1.0 - eff), 0.0, None) / n), np.nan)
    return eff, err


def sf_err_from_counts(k_pass_d, n_base_d, k_pass_m, n_base_m):
    """
    Statistical uncertainty on SF = eff_data / eff_mc, propagated assuming independence.
    sigma(SF) = SF * sqrt((1-eff_d)/k_pass_d + (1-eff_m)/k_pass_m).

    Returns (sf, sf_err) with NaN where any input bin has zero pass or baseline counts.
    """
    eff_d, _ = binomial_eff_err(k_pass_d, n_base_d)
    eff_m, _ = binomial_eff_err(k_pass_m, n_base_m)
    with np.errstate(divide='ignore', invalid='ignore'):
        sf  = np.where((eff_m > 0) & np.isfinite(eff_d), eff_d / eff_m, np.nan)
        rel = np.where(np.asarray(k_pass_d) > 0, (1.0 - eff_d) / np.asarray(k_pass_d), np.nan) \
            + np.where(np.asarray(k_pass_m) > 0, (1.0 - eff_m) / np.asarray(k_pass_m), np.nan)
        sf_err = np.where(np.isfinite(sf) & np.isfinite(rel), np.abs(sf) * np.sqrt(np.clip(rel, 0.0, None)), np.nan)
    return sf, sf_err


def get_pog_json(obj, year):
    """Return path to a POG json file for a given object/year."""
    group, filename = POG_JSONS[obj]
    year_key = POG_YEARS[year]
    return f"{POG_CORRECTION_PATH}/POG/{group}/{year_key}/{filename}"


def get_txbb_txcc_score(jets):
    """
    AN Eq. 4 TXbbcc = (Xbb + Xcc) / (Xbb + Xcc + QCD).

    Branch resolution (in priority):
      1. globalParT3 raw probabilities (newer Run 3 NanoAOD / private v14+):
            (globalParT3_Xbb + globalParT3_Xcc) / (Xbb + Xcc + globalParT3_QCD)
      2. particleNet NanoAODv12 normalized scores, mirroring main-hbb-run3 objects.py:
            (particleNet_XbbVsQCD + particleNet_XccVsQCD)
            / (particleNet_XbbVsQCD + particleNet_XccVsQCD + particleNet_QCD)
      3. Legacy: just particleNet_XbbVsQCD (single-flavor proxy).
    """
    global _TAGGER_SOURCE_CACHE
    if _TAGGER_SOURCE_CACHE is None:
        _TAGGER_SOURCE_CACHE = _detect_tagger_source(jets.fields)
        # one-time per worker, terse:
        print(f"[txbb tagger source] {_TAGGER_SOURCE_CACHE}")
    src = _TAGGER_SOURCE_CACHE

    if src == 'globalParT3':
        xbb = ak.fill_none(jets.globalParT3_Xbb, 0.0)
        xcc = ak.fill_none(jets.globalParT3_Xcc, 0.0)
        qcd = ak.fill_none(jets.globalParT3_QCD, 0.0)
        return (xbb + xcc) / (xbb + xcc + qcd + 1e-12)

    if src == 'particleNet_v12':
        xbb_vs = ak.fill_none(jets.particleNet_XbbVsQCD, 0.0)
        xcc_vs = ak.fill_none(jets.particleNet_XccVsQCD, 0.0)
        qcd    = ak.fill_none(jets.particleNet_QCD, 0.0)
        return (xbb_vs + xcc_vs) / (xbb_vs + xcc_vs + qcd + 1e-12)

    if src == 'pnetVsQCD_legacy':
        return ak.fill_none(jets.particleNet_XbbVsQCD, 0.0)

    return ak.zeros_like(jets.pt)


def get_xbb_score(jets):
    """Single-flavour Xbb score for the pass-bb vs pass-cc comparison."""
    f = set(jets.fields)
    if 'globalParT3_Xbb' in f:
        return ak.fill_none(jets.globalParT3_Xbb, 0.0)
    if 'particleNet_XbbVsQCD' in f:
        return ak.fill_none(jets.particleNet_XbbVsQCD, 0.0)
    return ak.zeros_like(jets.pt)


def get_xcc_score(jets):
    """Single-flavour Xcc score for the pass-bb vs pass-cc comparison."""
    f = set(jets.fields)
    if 'globalParT3_Xcc' in f:
        return ak.fill_none(jets.globalParT3_Xcc, 0.0)
    if 'particleNet_XccVsQCD' in f:
        return ak.fill_none(jets.particleNet_XccVsQCD, 0.0)
    return ak.zeros_like(jets.pt)


def _fallback_jetid_tight_mask(jets):
    """NanoAOD jetId bitmask when the branch exists (pre-v15 style).

    v15 / 2024 NanoAOD often lacks the jetId branch entirely. In that case the
    only correct behavior is to fail loudly, since silently passing every jet
    would bias the trigger probe. Callers must have a correctionlib evaluator
    available for years where jetId is absent.
    """
    if "jetId" in jets.fields:
        return jets.isTight
    raise RuntimeError(
        "Cannot compute tight JetID: no jetId branch and correctionlib evaluation failed."
    )


def _evaluate_jetid_tight(evaluator, corr_name, flat_jets):
    """Evaluate POG tight JetID; scalar-loop fallback for older correctionlib schemas."""
    def real_col(x):
        x = ak.fill_none(x, 0.0)
        a = np.asarray(ak.to_numpy(x), dtype=np.float64, order="C")
        np.nan_to_num(a, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        return a

    def int_col(x):
        x = ak.fill_none(x, 0)
        return np.asarray(ak.to_numpy(x), dtype=np.int64, order="C")

    eta = np.clip(real_col(flat_jets.eta), -4.7, 4.7)
    chHEF = np.clip(real_col(flat_jets.chHEF), 0.0, 1.0)
    neHEF = np.clip(real_col(flat_jets.neHEF), 0.0, 1.0)
    chEmEF = np.clip(real_col(flat_jets.chEmEF), 0.0, 1.0)
    neEmEF = np.clip(real_col(flat_jets.neEmEF), 0.0, 1.0)
    muEF = np.clip(real_col(flat_jets.muEF), 0.0, 1.0)
    chMult = np.clip(int_col(flat_jets.chMultiplicity), 0, 100)
    neMult = np.clip(int_col(flat_jets.neMultiplicity), 0, 100)
    nConst = chMult + neMult

    if nConst.size == 0:
        return np.array([], dtype=bool)

    corr = evaluator[corr_name]
    try:
        mask_flat = np.asarray(
            corr.evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, nConst),
            dtype=bool,
        )
    except RuntimeError:
        mask_flat = np.empty(nConst.size, dtype=bool)
        for i in range(nConst.size):
            mask_flat[i] = bool(
                corr.evaluate(
                    float(eta[i]),
                    float(chHEF[i]),
                    float(neHEF[i]),
                    float(chEmEF[i]),
                    float(neEmEF[i]),
                    float(muEF[i]),
                    int(chMult[i]),
                    int(neMult[i]),
                    int(nConst[i]),
                )
            )
    return mask_flat


def corrected_jetid_tight_mask(jets, jet_type, year, use_jetid_correction=True):
    """
    For v14+ / 2024 NanoAOD, compute Tight JetID from correctionlib (jetId often absent).
    Falls back to NanoAOD isTight only when the jetId branch is present.
    """
    global _JETID_FALLBACK_WARNED

    if not use_jetid_correction or year not in _JETID_CORRECTION_YEARS:
        return _fallback_jetid_tight_mask(jets)

    if jet_type == "AK8":
        corr_name = "AK8PUPPI_Tight"
    elif jet_type == "AK4":
        corr_name = "AK4PUPPI_Tight"
    else:
        return _fallback_jetid_tight_mask(jets)

    required_fields = {
        "eta", "chHEF", "neHEF", "chEmEF", "neEmEF", "muEF",
        "chMultiplicity", "neMultiplicity",
    }
    if not required_fields.issubset(set(jets.fields)):
        if not _JETID_FALLBACK_WARNED:
            print("Warning: missing JetID input branches; cannot apply correctionlib JetID.")
            _JETID_FALLBACK_WARNED = True
        return _fallback_jetid_tight_mask(jets)

    try:
        import correctionlib  # local import to keep script usable when unavailable

        evaluator = correctionlib.CorrectionSet.from_file(get_pog_json("jetid", year))
        flat_jets = ak.flatten(jets)
        njets = ak.num(jets)
        mask_flat = _evaluate_jetid_tight(evaluator, corr_name, flat_jets)
        return ak.values_astype(ak.unflatten(mask_flat, njets), bool)
    except Exception as exc:
        if not _JETID_FALLBACK_WARNED:
            print(f"Warning: JetID correction unavailable ({exc}).")
            _JETID_FALLBACK_WARNED = True
        return _fallback_jetid_tight_mask(jets)


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


# 2D trigger variables configuration.
# Bin edges are coarser than uniform 15x20 to keep per-bin stats high
# (target sigma_SF ~5-10%): finer in the populated band, coarser in the tails.
# Variable axes are drop-in compatible with downstream np.histogram2d usage
# via .edges.
# pT edges start at 450 (matches the new pT > 450 baseline cut — see
# create_baseline_selection_mask_new). Old edges 200..360 are now always empty.
_msd_edges = [0, 30, 60, 90, 120, 160, 220, 300]
_pt_edges  = [450, 525, 600, 700, 850, 1000]
trig_vars_2d = {
    'pt_vs_msd': {
        'label_x': "Leading Jet $m_{SD}$ [GeV]",
        'label_y': "Leading Jet $p_{T}$ [GeV]",
        'axis_x': hist.axis.Variable(_msd_edges, name="msd", label="Leading Jet $m_{SD}$ [GeV]"),
        'axis_y': hist.axis.Variable(_pt_edges,  name="pt",  label="Leading Jet $p_{T}$ [GeV]"),
        'proc_x': lambda events: ak.fill_none(ak.pad_none(events.FatJet.msoftdrop, 1, clip=True)[:, 0], np.nan),
        'proc_y': lambda events: ak.fill_none(ak.pad_none(events.FatJet.pt, 1, clip=True)[:, 0], np.nan)
    }
}

trigger_dict_periods = {
    # All years extended on 2026-05-22 with PFJet500, AK8PFJet425_SoftDropMass40,
    # and PFHT1050 (the "augmented OR"). study_extra_triggers.py q80 showed these
    # close out the residual SF gap; study_signal_eff.py q10 confirmed +2-5% on
    # ggF/VBF acceptance. scan_2022 / 2022EE / 2023 / 2023BPix verified all 5
    # triggers exist in every readable Muon NanoAOD file in each year (2023
    # restricted to v3+ in muon_keys_by_period to avoid the PNetBB0p06-missing
    # v1/v2 re-recos).
    '2022': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65',
        'PFJet500',
        'AK8PFJet425_SoftDropMass40',
        'PFHT1050',
    ],
    '2022EE': [
        'AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35',
        'QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65',
        'PFJet500',
        'AK8PFJet425_SoftDropMass40',
        'PFHT1050',
    ],
    '2023': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'PFJet500',
        'AK8PFJet425_SoftDropMass40',
        'PFHT1050',
        # VBF_DiPFJet125_45_Mjj720_Detajj3p0 uses a different calibration strategy.
        # Excluded from trigger SF OR and per-trigger table in this workflow.
    ],
    '2023BPix': [
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'PFJet500',
        'AK8PFJet425_SoftDropMass40',
        'PFHT1050',
        # VBF_DiPFJet125_45_Mjj720_Detajj3p0 uses a different calibration strategy.
        # It is also in ParkingVBF stream only and not stored in Muon0 PD for Run2023D.
    ],
    '2024': [
        # 2024 NanoAODv15 HLT menu: 4 of the 5 names match 2023/2023BPix.
        # AK8PFJet425_SoftDropMass40 is NOT in the 2024 menu; the equivalent
        # leg is AK8PFJet425_SoftDropMass30 (same pT threshold, looser mSD
        # cut). Verified via probe_2024_ak8_triggers.py on 2026-05-29 across
        # Muon0_Run2024C/F/H/I, TTtoLNu2Q, and QCD_MuEnriched_PT_300to470.
        'AK8PFJet250_SoftDropMass40_PNetBB0p06',
        'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
        'PFJet500',
        'AK8PFJet425_SoftDropMass30',
        'PFHT1050',
    ],
}

# Muon reference triggers for SF calculation
muon_triggers = {
    '2022': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2022EE': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2023BPix': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
    '2024': ["HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100", "HLT_IsoMu24"],
}

def create_baseline_selection_mask_new(events, tag, year, txbb_region='inclusive', use_jetid_correction=True, require_muon=True):
    """
    txbb_region: 'fail'      -> candidatejet (Xbb+Xcc) <  TXBB_WP (0.82) [AN Sec 5.1 Fail]
                 'pass_bb'   -> Xbb+Xcc >= 0.82 AND Xbb > Xcc           [AN Sec 5.1 Pass-bb]
                 'pass_cc'   -> Xbb+Xcc >= 0.82 AND Xcc > Xbb           [AN Sec 5.1 Pass-cc]
                 'pass'      -> Xbb+Xcc >= 0.82 (no bb/cc split, back-compat)
                 'fail_legacy' / 'inclusive' -> no region cut

    Also enforces the analysis baseline on the chosen candidate jet:
        candidatejet.pt > 300 GeV  AND  candidatejet.msoftdrop > 40 GeV
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

    # Per-event candidate observables (NaN/-1 where no candidate jet exists)
    candidate_score = ak.fill_none(get_txbb_txcc_score(candidatejet), -1)
    candidate_xbb   = ak.fill_none(get_xbb_score(candidatejet), -1)
    candidate_xcc   = ak.fill_none(get_xcc_score(candidatejet), -1)
    candidate_pt    = ak.fill_none(candidatejet.pt, -1)
    candidate_msd   = ak.fill_none(candidatejet.msoftdrop, -1)

    # Step 2a: analysis-baseline kinematic cuts on the candidate jet (AN signal-region preselection).
    # pT > 450 (was 300): puts the candidate inside the plateau of every trigger leg in
    # the augmented OR (AK8PFJet250 + AK8PFJet425 + PFJet500 etc.) so the per-bin SF is
    # not biased by the AK8 turn-on shape. q10 study on 2026-05-22 showed augmented_OR SF
    # rises 0.73->0.97 when this cut is applied — see study_extra_triggers.py results.
    kine_mask = (candidate_pt > 450) & (candidate_msd > 40)
    baseline_mask &= ak.to_numpy(kine_mask)

    # Step 2b: AN Section 5.1 region split on the candidate jet
    pass_overall = candidate_score >= TXBB_WP
    if txbb_region == 'pass_bb':
        baseline_mask &= ak.to_numpy(pass_overall & (candidate_xbb > candidate_xcc))
    elif txbb_region == 'pass_cc':
        baseline_mask &= ak.to_numpy(pass_overall & (candidate_xcc > candidate_xbb))
    elif txbb_region == 'pass':           # back-compat: combined pass region
        baseline_mask &= ak.to_numpy(pass_overall)
    elif txbb_region == 'fail':
        baseline_mask &= ak.to_numpy(~pass_overall)
    # 'inclusive' / 'fail_legacy' / anything else -> no extra cut

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

    # Step 4: Muon veto - require at least one loose muon far from leading FatJet.
    # Only relevant for the SF probe (Z(μμ)-orthogonal data + QCD/ttbar MC). For
    # signal MC (ggF/VBF/WH/ZH/ttH) we don't want this — pass require_muon=False
    # to skip it and measure trigger efficiency directly on the analysis-region
    # selection.
    if require_muon:
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

        # Initialize output using accumulators.
        # 'score' / 'pass_score' carry the leading-FatJet (txbb+txcc) per event,
        # used downstream for the AN Fig 9 efficiency-vs-tagger-score plot.
        self.output = dict_accumulator({
            'Baseline': dict_accumulator({
                var_name: dict_accumulator({
                    'x': list_accumulator(),
                    'y': list_accumulator(),
                    'score': list_accumulator(),
                })
                for var_name in trig_vars_2d
            }),
            'Numerator': dict_accumulator({
                var_name: dict_accumulator({
                    'pass_x': list_accumulator(),
                    'pass_y': list_accumulator(),
                    'pass_score': list_accumulator(),
                })
                for var_name in trig_vars_2d
            }),
            'PerTrigger': dict_accumulator({
                trigger: dict_accumulator({
                    var_name: dict_accumulator({
                        'pass_x': list_accumulator(),
                        'pass_y': list_accumulator(),
                        'pass_score': list_accumulator(),
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

        # Compute leading-FatJet (txbb+txcc) score per event for the AN Fig 9
        # efficiency-vs-tagger-score projection. Mirrors the candidate-jet
        # selection used inside create_baseline_selection_mask_new (without
        # the txbb pass/fail cut, since that is what we plot against).
        fatjets = events.FatJet
        fatjet_tight = corrected_jetid_tight_mask(
            fatjets, jet_type="AK8", year=self.year,
            use_jetid_correction=self.use_jetid_correction,
        )
        cand = fatjets[
            (fatjets.pt > 200) & (abs(fatjets.eta) < 2.5) & fatjet_tight
        ][:, :2]
        cand_scores = get_txbb_txcc_score(cand)
        cand_lead = ak.firsts(cand[ak.argmax(cand_scores, axis=1, keepdims=True)])
        score_arr = ak.to_numpy(ak.fill_none(get_txbb_txcc_score(cand_lead), np.nan))

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
            output['Baseline'][var_name]['score'].extend(score_arr[selection_baseline].tolist())

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
            output['Numerator'][var_name]['pass_score'].extend(score_arr[selection_pass].tolist())

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
                    output['PerTrigger'][trigger][var_name]['pass_score'].extend(
                        score_arr[selection_pass_single].tolist()
                    )

        return output

    def postprocess(self, accumulator):
        # Convert lists to arrays
        for key in accumulator:
            if key == 'Baseline':
                for var_name in accumulator[key]:
                    for sub in ('x', 'y', 'score'):
                        if sub in accumulator[key][var_name]:
                            accumulator[key][var_name][sub] = np.array(accumulator[key][var_name][sub])
            elif key == 'Numerator':
                for var_name in accumulator[key]:
                    for sub in ('pass_x', 'pass_y', 'pass_score'):
                        if sub in accumulator[key][var_name]:
                            accumulator[key][var_name][sub] = np.array(accumulator[key][var_name][sub])
            elif key == 'PerTrigger':
                for trigger in accumulator[key]:
                    for var_name in accumulator[key][trigger]:
                        for sub in ('pass_x', 'pass_y', 'pass_score'):
                            if sub in accumulator[key][trigger][var_name]:
                                accumulator[key][trigger][var_name][sub] = np.array(
                                    accumulator[key][trigger][var_name][sub]
                                )
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

        # Per-bin statistical uncertainty: binomial Gaussian on each efficiency,
        # propagated to SF assuming independence.
        _, sf_unc = sf_err_from_counts(
            hist_pass_data, hist_baseline_data,
            hist_pass_mc,   hist_baseline_mc,
        )
        sf_unc[~valid_sf] = np.nan

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

        # Companion figure: 2-panel (SF | SF statistical uncertainty)
        fig_unc, (ax_sf, ax_un) = plt.subplots(1, 2, figsize=(16, 7))

        im_sf = ax_sf.pcolormesh(
            xedges, yedges, sf.T,
            cmap='viridis', vmin=0, vmax=1.2, shading='flat',
        )
        cb_sf = plt.colorbar(im_sf, ax=ax_sf)
        cb_sf.set_label('Trigger efficiency scale factor', fontsize=12)
        ax_sf.set_xlabel(var_info['label_x'], fontsize=12)
        ax_sf.set_ylabel(var_info['label_y'], fontsize=12)
        ax_sf.set_title(f"Scale Factor\nOverall: {total_sf:.3f}", fontsize=12)
        hep.cms.label(ax=ax_sf, data=True, year=year, com="13.6", fontsize=10)

        im_un = ax_un.pcolormesh(
            xedges, yedges, sf_unc.T,
            cmap='viridis', vmin=0, vmax=0.1, shading='flat',
        )
        cb_un = plt.colorbar(im_un, ax=ax_un)
        cb_un.set_label('Trigger efficiency scale factor uncertainty', fontsize=12)
        ax_un.set_xlabel(var_info['label_x'], fontsize=12)
        ax_un.set_ylabel(var_info['label_y'], fontsize=12)
        ax_un.set_title("Scale Factor Statistical Uncertainty", fontsize=12)
        hep.cms.label(ax=ax_un, data=True, year=year, com="13.6", fontsize=10)

        plt.tight_layout()
        unc_path = os.path.join(save_dir, f"scale_factor_{var_name}_2d_unc.png")
        plt.savefig(unc_path, dpi=200, bbox_inches='tight')
        if show_img:
            plt.show()
        else:
            print(f"Saved {unc_path}")
        plt.close(fig_unc)


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
        """Return (bin_centers, efficiency, eff_err) arrays for one trigger key and one axis."""
        if key == 'OR':
            pass_x = np.array(output['Numerator'][var_name]['pass_x'])
            pass_y = np.array(output['Numerator'][var_name]['pass_y'])
        else:
            if 'PerTrigger' not in output or key not in output['PerTrigger']:
                return None, None, None
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

        eff, err = binomial_eff_err(n_pass, n_base)
        centers = 0.5 * (bins[:-1] + bins[1:])
        return centers, eff, err

    def _draw_panel(ax, output, is_data, axis, bins, plateau_mask_fn=None):
        for i, (key, label) in enumerate(zip(all_keys, all_labels)):
            centers, eff, err = _eff_arrays(output, key, axis, bins, plateau_mask_fn)
            if centers is None:
                continue
            mk  = '+' if label == 'Soup' else 'o'
            mks = 8   if label == 'Soup' else 4
            ax.errorbar(
                centers, eff, yerr=err, label=label, color=colors[i],
                linestyle='-', marker=mk, markersize=mks, linewidth=2,
                elinewidth=1, capsize=2,
            )
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


def plot_baseline_kinematic_distributions(output_data, output_mc, save_dir=None,
                                          year='2022', region='', show_img=False):
    """
    Three-panel data vs MC overlay (with Data/MC ratio) of the baseline-jet kinematics
    after the analysis preselection (pT > 300, mSD > 40) and the chosen AN region cut:
        - leading FatJet pT
        - leading FatJet mSD
        - leading FatJet (Xbb + Xcc) tagger score

    These are the per-event arrays already stored in the Baseline accumulator
    (no trigger requirement); used to validate that the region split lands
    where expected (matching the AN Fig 4 spirit).
    """
    if save_dir is None:
        save_dir = "./figures_sf_2d"
    os.makedirs(save_dir, exist_ok=True)

    var_name = 'pt_vs_msd'
    if var_name not in output_data.get('Baseline', {}):
        return

    pt_d   = np.asarray(output_data['Baseline'][var_name].get('y', []))
    msd_d  = np.asarray(output_data['Baseline'][var_name].get('x', []))
    sc_d   = np.asarray(output_data['Baseline'][var_name].get('score', []))
    pt_m   = np.asarray(output_mc['Baseline'][var_name].get('y', []))
    msd_m  = np.asarray(output_mc['Baseline'][var_name].get('x', []))
    sc_m   = np.asarray(output_mc['Baseline'][var_name].get('score', []))

    if pt_d.size == 0 or pt_m.size == 0:
        print(f"  plot_baseline_kinematic_distributions: empty data or MC; skipping.")
        return

    # MC normalised to data integral so shape comparison is meaningful regardless of cross-section
    mc_norm = pt_d.size / max(pt_m.size, 1)

    panels = [
        ('pt',  pt_d,  pt_m,  np.linspace(300, 1500, 25), "Leading jet $p_T$ [GeV]"),
        ('msd', msd_d, msd_m, np.linspace(40,  300, 27),  "Leading jet $m_{SD}$ [GeV]"),
        ('score', sc_d, sc_m, np.linspace(0,   1,   21),  "Leading jet $X_{bb} + X_{cc}$"),
    ]

    fig = plt.figure(figsize=(21, 7))
    for i, (vname, xd, xm, bins, xlabel) in enumerate(panels):
        # Two-row gridspec for histogram + ratio
        gs = fig.add_gridspec(2, 3, height_ratios=[3, 1], hspace=0.05)
        ax  = fig.add_subplot(gs[0, i])
        rax = fig.add_subplot(gs[1, i], sharex=ax)

        nd, edges = np.histogram(xd, bins=bins)
        nm, _     = np.histogram(xm, bins=bins)
        nm_scaled = nm * mc_norm
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths  = np.diff(edges)

        ax.errorbar(centers, nd, yerr=np.sqrt(nd), fmt='o', color='black',
                    label='Data', markersize=4, capsize=2)
        ax.bar(centers, nm_scaled, width=widths, color='tab:red', alpha=0.4,
               label='MC (norm. to data)', edgecolor='tab:red')
        ax.set_yscale('log')
        ax.set_ylabel('Events')
        ax.legend(fontsize=10, loc='upper right')
        ax.set_title(xlabel, fontsize=11)
        plt.setp(ax.get_xticklabels(), visible=False)
        hep.cms.label(ax=ax, data=True, year=year, com="13.6", fontsize=10,
                      label=f"Preliminary, region={region}")

        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(nm_scaled > 0, nd / nm_scaled, np.nan)
            ratio_err = np.where((nm_scaled > 0) & (nd > 0),
                                 ratio * np.sqrt(1.0 / np.clip(nd, 1, None) + 1.0 / np.clip(nm, 1, None)),
                                 np.nan)
        rax.errorbar(centers, ratio, yerr=ratio_err, fmt='o', color='black',
                     markersize=4, capsize=2)
        rax.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        rax.set_ylim(0, 2)
        rax.set_xlabel(xlabel)
        rax.set_ylabel('Data / MC')

    plt.tight_layout()
    save_path = os.path.join(save_dir, "baseline_kinematic_distributions.png")
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    if show_img:
        plt.show()
    else:
        print(f"Saved {save_path}")
    plt.close(fig)


def plot_efficiency_vs_score(output_data, output_mc, trig_vars_2d, triggers,
                             save_dir=None, year='2022', show_img=False,
                             pt_min_plateau=PT_MIN_AN, msd_min_plateau=MSD_MIN_AN,
                             score_bins=None):
    """
    AN Fig 9 equivalent: per-trigger efficiency in the kinematic plateau as a
    function of leading-FatJet (txbb+txcc) tagger score.

    Two-panel figure (data left, MC right). Each curve is one trigger; an extra
    "Soup" curve shows the OR. Plateau cuts pT >= pt_min_plateau and mSD >=
    msd_min_plateau are applied so the variation reflects tagger-correlation
    rather than kinematic turn-on.
    """
    if save_dir is None:
        save_dir = "./figures_sf_2d"
    os.makedirs(save_dir, exist_ok=True)

    var_name = 'pt_vs_msd'
    if var_name not in trig_vars_2d:
        return
    if score_bins is None:
        score_bins = np.linspace(0, 1, 21)

    trigger_names  = list(triggers)
    trigger_labels = [t.replace('HLT_', '') for t in trigger_names]
    all_keys   = trigger_names + ['OR']
    all_labels = trigger_labels + ['Soup']

    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = [prop_cycle[i % len(prop_cycle)] for i in range(len(all_keys))]

    def _curves(output, key):
        if key == 'OR':
            ps = np.asarray(output['Numerator'][var_name].get('pass_score', []))
            px = np.asarray(output['Numerator'][var_name].get('pass_x', []))
            py = np.asarray(output['Numerator'][var_name].get('pass_y', []))
        else:
            if 'PerTrigger' not in output or key not in output['PerTrigger']:
                return None, None, None
            ps = np.asarray(output['PerTrigger'][key][var_name].get('pass_score', []))
            px = np.asarray(output['PerTrigger'][key][var_name].get('pass_x', []))
            py = np.asarray(output['PerTrigger'][key][var_name].get('pass_y', []))
        bs = np.asarray(output['Baseline'][var_name].get('score', []))
        bx = np.asarray(output['Baseline'][var_name].get('x', []))
        by = np.asarray(output['Baseline'][var_name].get('y', []))

        if bs.size == 0:
            return None, None, None

        plateau_b = (by >= pt_min_plateau) & (bx >= msd_min_plateau) & np.isfinite(bs)
        plateau_p = (py >= pt_min_plateau) & (px >= msd_min_plateau) & np.isfinite(ps)

        n_b, _ = np.histogram(bs[plateau_b], bins=score_bins)
        n_p, _ = np.histogram(ps[plateau_p], bins=score_bins)
        eff, err = binomial_eff_err(n_p, n_b)
        centers = 0.5 * (score_bins[:-1] + score_bins[1:])
        return centers, eff, err

    def _draw_panel(ax, output, is_data):
        for i, (key, label) in enumerate(zip(all_keys, all_labels)):
            centers, eff, err = _curves(output, key)
            if centers is None:
                continue
            mk  = '+' if label == 'Soup' else 'o'
            mks = 8   if label == 'Soup' else 4
            ax.errorbar(
                centers, eff, yerr=err, label=label, color=colors[i],
                linestyle='-', marker=mk, markersize=mks, linewidth=2,
                elinewidth=1, capsize=2,
            )
        ax.axhline(1.0, color='gray', linestyle='--', linewidth=1)
        ax.set_ylim(0, 1.4)
        ax.set_xlim(score_bins[0], score_bins[-1])
        ax.legend(fontsize=8, loc='lower right', ncol=1)
        hep.cms.label(ax=ax, data=is_data, year=year, com="13.6", fontsize=10)

    fig, (ax_d, ax_mc) = plt.subplots(1, 2, figsize=(18, 7))
    _draw_panel(ax_d,  output_data, is_data=True)
    _draw_panel(ax_mc, output_mc,   is_data=False)
    plateau_caption = (f"jet $p_T \\geq {pt_min_plateau:.0f}$, "
                       f"$m_{{SD}} \\geq {msd_min_plateau:.0f}$")
    for ax, title in [(ax_d,  f"Triggers ($\\mu$ ref, {plateau_caption})"),
                      (ax_mc, f"Triggers (MC, {plateau_caption})")]:
        ax.set_xlabel("Leading jet $X_{bb} + X_{cc}$ score", fontsize=13)
        ax.set_ylabel("Trigger Efficiency", fontsize=13)
        ax.set_title(title, fontsize=11)
    plt.tight_layout()
    save_path = os.path.join(save_dir, "efficiency_vs_score.png")
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    if show_img:
        plt.show()
    else:
        print(f"Saved {save_path}")
    plt.close(fig)


def _era_from_key(file_key):
    """Map a MuonData file key to an AN-style era label.

    Examples:
      'Muon0_Run2023C-v3' -> '2023C'
      'Muon_Run2022E'     -> '2022E'
      'Muon_Run2022C'     -> '2022C'
    Falls back to the original key if no Run<YYYY><ERA> pattern is found.
    """
    import re
    m = re.search(r'Run(\d{4})([A-Z])', file_key)
    if m:
        return f"{m.group(1)}{m.group(2)}"
    return file_key


def merge_data_outputs(per_era_outputs):
    """
    Concatenate per-era data accumulator dicts into a single combined output,
    so the downstream 2D / 1D plotting functions can run unchanged.
    """
    if len(per_era_outputs) == 0:
        return None
    if len(per_era_outputs) == 1:
        return next(iter(per_era_outputs.values()))

    keys = list(per_era_outputs.keys())
    template = per_era_outputs[keys[0]]
    merged = {
        'Baseline':  {vn: {sub: [] for sub in template['Baseline'][vn]}  for vn in template['Baseline']},
        'Numerator': {vn: {sub: [] for sub in template['Numerator'][vn]} for vn in template['Numerator']},
        'PerTrigger': {
            trig: {vn: {sub: [] for sub in template['PerTrigger'][trig][vn]}
                   for vn in template['PerTrigger'][trig]}
            for trig in template['PerTrigger']
        },
    }
    for k in keys:
        out = per_era_outputs[k]
        for vn in out['Baseline']:
            for sub in out['Baseline'][vn]:
                merged['Baseline'][vn][sub].append(np.asarray(out['Baseline'][vn][sub]))
        for vn in out['Numerator']:
            for sub in out['Numerator'][vn]:
                merged['Numerator'][vn][sub].append(np.asarray(out['Numerator'][vn][sub]))
        for trig in out['PerTrigger']:
            for vn in out['PerTrigger'][trig]:
                for sub in out['PerTrigger'][trig][vn]:
                    merged['PerTrigger'][trig][vn][sub].append(
                        np.asarray(out['PerTrigger'][trig][vn][sub])
                    )
    # concatenate
    for vn in merged['Baseline']:
        for sub in merged['Baseline'][vn]:
            merged['Baseline'][vn][sub] = np.concatenate(merged['Baseline'][vn][sub]) \
                if merged['Baseline'][vn][sub] else np.array([])
    for vn in merged['Numerator']:
        for sub in merged['Numerator'][vn]:
            merged['Numerator'][vn][sub] = np.concatenate(merged['Numerator'][vn][sub]) \
                if merged['Numerator'][vn][sub] else np.array([])
    for trig in merged['PerTrigger']:
        for vn in merged['PerTrigger'][trig]:
            for sub in merged['PerTrigger'][trig][vn]:
                merged['PerTrigger'][trig][vn][sub] = np.concatenate(merged['PerTrigger'][trig][vn][sub]) \
                    if merged['PerTrigger'][trig][vn][sub] else np.array([])
    return merged


def generate_per_era_efficiency_csv(per_era_data, output_mc, trig_vars_2d, triggers,
                                     save_dir=None, year='2022',
                                     pt_min=450, pt_max=1000, msd_min=0, msd_max=300):
    """
    Write an AN-style table (Tables 2-4 layout): one row per trigger,
    one column per data era plus a final 'MC' column. Cells are efficiencies
    in the AN kinematic window. Companion '_Err' table is written alongside.

    per_era_data : {era_label: data accumulator dict}, e.g. {'2023C': {...}}
    """
    if save_dir is None:
        save_dir = "./output"
    os.makedirs(save_dir, exist_ok=True)

    import csv
    var_name = 'pt_vs_msd'
    if var_name not in trig_vars_2d:
        return

    eras = sorted(per_era_data.keys())

    def _eff_err_in_window(out, trigger):
        """Return (eff, err, k_pass, n_base) integrated over the AN window."""
        bx = np.asarray(out['Baseline'][var_name]['x'])
        by = np.asarray(out['Baseline'][var_name]['y'])
        base_cut = (by >= pt_min) & (by < pt_max) & (bx >= msd_min) & (bx < msd_max)
        n_base = int(np.sum(base_cut))

        if trigger == 'OR':
            px = np.asarray(out['Numerator'][var_name]['pass_x'])
            py = np.asarray(out['Numerator'][var_name]['pass_y'])
        else:
            if 'PerTrigger' not in out or trigger not in out['PerTrigger']:
                return 0.0, 0.0, 0, n_base
            px = np.asarray(out['PerTrigger'][trigger][var_name]['pass_x'])
            py = np.asarray(out['PerTrigger'][trigger][var_name]['pass_y'])
        pass_cut = (py >= pt_min) & (py < pt_max) & (px >= msd_min) & (px < msd_max)
        k_pass = int(np.sum(pass_cut))
        eff_arr, err_arr = binomial_eff_err(np.array([k_pass]), np.array([n_base]))
        eff = float(eff_arr[0]) if np.isfinite(eff_arr[0]) else 0.0
        err = float(err_arr[0]) if np.isfinite(err_arr[0]) else 0.0
        return eff, err, k_pass, n_base

    rows_eff = [['Trigger'] + eras + ['MC']]
    rows_err = [['Trigger'] + eras + ['MC']]

    trigger_keys = list(triggers) + ['OR']
    for trig in trigger_keys:
        label = 'Combined (OR)' if trig == 'OR' else trig.replace('HLT_', '')
        eff_row = [label]
        err_row = [label]
        for era in eras:
            eff, err, _, _ = _eff_err_in_window(per_era_data[era], trig)
            eff_row.append(f"{eff:.4f}")
            err_row.append(f"{err:.4f}")
        eff_mc, err_mc, _, _ = _eff_err_in_window(output_mc, trig)
        eff_row.append(f"{eff_mc:.4f}")
        err_row.append(f"{err_mc:.4f}")
        rows_eff.append(eff_row)
        rows_err.append(err_row)

    eff_path = os.path.join(save_dir, f"per_era_efficiency_{year}.csv")
    err_path = os.path.join(save_dir, f"per_era_efficiency_{year}_err.csv")
    with open(eff_path, 'w', newline='') as f:
        csv.writer(f).writerows(rows_eff)
    with open(err_path, 'w', newline='') as f:
        csv.writer(f).writerows(rows_err)
    print(f"Saved per-era efficiency table to: {eff_path}")
    print(f"Saved per-era efficiency error  to: {err_path}")


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
    csv_data.append(['Trigger', 'Data_Efficiency', 'Data_Eff_Err', 'Data_Pass', 'Data_Total',
                     'MC_Efficiency', 'MC_Eff_Err', 'MC_Pass', 'MC_Total',
                     'Scale_Factor', 'Scale_Factor_Err'])

    def _row(label, k_d, n_d, k_m, n_m):
        eff_d_arr, err_d_arr = binomial_eff_err(np.array([k_d]), np.array([n_d]))
        eff_m_arr, err_m_arr = binomial_eff_err(np.array([k_m]), np.array([n_m]))
        eff_d = float(eff_d_arr[0]) if np.isfinite(eff_d_arr[0]) else 0.0
        eff_m = float(eff_m_arr[0]) if np.isfinite(eff_m_arr[0]) else 0.0
        err_d = float(err_d_arr[0]) if np.isfinite(err_d_arr[0]) else 0.0
        err_m = float(err_m_arr[0]) if np.isfinite(err_m_arr[0]) else 0.0
        _, sf_err_arr = sf_err_from_counts(np.array([k_d]), np.array([n_d]),
                                           np.array([k_m]), np.array([n_m]))
        sf = eff_d / eff_m if eff_m > 0 else 1.0
        sf_err = float(sf_err_arr[0]) if np.isfinite(sf_err_arr[0]) else 0.0
        print(f"{label:60s} | Data: {eff_d:.2%} ± {err_d:.2%} ({k_d}/{n_d}) | "
              f"MC: {eff_m:.2%} ± {err_m:.2%} ({k_m}/{n_m}) | SF: {sf:.4f} ± {sf_err:.4f}")
        return [label,
                f"{eff_d:.4f}", f"{err_d:.4f}", f"{k_d}", f"{n_d}",
                f"{eff_m:.4f}", f"{err_m:.4f}", f"{k_m}", f"{n_m}",
                f"{sf:.4f}",   f"{sf_err:.4f}"]

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

        n_pass_data = int(np.sum(pass_cut_data))
        n_pass_mc   = int(np.sum(pass_cut_mc))

        csv_data.append(_row(trigger.replace('HLT_', ''),
                             n_pass_data, n_baseline_data,
                             n_pass_mc,   n_baseline_mc))

    # Calculate combined (OR) efficiency
    pass_x_data = np.array(output_data['Numerator'][var_name]['pass_x'])
    pass_y_data = np.array(output_data['Numerator'][var_name]['pass_y'])
    pass_x_mc = np.array(output_mc['Numerator'][var_name]['pass_x'])
    pass_y_mc = np.array(output_mc['Numerator'][var_name]['pass_y'])

    pass_cut_data = (pass_y_data >= pt_min) & (pass_y_data < pt_max) & \
                    (pass_x_data >= msd_min) & (pass_x_data < msd_max)
    pass_cut_mc = (pass_y_mc >= pt_min) & (pass_y_mc < pt_max) & \
                  (pass_x_mc >= msd_min) & (pass_x_mc < msd_max)

    n_pass_data = int(np.sum(pass_cut_data))
    n_pass_mc   = int(np.sum(pass_cut_mc))

    csv_data.append(_row('Combined (OR)',
                         n_pass_data, n_baseline_data,
                         n_pass_mc,   n_baseline_mc))

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
    parser.add_argument('--workers', type=int, default=4,
                        help="Number of FuturesExecutor workers in regular (non-quick) mode "
                             "(default: 4). Bumped from 2 after 48h timeouts on 2022/2022EE/2023.")
    parser.add_argument('--no-recoverable', action='store_true',
                        help="Debug-only: disable FuturesExecutor recoverable=True so per-chunk "
                             "exceptions propagate. Use to diagnose 'No chunks returned results' "
                             "failures (recoverable hides the underlying cause).")
    parser.add_argument('--txbb-region',
                        choices=['fail', 'pass_bb', 'pass_cc', 'pass', 'inclusive', 'all'],
                        default='all',
                        help="Candidate-jet region (AN Sec 5.1): 'fail' (Xbb+Xcc<0.82), "
                             "'pass_bb' (>=0.82 AND Xbb>Xcc), 'pass_cc' (>=0.82 AND Xcc>Xbb), "
                             "'pass' (>=0.82, back-compat combined), 'inclusive' (no region cut), "
                             "or 'all' which runs the three AN regions {fail, pass_bb, pass_cc}.")
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
    parser.add_argument('--mc-sample', choices=['ttbar', 'qcd_muenriched', 'both'], default='both',
                        help="MC sample for the SF denominator. 'ttbar' uses TTtoLNu2Q (may bias "
                             "SF low due to top-like jets); 'qcd_muenriched' uses the full "
                             "QCD_MuEnrichedPt5 PT-binned set (closer match to single-muon data "
                             "jet mix); 'both' (default) runs both in sequence and writes the two "
                             "result trees side-by-side under <mc_sample>/ subdirectories.")
    parser.add_argument('--outdir', default='.',
                        help="Base output directory. Results are saved under <outdir>/figures_sf_2d and <outdir>/output.")
    parser.add_argument('--mc-keys-filter', default=None,
                        help="Optional regex; when --mc-sample resolves to a list of dataset keys "
                             "(e.g. the 12 QCD_MuEnriched pT bins), keep only keys matching this "
                             "regex. Used to split a large MC sample across multiple condor jobs by "
                             "pT bin. Example: --mc-keys-filter '^QCD_MuEnriched_PT_300to470$'.")
    parser.add_argument('--mc-output-suffix', default='',
                        help="Suffix to append to coffea/CSV/figure output filenames so multiple "
                             "--mc-keys-filter sub-jobs of the same year/region/mc don't clobber each "
                             "other. Example: '_PT_300to470'. Empty string (default) leaves names "
                             "unchanged for backwards compatibility.")
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
    # Muon dataset keys per period (data side, unchanged by --mc-sample)
    # For 2023, only v3 and v4 re-recos contain the
    # HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06 branch (see scan_2023 study,
    # 2026-05-20). v1/v2 files lack it entirely (0/121 in Muon0, 0/106 in Muon1)
    # so including them silently halves the OR-eff in data and produces a
    # spurious ~30-point data/MC gap. Restrict to v3+ here.
    muon_keys_by_period = {
        '2022':     ["Muon_Run2022C", "Muon_Run2022D"],
        '2022EE':   ["Muon_Run2022E", "Muon_Run2022F", "Muon_Run2022G"],
        '2023':     [
            "Muon0_Run2023C-v3", "Muon0_Run2023C-v4",
            "Muon1_Run2023C-v3", "Muon1_Run2023C-v4",
        ],
        '2023BPix': [
            "Muon0_Run2023D-v1", "Muon0_Run2023D-v2",
            "Muon1_Run2023D-v1", "Muon1_Run2023D-v2",
        ],
        '2024': [
            # NB: 2024 has separate Muon0 and Muon1 PDs (like 2023). The coworker's
            # original JSON mashed both into the Muon0_* keys; we split them out
            # on 2026-06-03 to match the 2022/2023 pipeline. No Muon1_Run2024H
            # because that era's NanoAOD only carries Muon0 PD.
            "Muon0_Run2024C", "Muon0_Run2024D", "Muon0_Run2024E",
            "Muon0_Run2024F", "Muon0_Run2024G", "Muon0_Run2024H", "Muon0_Run2024I",
            "Muon1_Run2024C", "Muon1_Run2024D", "Muon1_Run2024E",
            "Muon1_Run2024F", "Muon1_Run2024G", "Muon1_Run2024I",
        ],
    }

    def _resolve_mc_source(mc_sample):
        """Map --mc-sample value -> (dataset_type for infile JSON, file_key string-or-list).

        When --mc-keys-filter is set and the resolved key set is a list (e.g. the
        12 QCD_MuEnriched pT bins), apply the regex filter to it. Used by
        submit_sf_2d_qcd_ptbin.sh to fan one slow QCD job into 12 sub-jobs.
        """
        if mc_sample == 'ttbar':
            return 'ttbar', "TTtoLNu2Q"
        if mc_sample == 'qcd_muenriched':
            qcd_json_path = os.path.join(f'infiles/{year}', f"{year}_QCD_MuEnriched.json")
            if not os.path.exists(qcd_json_path):
                raise FileNotFoundError(
                    f"--mc-sample qcd_muenriched needs {qcd_json_path}; "
                    "run make_qcd_muenriched_infiles.py after voms-proxy-init."
                )
            with open(qcd_json_path) as _f:
                keys = list(json.load(_f).keys())
            if args.mc_keys_filter:
                import re as _re
                pat = _re.compile(args.mc_keys_filter)
                kept = [k for k in keys if pat.search(k)]
                if not kept:
                    raise ValueError(
                        f"--mc-keys-filter {args.mc_keys_filter!r} matched none of {keys}"
                    )
                print(f"--mc-keys-filter kept {len(kept)} / {len(keys)} keys: {kept}")
                keys = kept
            return 'QCD_MuEnriched', keys
        raise ValueError(f"Unknown mc-sample: {mc_sample}")

    if args.mc_sample == 'both':
        mc_samples = ['ttbar', 'qcd_muenriched']
    else:
        mc_samples = [args.mc_sample]
    print(f"MC denominators to process: {mc_samples}")

    # Trigger SFs are measured inclusively (not separated by production mode)
    prod_modes = ['Inclusive']

    # AN Section 5.1 regions
    if args.txbb_region == 'all':
        txbb_regions = ['fail', 'pass_bb', 'pass_cc']
    else:
        txbb_regions = [args.txbb_region]

    for mc_sample in mc_samples:
      mc_dataset_type, mc_file_key = _resolve_mc_source(mc_sample)
      prod_modes_map = {
          mc_dataset_type: mc_file_key,
          'MuonData':      muon_keys_by_period[year],
      }
      print(f"\n{'#'*60}")
      print(f"# MC denominator: {mc_sample} ({mc_dataset_type})")
      print(f"{'#'*60}")

      for prod_mode in prod_modes:
        for txbb_region in txbb_regions:
            print(f"\n{'='*60}")
            print(f"Measuring trigger SF: {prod_mode}, txbb={txbb_region} (WP={TXBB_WP})")
            print(f"{'='*60}\n")

        # Load fileset once per period (shared across txbb regions).
        # MuonData with multiple keys is split by AN-style era (e.g. 2023C) so
        # the per-era CSV (Tables 2-4 layout) can be produced downstream.
        fileset = {}
        # Map fileset key -> era label, for per-era CSV (only populated for MuonData).
        era_lookup = {}

        for dataset_type, file_key in prod_modes_map.items():
            json_path = os.path.join(f'infiles/{year}', f"{year}_{dataset_type}.json")
            print(f"Reading: {json_path} (key: {file_key})")

            if not os.path.exists(json_path):
                print(f"  Warning: {json_path} not found. Skipping {dataset_type}.")
                continue

            with open(json_path, 'r') as file:
                data = json.load(file)

            # Group source-keys per era. For MuonData each AN-era (e.g. 2023C) is
            # one fileset; for everything else (ttbar, QCD_MuEnriched pT-bins) all
            # sources are merged into a single fileset entry. We keep the
            # *sub-lists* so test-mode truncation can round-robin across them
            # and avoid biasing toward the first source key (e.g. lowest pT bin).
            split_by_era = (dataset_type == 'MuonData')
            era_to_sublists = {}  # era_label -> list of file lists (one per source key)
            if isinstance(file_key, list):
                for k in file_key:
                    s = data.get(k, [])
                    era = _era_from_key(k) if split_by_era else file_key[0]
                    era_to_sublists.setdefault(era, []).append(s)
                    print(f"  + {k}: {len(s)} files (era={era})")
            else:
                era = _era_from_key(file_key) if split_by_era else dataset_type
                era_to_sublists = {era: [list(data.get(file_key, []))]}

            for era, sublists in era_to_sublists.items():
                if quick_test_mode:
                    samples = round_robin_sample_lists(sublists, args.quick_files_per_dataset)
                    print(f"QUICK TEST MODE [{era}]: round-robin {len(samples)} files "
                          f"across {len(sublists)} source key(s)")
                elif test_mode:
                    samples = round_robin_sample_lists(sublists, 10)
                    print(f"TEST MODE [{era}]: round-robin {len(samples)} files "
                          f"across {len(sublists)} source key(s)")
                else:
                    samples = [f for sub in sublists for f in sub]

                if split_by_era:
                    dataset_key = f"{dataset_type}_{year}_{era}"
                    era_lookup[dataset_key] = era
                else:
                    dataset_key = f"{dataset_type}_{year}"
                fileset[dataset_key] = samples
                print(f"Loaded {len(samples)} files for dataset '{dataset_key}'")

        if len(fileset) == 0:
            print(f"Error: No data loaded for production mode {prod_mode}")
            continue

        print("Constructed fileset:")
        print(json.dumps(fileset, indent=2))

        # ---- sanity check: every analysis trigger has a branch in the first file ----
        # Catches regressions like "added a new trigger name but the v1/v2 NanoAODs
        # don't actually carry that branch" (the bug that caused the original 2023
        # 30-pt SF gap). Soft warning: we don't kill the job, since transient xrootd
        # failures at startup shouldn't be fatal, but every missing branch is logged
        # loudly so a grep -E "TRIGGER-MISSING|TRIGGER-OK" in the .out catches it.
        _expected_trigs = list(trigger_dict_periods.get(year, [])) + list(muon_triggers.get(year, []))
        for _ds_key, _files in fileset.items():
            if not _files:
                continue
            try:
                with uproot.open(_files[0]) as _f:
                    _branches = set(_f["Events"].keys())
            except Exception as _exc:
                print(f"  TRIGGER-CHECK-SKIP {_ds_key}: could not open {_files[0]}: {_exc!r}")
                continue
            for _t in _expected_trigs:
                _tfield = _t[4:] if _t.startswith("HLT_") else _t
                _hlt_name = f"HLT_{_tfield}"
                if _hlt_name in _branches:
                    print(f"  TRIGGER-OK      {_ds_key}: {_hlt_name}")
                else:
                    print(f"  TRIGGER-MISSING {_ds_key}: {_hlt_name} (downstream OR will silently skip this leg!)")

        workers = args.quick_workers if quick_test_mode else args.workers
        runner_maxchunks = args.quick_maxchunks if quick_test_mode else None
        # xrootdtimeout=300 (vs coffea default 60s) — gives slow xrootd reads time
        # to succeed; preprocessing has been seen to silently drop ALL files when
        # transient timeouts hit, leading to empty chunks and a confusing
        # TypeError in FuturesExecutor.__call__.
        # recoverable=True swallows per-chunk exceptions and returns the partial
        # accumulator (or None if everything failed). Use --no-recoverable for a
        # debug run that lets the underlying chunk exception propagate so we can
        # actually see why chunks are dying.
        iterative_run = processor.Runner(
            executor=processor.FuturesExecutor(
                compression=None, workers=workers,
                recoverable=(not args.no_recoverable),
            ),
            schema=NanoAODSchema,
            skipbadfiles=True,
            savemetrics=True,
            chunksize=args.chunksize,
            maxchunks=runner_maxchunks,
            xrootdtimeout=300,
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

                # Guard against the coffea executor returning None when *all*
                # files fail metadata fetch (e.g. xrootd timeouts). Without this
                # the whole job dies on the first bad region and later regions
                # never run.
                try:
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
                except Exception as exc:
                    print(f"  ERROR: iterative_run raised for {dataset_key} [{txbb_region}]: {exc!r}")
                    print(f"  Continuing to next dataset/region.")
                    continue
                if out is None or not isinstance(out, tuple) or out[0] is None:
                    print(f"  WARNING: iterative_run returned no output for {dataset_key} [{txbb_region}] "
                          "(likely all files failed metadata fetch — check xrootd). Skipping.")
                    continue
                outputs[dataset_key] = out[0]

            data_keys = [k for k in outputs.keys() if 'Muon' in k]
            mc_key    = next((k for k in outputs.keys() if k.startswith(mc_dataset_type)), None)

            if not data_keys or mc_key is None:
                print("Error: Missing data or MC output. Skipping plots.")
                continue

            output_mc = outputs[mc_key]
            # Per-era data outputs (for AN-style table); fall back to a single
            # implicit "data" key when MuonData wasn't split by era.
            per_era_data = {era_lookup.get(k, 'data'): outputs[k] for k in data_keys}
            # Merge across eras for the 2D / 1D plots, which are integrated over the period.
            output_data = merge_data_outputs(per_era_data) if len(per_era_data) > 1 else next(iter(per_era_data.values()))

            # Output paths include txbb_region AND mc_sample so the ttbar and
            # qcd_muenriched results live side-by-side without clobbering.
            # --mc-output-suffix lets pT-bin sub-jobs land alongside each other.
            mc_suffix = getattr(args, 'mc_output_suffix', '') or ''
            FIGURES_DIR = os.path.join(figures_base_dir, year, prod_mode, txbb_region, mc_sample + mc_suffix)
            os.makedirs(FIGURES_DIR, exist_ok=True)

            output_path_data = os.path.join(output_base_dir, f"sf_2d_{year}_{prod_mode}_{txbb_region}_{mc_sample}{mc_suffix}_data.coffea")
            output_path_mc   = os.path.join(output_base_dir, f"sf_2d_{year}_{prod_mode}_{txbb_region}_{mc_sample}{mc_suffix}_mc.coffea")
            # Per-era dict so regenerate_csv_from_coffea.py can rebuild the
            # per-era CSVs from disk. The merged output_data alone loses the
            # era split; without this file the per-era CSVs would need to be
            # generated in-process (as they are here in generate_per_era_...).
            output_path_per_era = os.path.join(output_base_dir, f"sf_2d_{year}_{prod_mode}_{txbb_region}_{mc_sample}{mc_suffix}_per_era.coffea")
            os.makedirs(os.path.dirname(output_path_data), exist_ok=True)

            util.save(output_data, output_path_data)
            util.save(per_era_data, output_path_per_era)
            util.save(output_mc,   output_path_mc)
            print(f"Saved data output to {output_path_data}")
            print(f"Saved per-era coffea to {output_path_per_era}")
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

            # AN Fig 9 equivalent: per-trigger eff vs leading-jet (txbb+txcc) score
            plot_efficiency_vs_score(
                output_data, output_mc, trig_vars_2d,
                triggers=trigger_dict_periods[year],
                save_dir=FIGURES_DIR,
                year=year,
                show_img=False,
                pt_min_plateau=pt_min_csv,
                msd_min_plateau=msd_min_csv,
            )

            # 1D distributions of leading-jet pT, mSD, TXbb+TXcc after baseline cuts
            # and the chosen AN region split (data + MC overlay with Data/MC ratio).
            plot_baseline_kinematic_distributions(
                output_data, output_mc,
                save_dir=FIGURES_DIR, year=year, region=txbb_region, show_img=False,
            )
            print(f"\nAll plots saved to {FIGURES_DIR}")

            csv_output_dir = os.path.join(output_base_dir, year, prod_mode, txbb_region, mc_sample + mc_suffix)
            os.makedirs(csv_output_dir, exist_ok=True)

            generate_per_trigger_efficiency_csv(
                output_data, output_mc, trig_vars_2d,
                triggers=trigger_dict_periods[year],
                save_dir=csv_output_dir,
                year=f"{year}_{prod_mode}_{txbb_region}_{mc_sample}{mc_suffix}",
                pt_min=pt_min_csv,
                pt_max=pt_max_csv,
                msd_min=msd_min_csv,
                msd_max=msd_max_csv
            )

            # AN Tables 2-4 equivalent: per-era columns + MC column.
            # Only meaningful when MuonData was split by era (i.e., at least one era key).
            if len(per_era_data) >= 1:
                generate_per_era_efficiency_csv(
                    per_era_data, output_mc, trig_vars_2d,
                    triggers=trigger_dict_periods[year],
                    save_dir=csv_output_dir,
                    year=f"{year}_{prod_mode}_{txbb_region}_{mc_sample}{mc_suffix}",
                    pt_min=pt_min_csv,
                    pt_max=pt_max_csv,
                    msd_min=msd_min_csv,
                    msd_max=msd_max_csv,
                )

    print(f"\n{'='*60}")
    print(f"Completed processing for all production modes")
    print(f"{'='*60}")
