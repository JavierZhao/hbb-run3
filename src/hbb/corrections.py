"""
Collection of utilities for corrections and systematics in processors.

Most corrections retrieved from the cms-nanoAOD repo:
See https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/
"""

from __future__ import annotations

import pathlib
from pathlib import Path

import contextlib

import awkward as ak
import dask_awkward as dak
import numpy as np
import correctionlib
import correctionlib.schemav2
import pickle
from coffea.analysis_tools import Weights
from coffea.nanoevents.methods import vector
from coffea.nanoevents.methods.nanoaod import JetArray, MuonArray
from coffea.jetmet_tools import CorrectedJetsFactory, CorrectedMETFactory, JECStack
from coffea.lookup_tools import extractor

from hbb.MuonScaRe import pt_resol, pt_scale, pt_resol_var, pt_scale_var 
from hbb.jerc_eras import jec_eras,jer_eras, jec_mc, jer_mc, jec_data, fatjet_jerc_keys, jet_jerc_keys
from hbb.taggers import b_taggers
from hbb.EWHiggs_corrections import theory_xs, xs_ewkcorr, ewh_ptbin

ak.behavior.update(vector.behavior)
package_path = str(pathlib.Path(__file__).parent.parent.resolve())

# Important Run3 start of Run
FirstRun_2022C = 355794
FirstRun_2022D = 357487
LastRun_2022D = 359021
FirstRun_2022E = 359022
LastRun_2022F = 362180

"""
CorrectionLib files are available from: /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration - synced daily
"""
pog_correction_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/"
pog_jsons = {
    "muon": ["MUO", "muon_Z.json.gz"],
    "muon_pt" : ["MUO", "muon_scalesmearing.json.gz"],
    "electron": ["EGM", "electron.json.gz"],
    "photon": ["EGM", "photon.json.gz"],
    "electronSS": ["EGM", "electronSS_EtDependent.json.gz"],
    "photonSS": ["EGM", "photonSS_EtDependent.json.gz"],
    "pileup": ["LUM", "puWeights.json.gz"],
    "pileup2024": ["LUM", "puWeights_CDEFGHI.json.gz"], # file excludes 2024B 
    "fatjet_jec": ["JME", "fatJet_jerc.json.gz"],
    "jet_jec": ["JME", "jet_jerc.json.gz"],
    "jetveto": ["JME", "jetvetomaps.json.gz"],
    "btagging": ["BTV", "btagging.json.gz"],
    "jetid" : ["JME", "jetid.json.gz"],
}

years = {
    "2022": "Run3-22CDSep23-Summer22-NanoAODv12",
    "2022EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
    "2023": "Run3-23CSep23-Summer23-NanoAODv12",
    "2023BPix": "Run3-23DSep23-Summer23BPix-NanoAODv12",
    "2024": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
}


def ak_clip(arr: ak.Array, min_value: float, max_value: float):
    """
    Clip the values of an awkward array using where
    """
    return ak.where(arr < min_value, min_value, ak.where(arr > max_value, max_value, arr))


def get_pog_json(obj: str, year: str) -> str:
    try:
        pog_json = pog_jsons[obj]
    except:
        print(f"No json for {obj}")

    year = years[year]

    return f"{pog_correction_path}/{pog_json[0]}/{year}/latest/{pog_json[1]}"


def build_lumimask(filename):
    from coffea.lumi_tools import LumiMask

    path = Path(f"{package_path}/hbb/data/{filename}")
    return LumiMask(path)


lumiMasks = {
    "2022": build_lumimask("Cert_Collisions2022_355100_362760_Golden.json"),
    "2022EE": build_lumimask("Cert_Collisions2022_355100_362760_Golden.json"),
    "2023": build_lumimask("Cert_Collisions2023_366442_370790_Golden.json"),
    "2023BPix": build_lumimask("Cert_Collisions2023_366442_370790_Golden.json"),
    "2024": build_lumimask("Cert_Collisions2024_378981_386951_Golden.json"),
}


def add_pileup_weight(weights: Weights, year: str, nPU):
    # clip nPU from 0 to 150
    # nPU = ak_clip(nPU, 0, 99)

    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3
    values = {}

    if not year == "2024":
        cset = correctionlib.CorrectionSet.from_file(get_pog_json("pileup", year))
    else:
        pog_json_file = f"{package_path}/hbb/data/puWeights_2024.json"
        cset = correctionlib.CorrectionSet.from_file(get_pog_json("pileup2024", year))

    corr = {
        "2018": "Collisions18_UltraLegacy_goldenJSON",
        "2022": "Collisions2022_355100_357900_eraBCD_GoldenJson",
        "2022EE": "Collisions2022_359022_362760_eraEFG_GoldenJson",
        "2023": "Collisions2023_366403_369802_eraBC_GoldenJson",
        "2023BPix": "Collisions2023_369803_370790_eraD_GoldenJson",
        "2024": "Collisions24_CDEFGHI_goldenJSON",
    }[year]
    # evaluate and clip up to 4 to avoid large weights
    values["nominal"] = ak_clip(cset[corr].evaluate(nPU, "nominal"), 0, 4)
    values["up"] = ak_clip(cset[corr].evaluate(nPU, "up"), 0, 4)
    values["down"] = ak_clip(cset[corr].evaluate(nPU, "down"), 0, 4)

    weights.add("pileup", values["nominal"], values["up"], values["down"])

def add_pdf_weight(gen_weights, pdf_weights, output):
    """
    Apply pdf weight variation for standard Hessian set
    """

    out_pdf = {}
    for i in range(103):
        out_pdf[f"weight_pdf_{i}"] = pdf_weights[:, i]
        output["sumw_pdf"][f"sumweight_pdf_{i}"] = ak.sum(pdf_weights[:, i] * gen_weights)

    return out_pdf

def add_ps_weight(weights: Weights, ps_weights):
    """
    Parton Shower Weights (FSR and ISR)
    """
    nom = ak.ones_like(weights.weight())

    up_isr = ak.ones_like(nom)
    down_isr = ak.ones_like(nom)
    up_fsr = ak.ones_like(nom)
    down_fsr = ak.ones_like(nom)

    if ak.num(ps_weights[0], axis=0) == 4:
        up_isr = ps_weights[:, 0]  # ISR=2, FSR=1
        down_isr = ps_weights[:, 2]  # ISR=0.5, FSR=1

        up_fsr = ps_weights[:, 1]  # ISR=1, FSR=2
        down_fsr = ps_weights[:, 3]  # ISR=1, FSR=0.5

    elif ak.num(ps_weights[0], axis=0) > 1:
        print("PS weight vector has length ", ak.num(ps_weights[0]))

    weights.add("ISRPartonShower", nom, up_isr, down_isr)
    weights.add("FSRPartonShower", nom, up_fsr, down_fsr)

def add_scalevar(gen_weights, var_weights, output, structure = "7pt"):
    """
    QCD scale variations according to recommendations by the LHCXSWG
    For application to:
         high pt ggf and ttH higgs production mc : muF = muR
         high pt VBF and VH higgs production mc : muF^2 = muR^2
    Recommendation by LHCXSWG cds.cern.ch/record/2669113
    """
    var_map = {
        "3pt" : [0, 4, 8],             # case where muF^2 = muR^2
        "7pt" : [0, 1, 3, 4, 5, 7, 8]  # case where muF = muR 
    }

    out_lhe = {}
    try:
        for var in var_map[structure]:
            out_lhe[f"weight_scalevar_{structure}_{var}"] = var_weights[:, var]
            output["sumw_scalevar"][f"sumweight_scalevar_{structure}_{var}"] = ak.sum(var_weights[:, var] * gen_weights)
            
    except Exception as e:
        print("Scale variation structure unexpected:", e)

    return out_lhe

def get_EWHiggs_corrector(prodmode: str):
    #Create the corrector for the EW Higgs xs corrections based on selected production mode

    # make the bins exclusive
    def make_excl(th):
        excl = -1*np.diff(th)
        return np.append(excl,[th[-1]])

    noewcorr = make_excl(theory_xs[prodmode])
    yesewcorr = make_excl(np.multiply(theory_xs[prodmode], 1 + xs_ewkcorr[prodmode]))

    # xs scale factors
    weights = np.divide(yesewcorr,noewcorr)

    # correction input variable
    hptvar = correctionlib.schemav2.Variable(
        name="hpt",
        type="real",
        description="Generated Higgs boson pT"
    )

    corr = correctionlib.schemav2.Correction(
        name=prodmode,
        version=1,
        description=f"Electroweak correction to {prodmode} Higgs production",
        inputs=[hptvar],
        output=correctionlib.schemav2.Variable(name="out", type="real", description="Multiplicative k-factor"),
        data={
            "nodetype": "binning",
            "input": "hpt",
            "edges": list(ewh_ptbin),
            "content": list(weights),
            "flow": "clamp",
        },
    )

    return corr

def add_EWHiggs_weight(weights: Weights, dataset: str, genpart):
    # Apply EW Higgs xs corrections

    boson = ak.firsts(genpart[
        (genpart.pdgId == 25)
        & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
    ])
    boson_pt = ak.fill_none(boson.pt, 0.)

    if "VBFH" in dataset:
        prodmode = "VBF"
    elif "WplusH" in dataset or "WminusH" in dataset or "ZH" in dataset:
        prodmode = "VH"
    elif "ttH" in dataset:
        prodmode = "ttH"
    else:
        return
    
    corr = get_EWHiggs_corrector(prodmode)
    ewk_nominal = corr.to_evaluator().evaluate(boson_pt)

    weights.add(f"{prodmode}_EW", ewk_nominal)

# Jet Veto Maps
# the JERC group recommends ALL analyses use these maps, as the JECs are derived excluding these zones.
# apply to both Data and MC
# https://cms-talk.web.cern.ch/t/jet-veto-maps-for-run3-data/18444?u=anmalara
# https://cms-talk.web.cern.ch/t/jes-for-2022-re-reco-cde-and-prompt-fg/32873
def get_jetveto_event(jets: JetArray, year: str):
    """
    Get event selection that rejects events with jets in the veto map
    """

    # correction: Non-zero value for (eta, phi) indicates that the region is vetoed
    cset = correctionlib.CorrectionSet.from_file(get_pog_json("jetveto", year))
    j, nj = ak.flatten(jets), ak.num(jets)

    def get_veto(j, nj, csetstr):
        j_phi = ak_clip(j.phi, -3.1415, 3.1415)
        j_eta = ak_clip(j.eta, -4.7, 4.7)
        veto = cset[csetstr].evaluate("jetvetomap", j_eta, j_phi)
        return ak.unflatten(veto, nj)

    corr_str = {
        "2022": "Summer22_23Sep2023_RunCD_V1",
        "2022EE": "Summer22EE_23Sep2023_RunEFG_V1",
        "2023": "Summer23Prompt23_RunC_V1",
        "2023BPix": "Summer23BPixPrompt23_RunD_V1",
        "2024": "Summer24Prompt24_RunBCDEFGHI_V1",
    }[year]

    jet_veto = get_veto(j, nj, corr_str) > 0

    event_sel = ~(ak.any((jets.pt > 15) & (jets.jetidtightlepveto) & jet_veto, axis=1))
    return event_sel

def correct_jetid(jets, jet_type: str, year: str):
    """
    Apply jetid correction for v14+
    https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#nanoAOD_Flags
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/199ba071f68176a815615651f8a5ae939ef0793e/examples/jetidExample.py
    """
    evaluator = correctionlib.CorrectionSet.from_file(get_pog_json("jetid", year))

    if jet_type == "AK8":
        name_tight = "AK8PUPPI_Tight"
        name_tightlv = "AK8PUPPI_TightLeptonVeto"
    elif jet_type == "AK4":
        name_tight = "AK4PUPPI_Tight"
        name_tightlv = "AK4PUPPI_TightLeptonVeto"

    flat_j = ak.flatten(jets)

    def get_jetid(j, nj, ev_str):
 
        eta    = ak_clip(j.eta, -4.7, 4.7)
        chHEF  = ak_clip(j.chHEF, 0., 1.)
        neHEF  = ak_clip(j.neHEF, 0., 1.)
        chEmEF = ak_clip(j.chEmEF, 0., 1.)
        neEmEF = ak_clip(j.neEmEF, 0., 1.)
        muEF   = ak_clip(j.muEF, 0., 1.)
        chMult = ak_clip(j.chMultiplicity, 0., 100.)
        neMult = ak_clip(j.neMultiplicity, 0., 100.)
        nConst = chMult + neMult

        jetid = evaluator[ev_str].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, nConst)
        return ak.unflatten(jetid, nj)

    jets["jetidtight"] = ak.values_astype(get_jetid(flat_j, ak.num(jets), name_tight), bool)
    jets["jetidtightlepveto"] = ak.values_astype(get_jetid(flat_j, ak.num(jets), name_tightlv), bool)

    return jets

jec_name_map = {
    "JetPt": "pt",
    "JetMass": "mass",
    "JetEta": "eta",
    "JetA": "area",
    "ptGenJet": "pt_gen",
    "ptRaw": "pt_raw",
    "massRaw": "mass_raw",
    "Rho": "event_rho",
    "METpt": "pt",
    "METphi": "phi",
    "JetPhi": "phi",
    "UnClusteredEnergyDeltaX": "MetUnclustEnUpDeltaX",
    "UnClusteredEnergyDeltaY": "MetUnclustEnUpDeltaY",
}

def apply_jerc(jets, jet_type: str, year: str, runkey: str):
    #Create CorrectedJetFactory and apply jercs+variations to JetArray or FatJetArray

    jerc_path =f"{package_path}/hbb/data/jerc"
    jec_path = f"{jerc_path}/{jec_eras[runkey]}"

    if jet_type == "AK8":
        jet_key = fatjet_jerc_keys[year]
    elif jet_type == "AK4":
        jet_key = jet_jerc_keys[year]

    #build filelist
    files = []
    if "mc" in runkey:
        jer_path = f"{jerc_path}/{jer_eras[runkey]}"
        for jec, ps in jec_mc.items():
            files.append(f"{jec_path}/{jec_eras[runkey]}_{jec}_{jet_key}{ps}")
        for jer, ps in jer_mc.items():
            files.append(f"{jer_path}/{jer_eras[runkey]}_{jer}_{jet_key}{ps}")
    else:
        for jec, ps in jec_data.items():
            files.append(f"{jec_path}/{jec_eras[runkey]}_{jec}_{jet_key}{ps}")

    ext = extractor()
    with contextlib.ExitStack() as stack:
        real_files = [stack.enter_context(Path(f)) for f in files]
        ext.add_weight_sets([f"* * {file}" for file in real_files])
        ext.finalize()

    jec_stack = JECStack(ext.make_evaluator())
    jet_factory = CorrectedJetsFactory(jec_name_map, jec_stack)

    corrected_jets = jet_factory.build(jets)
    return corrected_jets

def correct_met(met, jets: JetArray):
    #Create CorrectedMETFactory and recluster met

    dX_up = met.ptUnclusteredUp * np.cos(met.phiUnclusteredUp)
    dY_up = met.ptUnclusteredUp * np.sin(met.phiUnclusteredUp)
    dX_nom = met.pt * np.cos(met.phi)
    dY_nom = met.pt * np.sin(met.phi)
    met["MetUnclustEnUpDeltaX"] = dX_up - dX_nom
    met["MetUnclustEnUpDeltaY"] = dY_up - dY_nom

    met_factory = CorrectedMETFactory(jec_name_map)
    corrected_met = met_factory.build(met, jets)

    return corrected_met

def add_btag_weights(weights: Weights, jets: JetArray, btagger: str, wp: str, year: str, alt_str: str):
    """
    Apply btag event scale factor for AK4 jets queried
    Using BTV fixed WP recommendations
    https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/
    """
    sys_name = ""
    if "PNet" in btagger:
        sys_name = "particleNet"
    elif "RobustParT" in btagger:
        sys_name = "robustParticleTransformer"
    elif "DeepFlav" in btagger:
        sys_name = "deepJet"
    elif "UParT" in btagger:
        sys_name = "UParTAK4"

    cset = correctionlib.CorrectionSet.from_file(get_pog_json("btagging", year))
    btag_cut = b_taggers[year]["AK4"][btagger][wp]

    eff_file = f"{package_path}/hbb/data/btag/mc_eff_{btagger}_{year}.pkl"
    with open(eff_file, 'rb') as f:
        lookup_dict = pickle.load(f)

    eff_opt = "TTbar+QCD"  
        #options = "TTbar+QCD", "TTbar", "QCD"
        #defined in src/hbb/data/btag/compile_btag_eff.py

    def eff_lookup(x, y, z): return lookup_dict[eff_opt](x, y, z)

    jets_l = jets[(jets.hadronFlavour == 0) & (abs(jets.eta)<2.5)]
    jets_b = jets[(jets.hadronFlavour == 4) & (abs(jets.eta)<2.5)]
    jets_c = jets[(jets.hadronFlavour == 5) & (abs(jets.eta)<2.5)]

    pass_l = getattr(jets_l, btagger) > btag_cut
    pass_b = getattr(jets_b, btagger) > btag_cut
    pass_c = getattr(jets_c, btagger) > btag_cut

    eff_l = eff_lookup(jets_l.hadronFlavour, jets_l.pt, abs(jets_l.eta))
    eff_b = eff_lookup(jets_b.hadronFlavour, jets_b.pt, abs(jets_b.eta))
    eff_c = eff_lookup(jets_c.hadronFlavour, jets_c.pt, abs(jets_c.eta))

    def calc_weight(eff, sf, pass_tag):
        tagged = ak.prod(sf[pass_tag], axis=-1)
        untagged = ak.prod(((1 - sf*eff) / (1 - eff))[~pass_tag], axis=-1)
        return ak.fill_none(tagged * untagged, 1.)
    
    def get_sf(jets, j_flav, syst):
        j, nj = ak.flatten(jets), ak.num(jets)
        sf = cset[f"{sys_name}_{j_flav}"].evaluate(syst, wp, j.hadronFlavour, abs(j.eta), j.pt)
        return ak.unflatten(sf, nj)

    weight_l = calc_weight( eff_l, get_sf(jets_l, "light", "central"), pass_l )
    weight_b = calc_weight( eff_b, get_sf(jets_b, "comb", "central"), pass_b )
    weight_c = calc_weight( eff_c, get_sf(jets_c, "comb", "central"), pass_c )

    weights.add(f'{alt_str}btagLightSF', weight_l)
    weights.add(f'{alt_str}btagBSF', weight_b)
    weights.add(f'{alt_str}btagCSF', weight_c)
    
    nominal = weight_l * weight_b * weight_c

    weights.add(
        f"{alt_str}btagSFlight_{year}",
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_l, get_sf(jets_l, "light", "up"), pass_l),
        weightDown=calc_weight(eff_l, get_sf(jets_l, "light", "down"), pass_l),
    )
    weights.add(
        f"{alt_str}btagSFb_{year}",
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_b, get_sf(jets_b, "comb", "up"), pass_b),
        weightDown=calc_weight(eff_b, get_sf(jets_b, "comb", "down"), pass_b),
    )
    weights.add(
        f"{alt_str}btagSFc_{year}",
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_c, get_sf(jets_c, "comb", "up"), pass_c),
        weightDown=calc_weight(eff_c, get_sf(jets_c, "comb", "down"), pass_c),
    )
    weights.add(
        f'{alt_str}btagSFlight_correlated', 
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_l, get_sf(jets_l, "light", "up_correlated"), pass_l),
        weightDown=calc_weight(eff_l, get_sf(jets_l, "light", "down_correlated"), pass_l),
    )
    weights.add(
        f'{alt_str}btagSFb_correlated', 
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_b, get_sf(jets_b, "comb", "up_correlated"), pass_b),
        weightDown=calc_weight(eff_b, get_sf(jets_b, "comb", "down_correlated"), pass_b),
    )
    weights.add(
        f'{alt_str}btagSFc_correlated', 
        ak.ones_like(nominal),
        weightUp=calc_weight(eff_c, get_sf(jets_c, "comb", "up_correlated"), pass_c),
        weightDown=calc_weight(eff_c, get_sf(jets_c, "comb", "down_correlated"), pass_c),
    )
    return nominal

def add_muon_weights(weights: Weights, year: str, muons: MuonArray, pt_type: str, muon_type: str, alt_str: str):
    """
    Corrections for medium pt GeV muons
    https://muon-wiki.docs.cern.ch/guidelines/corrections/#medium-pt-30-gev-pt-200-gev

    Run 3 HLT Repo:
    https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run3
    """
    if muon_type == "loose":
        id_key = "NUM_LooseID_DEN_TrackerMuons"
        iso_key = "NUM_LoosePFIso_DEN_LooseID"
    elif muon_type == "highpt":
        id_key = "NUM_HighPtID_DEN_TrackerMuons"
        iso_key = "NUM_LooseRelTkIso_DEN_HighPtID"  #TODO double check

    #TODO add trigger SFs
    
    cset = correctionlib.CorrectionSet.from_file(get_pog_json("muon", year))
    m, nm = ak.flatten(muons), ak.num(muons, axis=1)
    
    def get_sf(cset_key, syst):
        sf = cset[cset_key].evaluate(abs(m.eta), getattr(m, pt_type), syst)
        return ak.prod(ak.unflatten(sf, nm), axis=-1)

    id_nom = get_sf(id_key, "nominal")
    id_up = get_sf(id_key, "systup") 
    id_down = get_sf(id_key, "systdown") 

    iso_nom = get_sf(iso_key, "nominal")
    iso_up = get_sf(iso_key, "systup")
    iso_down = get_sf(iso_key, "systdown")

    weights.add(f"{alt_str}muon_ID", id_nom, id_up, id_down)
    weights.add(f"{alt_str}muon_ISO", iso_nom, iso_up, iso_down)
    
    return

def add_photon_weights(weights: Weights, year: str, photons, alt_str: str):
    """
    Corrections for tight ID photons
    https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3
    """
    id_key = "Photon-ID-SF"

    year_map = {
        "2022" : "2022Re-recoBCD",
        "2022EE" : "2022Re-recoE+PromptFG",
        "2023" : "2023PromptC",
        "2023BPix" : "2023PromptD",
        "2024" : "2024Prompt",
    }

    cset = correctionlib.CorrectionSet.from_file(get_pog_json("photon", year))

    if "2023" in year:   
        #json format is different for 23 and 23BPix
        #https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3#Photon_JSON_format_AN1
        id_nom = cset[id_key].evaluate(year_map[year], "sf", "Tight", photons.eta, photons.pt, photons.phi)
        id_up = cset[id_key].evaluate(year_map[year], "sfup", "Tight", photons.eta, photons.pt, photons.phi)
        id_down = cset[id_key].evaluate(year_map[year], "sfdown", "Tight", photons.eta, photons.pt, photons.phi)
    else:
        id_nom = cset[id_key].evaluate(year_map[year], "sf", "Tight", photons.eta, photons.pt)
        id_up = cset[id_key].evaluate(year_map[year], "sfup", "Tight", photons.eta, photons.pt)
        id_down = cset[id_key].evaluate(year_map[year], "sfdown", "Tight", photons.eta, photons.pt)

    weights.add(f"{alt_str}photon_ID", id_nom, id_up, id_down)

    return

mupt_variations = {
    "MuonPTScale" : "ptscalecorr",
    "MuonPTRes" : "ptcorr_resol"
}

def correct_muons(muons, events, year: str, isRealData: bool):
    """
    Central corrections maintained by MUON POG
    https://muon-wiki.docs.cern.ch/guidelines/corrections/#medium-pt-scale-and-resolution
    https://gitlab.cern.ch/cms-muonPOG/muonscarekit
    src/hbb/MuonScaRe.py refactored to work with dask+awkward by Lara
    """
    cset = correctionlib.CorrectionSet.from_file(get_pog_json("muon_pt", year))

    if isRealData:
        muons["ptcorr"] = pt_scale(1, muons.pt, muons.eta, muons.phi, muons.charge, cset, nested=True)

    else:
        muons["ptscalecorr"] = pt_scale(0, muons.pt, muons.eta,  muons.phi, muons.charge, cset, nested=True)
        muons["ptcorr"] = pt_resol( muons.ptscalecorr, muons.eta, muons.phi, muons.nTrackerLayers, 
                                   events.event, events.luminosityBlock, cset, nested=True)

        muons["ptscalecorr_up"] = pt_scale_var(muons.ptcorr, muons.eta, muons.phi, muons.charge, "up", cset, nested=True)
        muons["ptscalecorr_down"] = pt_scale_var(muons.ptcorr, muons.eta, muons.phi, muons.charge, "dn", cset, nested=True)

        muons["ptcorr_resol_up"] = pt_resol_var(muons.ptscalecorr, muons.ptcorr, muons.eta, "up", cset, nested=True)
        muons["ptcorr_resol_down"] = pt_resol_var(muons.ptscalecorr, muons.ptcorr, muons.eta, "dn", cset, nested=True)

    return muons

egm_variations = {
    "ElectronPTScale": "ptcorr_scale",
    "ElectronPTRes": "ptcorr_smear",
    "PhotonPTScale": "ptcorr_scale",
    "PhotonPTRes": "ptcorr_smear"
}

def correct_egm_obj(egm_obj_array, egm_obj_type: str, year: str, isRealData: bool, nano_version: str, event_run):
    
    cset = correctionlib.CorrectionSet.from_file(get_pog_json(f"{egm_obj_type}SS", year))
    scale_ev = cset.compound["Scale"]
    smear_and_syst_ev = cset.compound["SmearAndSyst"]
    
    egm_obj_array["event_run"] = ak.broadcast_arrays(event_run, egm_obj_array.pt)[0]
    if "v15" in nano_version:
        egm_obj_array["sc_eta"] = egm_obj_array.superclusterEta
    else:
        if egm_obj_type == "electron":
            egm_obj_array["sc_eta"] = egm_obj_array.eta + egm_obj_array.deltaEtaSC
        elif egm_obj_type == "photon":
            egm_obj_array["sc_eta"] = egm_obj_array.eta
            
    e_flat, ne = ak.flatten(egm_obj_array), ak.num(egm_obj_array)
    
    def get_ss_sf(var):
        sf = smear_and_syst_ev.evaluate(var, e_flat.pt, e_flat.r9, e_flat.sc_eta)
        return ak.unflatten(sf, ne)
    
    if isRealData:
        scale_flat = scale_ev.avaluate("scale", e_flat.event_run, e_flat.sc_eta, e_flat.r9, e_flat.pt, e_flat.seedGain)
        scale = ak.unflatten(scale_flat, ne)
        
        egm_obj_array["ptcorr"] = egm_obj_array.pt * scale
        
    else:
        rng = np.random.default_rng(seed=125)
        random_numbers = rng.normal(loc=0.0, scale=1.0, size=len(egm_obj_array.pt))
        
        smearing = 1. + random_numbers * get_ss_sf("smear")
        smearing_up = 1. + random_numbers * get_ss_sf("smear_up")
        smearing_dn = 1. + random_numbers * get_ss_sf("smear_down")
        smearing_scale_up = 1. + random_numbers * get_ss_sf("scale_up")
        smearing_scale_dn = 1. + random_numbers * get_ss_sf("scale_down")
        
        egm_obj_array["ptcorr"] = egm_obj_array.pt * smearing
        
        egm_obj_array["ptcorr_smear_up"] = egm_obj_array.pt * smearing_up
        egm_obj_array["ptcorr_smear_dn"] = egm_obj_array.pt * smearing_dn
        
        egm_obj_array["ptcorr_scale_up"] = egm_obj_array.pt * smearing_scale_up
        egm_obj_array["ptcorr_scale_dn"] = egm_obj_array.pt * smearing_scale_dn
        
    return egm_obj_array

def add_VJets_corrections(weights: Weights, dataset: str, genpart):

    boson = ak.firsts(genpart[
            ((genpart.pdgId == 23)|(abs(genpart.pdgId) == 24))
            & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
        ])
    vpt = ak.fill_none(boson.pt, 0.)

    Zjets_pref = ["DYto2L-2Jets_MLL-50", "DYto2E-2Jets_MLL-50", "DYto2Mu-2Jets_MLL-50", "DYto2Tau-2Jets_MLL-50", "Zto2Q-4Jets_Bin-HT", "Zto2Q-4Jets_HT"]
    Wjets_pref = ["WtoLNu-2Jets", "Wto2Q-3Jets_Bin-HT", "Wto2Q-3Jets_HT"]
    isZ_dataset = any(ds in dataset for ds in Zjets_pref)
    isW_dataset = any(ds in dataset for ds in Wjets_pref)
    
    common_systs = [
        "d1K_NLO",
        "d2K_NLO",
        "d3K_NLO",
        "d1kappa_EW",
    ]
    zsysts = common_systs + [
        "Z_d2kappa_EW",
        "Z_d3kappa_EW",
    ]
    wsysts = common_systs + [
        "W_d2kappa_EW",
        "W_d3kappa_EW",
    ]

    def add_systs(systlist, qcdcorr, ewkcorr):
        ewknom = ewkcorr.evaluate("nominal", vpt)
        weights.add("vjets_nominal", qcdcorr * ewknom if qcdcorr is not None else ewknom)
        ones = ak.ones_like(vpt)
        for syst in systlist:
            weights.add(syst, ones, ewkcorr.evaluate(syst + "_up", vpt) / ewknom, ewkcorr.evaluate(syst + "_down", vpt) / ewknom)

    file = f"{package_path}/hbb/data/vjets/vjets_corrections_2flavorDY.json"
    vjets_kfactors = correctionlib.CorrectionSet.from_file(file)

    if isZ_dataset:
        qcdcorr = vjets_kfactors["Z_MLMtoFXFX"].evaluate(vpt)
        ewkcorr = vjets_kfactors["Z_FixedOrderComponent"]
        add_systs(zsysts, qcdcorr, ewkcorr)

    elif isW_dataset:
        qcdcorr = vjets_kfactors["W_MLMtoFXFX"].evaluate(vpt)
        ewkcorr = vjets_kfactors["W_FixedOrderComponent"]
        add_systs(wsysts, qcdcorr, ewkcorr)

    return