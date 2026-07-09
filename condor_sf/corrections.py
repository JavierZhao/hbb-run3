"""
POG JetID and jet veto map helpers for condor_sf scale-factor workflows.

Used when NanoAOD lacks ``jetId`` (e.g. v15 / 2024) and for Run-3 jet veto maps.
"""

import awkward as ak
import correctionlib
import numpy as np
from functools import lru_cache

POG_CORRECTION_PATH = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration"

POG_JSONS = {
    "jetid": ("JME", "jetid.json.gz"),
    "jetvetomap": ("JME", "jetvetomaps.json.gz"),
}

YEARS = {
    "2022": "2022_Summer22",
    "2022EE": "2022_Summer22EE",
    "2023": "2023_Summer23",
    "2023BPix": "2023_Summer23BPix",
    "2024": "2024_Summer24",
}

JET_VETOMAP_KEYS = {
    "2022": "Summer22_23Sep2023_RunCD_V1",
    "2022EE": "Summer22EE_23Sep2023_RunEFG_V1",
    "2023": "Summer23Prompt23_RunC_V1",
    "2023BPix": "Summer23BPixPrompt23_RunD_V1",
    "2024": "Summer24Prompt24_RunBCDEFGHI_V1",
}


def get_pog_json(obj, year):
    if obj not in POG_JSONS:
        raise KeyError(f"Unsupported correction object: {obj}")
    if year not in YEARS:
        raise KeyError(f"Unsupported year: {year}")
    group, filename = POG_JSONS[obj]
    return f"{POG_CORRECTION_PATH}/POG/{group}/{YEARS[year]}/{filename}"


@lru_cache(maxsize=32)
def _load_correctionset(obj, year):
    return correctionlib.CorrectionSet.from_file(get_pog_json(obj, year))


def _jetid_bit_mask(jets, bit_name, fallback_mask):
    bit = getattr(jets, bit_name, None)
    if bit is None:
        return fallback_mask
    return 1 << bit


def correct_jetid(jets, jet_type, year):
    """Add ``jetidtight`` and ``jetidtightlepveto`` via POG jetid.json.gz."""
    evaluator = _load_correctionset("jetid", year)

    if jet_type == "AK8":
        name_tight = "AK8PUPPI_Tight"
        name_tightlv = "AK8PUPPI_TightLeptonVeto"
    elif jet_type == "AK4":
        name_tight = "AK4PUPPI_Tight"
        name_tightlv = "AK4PUPPI_TightLeptonVeto"
    else:
        raise ValueError(f"Unsupported jet_type: {jet_type}. Expected 'AK4' or 'AK8'.")

    flat_j = ak.flatten(jets)
    counts = ak.num(jets)

    def get_jetid(j, nj, corr_name):
        def real_col(x):
            x = ak.fill_none(x, 0.0)
            a = np.asarray(ak.to_numpy(x), dtype=np.float64, order="C")
            np.nan_to_num(a, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
            return a

        def int_col(x):
            x = ak.fill_none(x, 0)
            return np.asarray(ak.to_numpy(x), dtype=np.int64, order="C")

        eta = np.clip(real_col(j.eta), -4.7, 4.7)
        chHEF = np.clip(real_col(j.chHEF), 0.0, 1.0)
        neHEF = np.clip(real_col(j.neHEF), 0.0, 1.0)
        chEmEF = np.clip(real_col(j.chEmEF), 0.0, 1.0)
        neEmEF = np.clip(real_col(j.neEmEF), 0.0, 1.0)
        muEF = np.clip(real_col(j.muEF), 0.0, 1.0)
        chMult = np.clip(int_col(j.chMultiplicity), 0, 100)
        neMult = np.clip(int_col(j.neMultiplicity), 0, 100)
        multiplicity = chMult + neMult

        if multiplicity.size == 0:
            return ak.unflatten(np.array([], dtype=np.float64), nj)

        corr = evaluator[corr_name]
        try:
            jetid = np.asarray(
                corr.evaluate(
                    eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMult, neMult, multiplicity
                ),
                dtype=np.float64,
            )
        except RuntimeError:
            n = int(multiplicity.size)
            jetid = np.empty(n, dtype=np.float64)
            for i in range(n):
                jetid[i] = corr.evaluate(
                    float(eta[i]),
                    float(chHEF[i]),
                    float(neHEF[i]),
                    float(chEmEF[i]),
                    float(neEmEF[i]),
                    float(muEF[i]),
                    int(chMult[i]),
                    int(neMult[i]),
                    int(multiplicity[i]),
                )
        return ak.unflatten(jetid, nj)

    jets["jetidtight"] = ak.values_astype(get_jetid(flat_j, counts, name_tight), bool)
    jets["jetidtightlepveto"] = ak.values_astype(get_jetid(flat_j, counts, name_tightlv), bool)
    return jets


def get_minimal_jet_mask(jets):
    """
    Run-3 minimal AK4 jet selection for jet veto maps:
      - pt > 15 GeV
      - tightLepVeto JetID
      - (chEmEF + neEmEF) < 0.9
    """
    if "jetId" in ak.fields(jets):
        mask_tightlv = _jetid_bit_mask(jets, "TIGHTLEPVETO", 4)
        pass_id = (jets.jetId & mask_tightlv) != 0
    else:
        pass_id = jets.jetidtightlepveto

    em_frac = ak.fill_none(jets.chEmEF, 0.0) + ak.fill_none(jets.neEmEF, 0.0)
    return (jets.pt > 15.0) & pass_id & (em_frac < 0.9)


def event_passes_jet_veto_map(jets, year, map_type="jetvetomap"):
    """
    Event-level mask: True if the event passes (no minimal jet in a veto region).
    """
    if year not in JET_VETOMAP_KEYS:
        raise KeyError(f"Unsupported year for jet veto map: {year}")

    corr = _load_correctionset("jetvetomap", year)[JET_VETOMAP_KEYS[year]]
    selected = jets[get_minimal_jet_mask(jets)]
    counts = ak.num(selected)
    flat = ak.flatten(selected)

    if len(flat) == 0:
        return np.ones(len(jets), dtype=bool)

    eta = np.asarray(ak.to_numpy(ak.fill_none(flat.eta, 0.0)), dtype=np.float64)
    phi = np.asarray(ak.to_numpy(ak.fill_none(flat.phi, 0.0)), dtype=np.float64)
    try:
        veto_val = np.asarray(corr.evaluate(map_type, eta, phi), dtype=np.float64)
    except (ValueError, RuntimeError):
        veto_val = np.empty(len(eta), dtype=np.float64)
        for i in range(len(eta)):
            veto_val[i] = float(corr.evaluate(map_type, float(eta[i]), float(phi[i])))

    in_veto = ak.unflatten(veto_val > 0.0, counts)
    has_veto_jet = ak.fill_none(ak.any(in_veto, axis=1), False)
    return ~ak.to_numpy(has_veto_jet)
