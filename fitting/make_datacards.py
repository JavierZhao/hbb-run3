"""
Datacard Maker - Fully Configuration-Driven
Supports: VBF Hbb Analysis, ZGamma Validation Region

Author(s): Gabi Hamilton, Lara Zygala, Cristina Mantilla
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import warnings
from pathlib import Path

import numpy as np
import rhalphalib as rl
import ROOT
from card_utils import (
    add_systematics,
    badtemplate,
    get_merged_template,
    get_template,
    one_bin,
    plot_mctf,
    safe_ratio,
)
from template_utils import (
    year_systs,
    Zjets_thsysts,
    Wjets_thsysts,
    sig_th_systs
)
from scalesmear import MorphHistW2

from hbb.common_vars import LUMI

ROOT.gROOT.SetBatch(True)
warnings.filterwarnings("ignore")
rl.util.install_roofit_helpers()

lumi_err = {"2022": 1.01, "2023": 1.02, "2024": 1.02}  # 2024: TODO get official CMS value
eps = 0.001


def rhalphabet(args):
    # ---------------------------------------------------------
    # 1. SETUP & LOAD CONFIG
    # ---------------------------------------------------------
    year = args.year
    tag = args.tag
    analysis = args.analysis
    print(f"Running Card Maker for {year} (Analysis: {analysis})")

    # Define Paths
    working_dir = Path(f"{args.outdir}/{tag}/{year}")
    datacard_dir = working_dir / "datacards"
    initvals_dir = working_dir / "initial_vals"

    datacard_dir.mkdir(parents=True, exist_ok=True)
    initvals_dir.mkdir(parents=True, exist_ok=True)

    # Load Configuration
    json_name = f"setup_{analysis}.json"
    with Path(json_name).open() as f:
        config = json.load(f)

    # ---------------------------------------------------------
    # 2. READ SETTINGS FROM JSON
    # ---------------------------------------------------------

    # Files & naming
    root_file_name = config.get("root_filename", "signalregion.root").replace("{year}", year)

    # Define infile_path here!
    infile_path = Path(args.indir) / root_file_name if args.indir else working_dir / root_file_name

    data_obs_name = config.get("data_obs_name", "Jetdata")
    qcd_tf_proc = config.get("qcd_proc", "QCD")
    pt_min_scale = config.get("pt_min_scale", 450.0)
    regions_to_fit = config.get("regions_to_fit", ["bb"])
    rho_scaling_max = config.get("rho_scaling_max", -2.1)

    # Process Definitions
    sample_dict = config.get("process_groups")
    for group in sample_dict.values():
        group["components"] = [tuple(c) for c in group["components"]]

    # Jet Mass Scale and Resolution
    # shift the mass axis up/down by JMSR_SCALE GeV,
    JMSR_SCALE = config.get("jmsr_scale", 1.0)  # read from JSON or hardcode
    # broaden/narrow the distribution by JMSR_SMEAR fraction
    JMSR_SMEAR = config.get("jmsr_smear", 0.1)  # read from JSON or hardcode
    jmsr_processes = config.get("jmsr_processes", [])  # e.g. ["zgammabb"] in JSON

    # ---------------------------------------------------------
    # 3. OBSERVABLES & SYSTEMATICS
    # ---------------------------------------------------------
    msd_cfg = config["observable"]
    msdbins = np.linspace(msd_cfg["min"], msd_cfg["max"], msd_cfg["nbins"] + 1)
    msd = rl.Observable(msd_cfg["name"], msdbins)

    cats_cfg = config["categories"]
    cats = list(cats_cfg.keys())
    if "mucr" in cats:
        cats.remove("mucr")  # necessary since all categories get the same treatment

    # TT Independent Parameters
    tqqeffSF = rl.IndependentParameter(f"tqqeffSF_{year}", 1.0, -50, 50)
    tqqeffBCSF = rl.IndependentParameter(f"tqqeffBCSF_{year}", 1.0, -50, 50)
    tqqnormSF = rl.IndependentParameter(f"tqqnormSF_{year}", 1.0, -50, 50)

    do_muon_CR = config.get("do_muon_CR", False)

    # Standard Luminosity Uncertainty
    sys_lumi_uncor = rl.NuisanceParameter(f"CMS_lumi_13TeV_{year[:4]}", "lnN")

    do_systematics = config.get("do_systematics", False)
    syst_map = {}

    if do_systematics:
        # --- A. Experimental Systematics (from sys_dict) ---
        available_exp_systs = {
            "pileup": rl.NuisanceParameter(f"CMS_PU_{year}", "lnN"),
            "JES": rl.NuisanceParameter(f"CMS_scale_j_{year}", "lnN"),
            "JER": rl.NuisanceParameter(f"CMS_res_j_{year}", "lnN"),
            "UES": rl.NuisanceParameter(f"CMS_ues_j_{year}", "lnN"),
            "MuonPTScale": rl.NuisanceParameter(f"CMS_scale_m_{year}", "lnN"),
            "MuonPTRes": rl.NuisanceParameter(f"CMS_res_m_{year}", "lnN"),
            f"btagSFb_{year}": rl.NuisanceParameter(f"CMS_btagSFb_{year}", "lnN"),
            f"btagSFc_{year}": rl.NuisanceParameter(f"CMS_btagSFc_{year}", "lnN"),
            f"btagSFlight_{year}": rl.NuisanceParameter(f"CMS_btagSFlight_{year}", "lnN"),
            "btagSFb_correlated": rl.NuisanceParameter(f"CMS_btagSFb_correlated_{year}", "lnN"),
            "btagSFc_correlated": rl.NuisanceParameter(f"CMS_btagSFc_correlated_{year}", "lnN"),
            "btagSFlight_correlated": rl.NuisanceParameter(f"CMS_btagSFlight_correlated_{year}", "lnN"),
            "JMSunconstrained": rl.NuisanceParameter(f"CMS_jms_{year}", "shapeU", lo=-5, hi=5),
            "JMRunconstrained": rl.NuisanceParameter(f"CMS_jmr_{year}", "shapeU", lo=-5, hi=5),
            "JMS": rl.NuisanceParameter(f"CMS_jms_{year}", "shape"),
            "JMR": rl.NuisanceParameter(f"CMS_jmr_{year}", "shape"),
        }

        # --- B. Theory Systematics (PDF, Scale, ISR/FSR) ---
        theory_systs = {}

        sig_mode = ["ggF", "VBF", "VH", "ttH"]

        for mode in sig_mode:
            for th in sig_th_systs:
                theory_systs[f"{th}_{mode}"] = rl.NuisanceParameter(f"{th}_{mode}", 'lnN') 
        
        for sys in set(Zjets_thsysts + Wjets_thsysts):
            theory_systs[sys] = rl.NuisanceParameter(f"CMS_hbb_{sys}", 'lnN') 


        # Combine all available systematics into one map
        all_available = {**available_exp_systs, **theory_systs}

        # Pull active list from JSON config
        active_list = config.get("active_systematics", [])

        def sys_check(sys_name):
            if sys_name in all_available:
                syst_map[sys_name] = all_available[sys_name]
            else:
                print(f"Warning: Systematic {sys_name} requested in JSON but not defined in script.")

        #Check all systematic special groupings (held in template_utils.py)
        for name in active_list:
            if name == "VJets":
                for sys in set(Zjets_thsysts + Wjets_thsysts):
                    sys_check(sys)
            elif name in year_systs:
                sys_check(f"{sys}_{year}")
            elif name in sig_th_systs:
                for mode in sig_mode:
                    sys_check(f"{name}_{mode}")
            else:
                sys_check(name)
    # ---------------------------------------------------------
    # 4. QCD ESTIMATION LOOP
    # ---------------------------------------------------------
    tf_params = {}
    validbins = {}

    for cat in cats:
        if "bins_pt" in cats_cfg[cat]:
            ptbins = np.array(cats_cfg[cat]["bins_pt"])
        else:
            ptbins = np.array(cats_cfg[cat]["bins"])
        npt = len(ptbins) - 1

        # Grid Setup
        ptpts, msdpts = np.meshgrid(
            ptbins[:-1] + 0.3 * np.diff(ptbins),
            msdbins[:-1] + 0.5 * np.diff(msdbins),
            indexing="ij",
        )
        rhopts = 2 * np.log(msdpts / ptpts)
        # per-category pT range and rho max ("pt_min"/"pt_max"/"rho_scaling_max"), else global
        cat_cfg = cats_cfg[cat]
        cat_pt_min = float(cat_cfg.get("pt_min", pt_min_scale))
        cat_pt_max = float(cat_cfg.get("pt_max", 1200.0))
        cat_rho_max = float(cat_cfg.get("rho_scaling_max", rho_scaling_max))
        ptscaled = (ptpts - cat_pt_min) / (cat_pt_max - cat_pt_min)
        rhoscaled = (rhopts - (-6.0)) / (cat_rho_max - (-6.0))

        validbins[cat] = (rhoscaled >= 0.0) & (rhoscaled <= 1.0)
        rhoscaled[~validbins[cat]] = 1

        tf_params[cat] = {}
        fitfailed_qcd = {}

        for reg in regions_to_fit:
            fitfailed_qcd[reg] = 0
            while fitfailed_qcd[reg] < 5:
                qcdmodel = rl.Model(f"qcdmodel_{cat}_{reg}")
                qcdpass, qcdfail = 0.0, 0.0

                for ptbin in range(npt):
                    binindex = ptbin
                    if analysis == "vbf" and "hi" in cat:
                        binindex = 1

                    failCh = rl.Channel(f"ptbin{ptbin}{cat}fail{year}{reg}")
                    passCh = rl.Channel(f"ptbin{ptbin}{cat}pass{year}{reg}")
                    qcdmodel.addChannel(failCh)
                    qcdmodel.addChannel(passCh)

                    failTempl = get_template(
                        infile_path, qcd_tf_proc, "fail_", binindex + 1, cat, msd, "nominal"
                    )
                    passTempl = get_template(
                        infile_path, qcd_tf_proc, f"pass_{reg}_", binindex + 1, cat, msd, "nominal"
                    )

                    failCh.setObservation(failTempl, read_sumw2=True)
                    passCh.setObservation(passTempl, read_sumw2=True)
                    qcdfail += failCh.getObservation()[0].sum()
                    qcdpass += passCh.getObservation()[0].sum()

                qcdeff = qcdpass / qcdfail
                print(f"Inclusive P/F ({cat} {reg}) = {qcdeff:.4f}")

                # Initial Values Loading
                initF = initvals_dir / f"initial_vals_{cat}_{reg}.json"
                initial_vals = None
                mc_pt_order = cats_cfg[cat]["tfmc_order"][reg]["pt"]
                mc_rho_order = cats_cfg[cat]["tfmc_order"][reg]["rho"]
                if initF.exists():
                    with initF.open() as f:
                        loaded = np.array(json.load(f)["initial_vals"])
                    if (loaded.shape[0] - 1 == mc_pt_order) and (
                        loaded.shape[1] - 1 == mc_rho_order
                    ):
                        initial_vals = loaded
                    else:
                        print(f"Order Mismatch for {reg}. Resetting.")

                if initial_vals is None:
                    initial_vals = np.full((mc_pt_order + 1, mc_rho_order + 1), 1.0)

                tf_MCtempl = rl.BasisPoly(
                    f"tf_MCtempl_{cat}{reg}{year}",
                    (mc_pt_order, mc_rho_order),
                    ["pt", "rho"],
                    basis="Bernstein",
                    init_params=initial_vals,
                    limits=(0, 10),
                )

                tf_MCtempl_params = qcdeff * tf_MCtempl(ptscaled, rhoscaled)

                for ptbin in range(npt):
                    failCh = qcdmodel[f"ptbin{ptbin}{cat}fail{year}{reg}"]
                    passCh = qcdmodel[f"ptbin{ptbin}{cat}pass{year}{reg}"]

                    failObs = failCh.getObservation()[0]
                    qcdparams = np.array(
                        [
                            rl.IndependentParameter(f"qcdparam_ptbin{ptbin}{cat}{year}{reg}_{i}", 0)
                            for i in range(msd.nbins)
                        ]
                    )
                    scaledparams = (
                        failObs * (1 + 10.0 / np.maximum(1.0, np.sqrt(failObs))) ** qcdparams
                    )

                    fail_qcd = rl.ParametericSample(
                        f"ptbin{ptbin}{cat}fail{year}{reg}_qcd",
                        rl.Sample.BACKGROUND,
                        msd,
                        scaledparams,
                    )
                    failCh.addSample(fail_qcd)

                    pass_qcd = rl.TransferFactorSample(
                        f"ptbin{ptbin}{cat}pass{year}{reg}_qcd",
                        rl.Sample.BACKGROUND,
                        tf_MCtempl_params[ptbin, :],
                        fail_qcd,
                    )
                    passCh.addSample(pass_qcd)

                    failCh.mask = validbins[cat][ptbin]
                    passCh.mask = validbins[cat][ptbin]

                # Fit
                qcdfit_ws = ROOT.RooWorkspace("w")
                simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
                qcdfit = simpdf.fitTo(
                    obs,
                    ROOT.RooFit.Extended(True),
                    ROOT.RooFit.SumW2Error(True),
                    ROOT.RooFit.Strategy(2),
                    ROOT.RooFit.Save(),
                    ROOT.RooFit.PrintLevel(-1),
                )

                if qcdfit.status() != 0:
                    #want to save the values every time so that you don't start from scratch next time you rerun the script.
                    fitfailed_qcd[reg] += 1
                    allparams = dict(zip(qcdfit.nameArray(), qcdfit.valueArray()))
                    pvalues = [allparams[p.name] for p in tf_MCtempl.parameters.reshape(-1)]
                    new_values = np.array(pvalues).reshape(tf_MCtempl.parameters.shape)
                    with initF.open("w") as outfile:
                        json.dump({"initial_vals": new_values.tolist()}, outfile)
                else:
                    break

            if fitfailed_qcd[reg] >= 5:
                print(f"\n[FIT] All 5 attempts failed for {cat} {reg}.")
                print(f"  Last status={qcdfit.status()}, covQual={qcdfit.covQual()}")
                print(f"  covQual meanings: -1=not calc, 0=not pos-def, 1=forced pos-def, 2=approx, 3=full accurate")
                print(f"  MC template order: pt={mc_pt_order}, rho={mc_rho_order}")
                print(f"  Inclusive P/F = {qcdeff:.4f} — if very small, model may be overparameterized")
                raise RuntimeError(f"Could not fit QCD for {cat} {reg} after 5 tries!")

            plot_mctf(
                tf_MCtempl,
                msdbins,
                f"{cat}_{reg}",
                year,
                tag,
                str(working_dir),
                pt_min=cat_pt_min,  # <--- Pass from config
                rho_max=cat_rho_max,  # <--- Pass from config (-1.0 for ZG)
                pt_max=cat_pt_max,
            )

            param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
            decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(
                tf_MCtempl.name + "_deco", qcdfit, param_names
            )
            tf_MCtempl.parameters = decoVector.correlated_params.reshape(
                tf_MCtempl.parameters.shape
            )

            # Residual
            res_pt_order = cats_cfg[cat]["tfres_order"][reg]["pt"]
            res_rho_order = cats_cfg[cat]["tfres_order"][reg]["rho"]
            resid_init = np.full((res_pt_order + 1, res_rho_order + 1), 1.0)
            tf_dataResidual = rl.BasisPoly(
                f"tf_dataResidual_{year}{cat}{reg}",
                (res_pt_order, res_rho_order),
                ["pt", "rho"],
                basis="Bernstein",
                init_params=resid_init,
                limits=(0, 20),
            )
            tf_params[cat][reg] = (
                qcdeff * tf_MCtempl(ptscaled, rhoscaled) * tf_dataResidual(ptscaled, rhoscaled)
            )

    # ---------------------------------------------------------
    # 5. MAIN MODEL BUILDING
    # ---------------------------------------------------------
    pm_config = config.get(
        "physics_model",
        {
            "main": {
                "model": "HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel",
                "maps": ["map=.*/.*:r[1,-20,20]"],
            }
        },
    )
    model_dict = {}
    for model_name in pm_config:
        model_dict[model_name] = rl.Model(f"{analysis}_{model_name}Model_{year}")
        alt_model_pnames = pm_config[model_name].get("alt_pnames", {})

        for cat in cats:
            if "bins_pt" in cats_cfg[cat]:
                ptbins = np.array(cats_cfg[cat]["bins_pt"])
            else:
                ptbins = np.array(cats_cfg[cat]["bins"])

            for ptbin in range(len(ptbins) - 1):
                binindex = ptbin
                # Handle the VBF hi/lo binning logic
                if analysis == "vbf" and "hi" in cat:
                    binindex = 1

                regions = [f"pass_{r}_" for r in regions_to_fit] + ["fail_"]

                for region in regions:
                    ch_name = f"ptbin{ptbin}{cat}{region.replace('_', '')}{year}"
                    ch = rl.Channel(ch_name)
                    model_dict[model_name].addChannel(ch)

                    for proc_name, info in sample_dict.items():
                        # proc_name is e.g., 'ggF', 'VBF', 'ttbar'
                        # this is a (sumw, edges, name, sumw2) tuple
                        templ = get_merged_template(
                            infile_path, info["components"], region, binindex + 1, cat, msd
                        )
                        nominal = templ[0]

                        if badtemplate(nominal):
                            print(
                                f"Warning: Skipping template for {proc_name} in {ch_name} (failed badtemplate check)"
                            )
                            continue

                        stype = rl.Sample.SIGNAL if info["is_signal"] else rl.Sample.BACKGROUND
                        datacard_pname = proc_name if proc_name not in alt_model_pnames else alt_model_pnames[proc_name]
                        sample = rl.TemplateSample(ch.name + "_" + datacard_pname, stype, templ)

                        # Apply Luminosity Uncertainty
                        sample.setParamEffect(
                            sys_lumi_uncor, lumi_err[year[:4]] ** (LUMI[year[:4]] / LUMI["2022-2024"])
                        )

                        if do_systematics:
                            # 1. Automatic MC Statistical Uncertainties (Barlow-Beeston Lite)
                            # (Already handled inside add_systematics in card_utils.py)

                            # 2. Experimental Systematics (Shapes from ROOT file)
                            # Filter out specific systematics so they aren't double-applied
                            exp_syst_map = {
                                k: v
                                for k, v in syst_map.items()
                                if not any(ts in k for ts in (sig_th_systs + ["JMS","JMR"] + Zjets_thsysts + Wjets_thsysts))
                            }

                            add_systematics(
                                sample,
                                nominal,
                                exp_syst_map,
                                infile_path,
                                info["components"],
                                region,
                                binindex + 1,
                                cat,
                                msd,
                            )

                            # JMS/R systematics
                            jmsr_syst_map = {
                                k: v for k, v in syst_map.items() if k.startswith(("JMS", "JMR"))
                            }

                            if proc_name in jmsr_processes:
                                # Build a MorphHistW2 from the already-loaded nominal template.
                                sumw, edges, _name, sumw2 = templ
                                morph = MorphHistW2((sumw, edges, sumw2))

                                # Morphed Nominal
                                morphed_nom, _, _ = morph.get(shift=0.0, scale=1.0)

                                # Each pair shares the same morph — only one variant should be active at a time.
                                jmsr_pairs = [
                                    (
                                        "JMS",
                                        "JMSunconstrained",
                                        morph.get(shift=+JMSR_SCALE),
                                        morph.get(shift=-JMSR_SCALE),
                                    ),
                                    (
                                        "JMR",
                                        "JMRunconstrained",
                                        morph.get(scale=1.0 + JMSR_SMEAR),
                                        morph.get(scale=1.0 - JMSR_SMEAR),
                                    ),
                                ]

                                for key_a, key_b, (up_vals, _, _), (dn_vals, _, _) in jmsr_pairs:
                                    active_key = next(
                                        (k for k in (key_a, key_b) if k in jmsr_syst_map), None
                                    )
                                    if active_key is None:
                                        continue
                                    print(
                                        f"Adding {jmsr_syst_map[active_key].name} to {proc_name}, {region}"
                                    )

                                    if not np.allclose(morphed_nom, nominal, rtol=1e-3, atol=1e-6):
                                        print(
                                            f"  Warning: MorphHistW2 nominal mismatch for {proc_name} in {ch_name}. "
                                            f"Max relative diff: {np.max(np.abs(morphed_nom - nominal) / np.maximum(nominal, 1e-10)):.4f}"
                                        )

                                    sample.setParamEffect(
                                        jmsr_syst_map[active_key],
                                        safe_ratio(
                                            up_vals, morphed_nom
                                        ),  # take the ratio between up and nominal
                                        safe_ratio(dn_vals, morphed_nom),
                                        scale=1,  # this is a rescaling effect,
                                        # most useful for shape effects where the nuisance parameter effect
                                        # needs to be magnified to ensure good vertical interpolation
                                        # it is better left at 1, to avoid confusion later
                                    )

                            # 3. Theory Systematics (Process-Specific Logic)

                            # --- V+Jets ---
                            if proc_name in ["Wjets"]:
                                for s_name in set(Wjets_thsysts):
                                    s_up = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Up",
                                    )[0]
                                    s_do = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Down",
                                    )[0]
                                    # Look up using the mapped name (e.g., pdf_VH)
                                    syst_obj = syst_map.get(s_name)
                                    if syst_obj:
                                        sample.setParamEffect(
                                            syst_obj,
                                            np.sum(s_up) / np.sum(nominal),
                                            np.sum(s_do) / np.sum(nominal),
                                        )
                            if proc_name in ["Zjets", "Zjetsbb", "Zjetsc", "Zjetslight"]:
                                for s_name in set(Zjets_thsysts):
                                    s_up = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Up",
                                    )[0]
                                    s_do = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Down",
                                    )[0]
                                    # Look up using the mapped name (e.g., pdf_VH)
                                    syst_obj = syst_map.get(s_name)
                                    if syst_obj:
                                        sample.setParamEffect(
                                            syst_obj,
                                            np.sum(s_up) / np.sum(nominal),
                                            np.sum(s_do) / np.sum(nominal),
                                        )

                            # --- Higgs Signal Theory (PDF, ISR/FSR, Scale) ---
                            if any(s in proc_name for s in ["ggF", "VBF", "WH", "ZH", "ggZH", "ttH"]):
                                # Mapping logic: if proc is WH/ZH/ggZH, use "VH" for the nuisance name
                                # Extended to the other signal processes, since we separated them by bb/cc truth generated mode
                                proc_map_name = proc_name
                                if any(s in proc_name for s in ["WH", "ZH", "ggZH"]):
                                    proc_map_name = "VH" 
                                elif any(s in proc_name for s in ["ggF"]):
                                    proc_map_name = "ggF" 
                                elif any(s in proc_name for s in ["VBF"]):
                                    proc_map_name = "VBF" 

                                for s_name in [ "pdf_Higgs", "FSRPartonShower", "ISRPartonShower", ]:
                                    s_up = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Up",
                                    )[0]
                                    s_do = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst=f"{s_name}Down",
                                    )[0]

                                    # Look up using the mapped name (e.g., pdf_VH)
                                    syst_obj = syst_map.get(f"{s_name}_{proc_map_name}")
                                    if syst_obj:
                                        sample.setParamEffect(
                                            syst_obj,
                                            np.sum(s_up) / np.sum(nominal),
                                            np.sum(s_do) / np.sum(nominal),
                                        )

                                # ggF specific Scale (7pt)
                                if proc_map_name in ["ggF", "ttH"]:
                                    sc_up = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst="scalevar7ptUp",
                                    )[0]
                                    sc_do = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst="scalevar7ptDown",
                                    )[0]
                                    sample.setParamEffect(
                                        syst_map[f"QCDScale_{proc_map_name}"],
                                        np.sum(sc_up) / np.sum(nominal),
                                        np.sum(sc_do) / np.sum(nominal),
                                    )
                                elif proc_map_name in ["VBF", "VH"]:
                                    sc_up = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst="scalevar3ptUp",
                                    )[0]
                                    sc_do = get_merged_template(
                                        infile_path,
                                        info["components"],
                                        region,
                                        binindex + 1,
                                        cat,
                                        msd,
                                        syst="scalevar3ptDown",
                                    )[0]
                                    sample.setParamEffect(
                                        syst_map[f"QCDScale_{proc_map_name}"],
                                        np.sum(sc_up) / np.sum(nominal),
                                        np.sum(sc_do) / np.sum(nominal),
                                    )

                        ch.addSample(sample)

                    # Data
                    data_obs = get_template(
                        infile_path, data_obs_name, region, binindex + 1, cat, msd, syst="nominal"
                    )
                    ch.setObservation(data_obs[0:3])

        # ---------------------------------------------------------
        # 6. ADD DATA-DRIVEN QCD
        # ---------------------------------------------------------
        print("Adding Data-Driven QCD Background...")
        for cat in cats:
            if "bins_pt" in cats_cfg[cat]:
                ptbins = np.array(cats_cfg[cat]["bins_pt"])
            else:
                ptbins = np.array(cats_cfg[cat]["bins"])

            for ptbin in range(len(ptbins) - 1):
                failCh = model_dict[model_name][f"ptbin{ptbin}{cat}fail{year}"]

                initial_qcd = failCh.getObservation().astype(float)
                for sample in failCh:
                    initial_qcd -= sample.getExpectation(nominal=True)
                initial_qcd[initial_qcd < 0] = 0

                if args.debug:
                    print("FAIL QCD", failCh.name, "min", initial_qcd.min(), "zeros", np.sum(initial_qcd == 0))
                    print(initial_qcd)

                qcdparams = np.array(
                    [
                        rl.IndependentParameter(f"qcdparam_ptbin{ptbin}{cat}{year}_{i}", 0)
                        for i in range(msd.nbins)
                    ]
                )
                scaledparams = (
                    initial_qcd * (1 + 10.0 / np.maximum(1.0, np.sqrt(initial_qcd))) ** qcdparams
                )
                fail_qcd = rl.ParametericSample(
                    f"ptbin{ptbin}{cat}fail{year}_qcd", rl.Sample.BACKGROUND, msd, scaledparams
                )
                failCh.addSample(fail_qcd)

                # Add QCD to all pass regions defined in JSON
                for reg in regions_to_fit:
                    passCh = model_dict[model_name][f"ptbin{ptbin}{cat}pass{reg}{year}"]
                    pass_qcd = rl.TransferFactorSample(
                        f"ptbin{ptbin}{cat}pass{reg}{year}_qcd",
                        rl.Sample.BACKGROUND,
                        tf_params[cat][reg][ptbin, :],
                        fail_qcd,
                    )
                    passCh.addSample(pass_qcd)

                if do_muon_CR:
                    passChbb = model_dict[model_name][f"ptbin{ptbin}{cat}passbb{year}"]
                    passChcc = model_dict[model_name][f"ptbin{ptbin}{cat}passcc{year}"]

                    tqqpassbb = passChbb["ttbar"]
                    tqqpasscc = passChcc["ttbar"]
                    tqqfail = failCh["ttbar"]

                    sumPass = (
                        tqqpassbb.getExpectation(nominal=True).sum()
                        + tqqpasscc.getExpectation(nominal=True).sum()
                    )
                    sumFail = tqqfail.getExpectation(nominal=True).sum()

                    sumPassbb = tqqpassbb.getExpectation(nominal=True).sum()
                    sumPasscc = tqqpasscc.getExpectation(nominal=True).sum()

                    if any(s.name == f'ptbin{ptbin}{cat}passbb{year}_singlet' for s in passChbb.samples) or any(s.name == f'ptbin{ptbin}{cat}passcc{year}_singlet' for s in passChcc.samples):
                        stqqpassbb = passChbb["singlet"]
                        stqqpasscc = passChcc["singlet"]
                        stqqfail = failCh["singlet"]

                        sumPass += stqqpassbb.getExpectation(nominal=True).sum()
                        sumPass += stqqpasscc.getExpectation(nominal=True).sum()

                        sumPassbb += stqqpassbb.getExpectation(nominal=True).sum()
                        sumPasscc += stqqpasscc.getExpectation(nominal=True).sum()

                        sumFail += stqqfail.getExpectation(nominal=True).sum()

                        tqqPF = sumPass / sumFail
                        tqqBC = sumPassbb / sumPasscc

                        stqqpassbb.setParamEffect(tqqeffSF, 1 * tqqeffSF)
                        stqqpasscc.setParamEffect(tqqeffSF, 1 * tqqeffSF)
                        stqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)

                        stqqpassbb.setParamEffect(tqqeffBCSF, 1 * tqqeffBCSF)
                        stqqpasscc.setParamEffect(tqqeffBCSF, (1 - tqqeffBCSF) * tqqBC + 1)

                        stqqpassbb.setParamEffect(tqqnormSF, 1 * tqqnormSF)
                        stqqpasscc.setParamEffect(tqqnormSF, 1 * tqqnormSF)
                        stqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

                    tqqPF = sumPass / sumFail
                    tqqBC = sumPassbb / sumPasscc

                    tqqpassbb.setParamEffect(tqqeffSF, 1 * tqqeffSF)
                    tqqpasscc.setParamEffect(tqqeffSF, 1 * tqqeffSF)
                    tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)

                    tqqpassbb.setParamEffect(tqqeffBCSF, 1 * tqqeffBCSF)
                    tqqpasscc.setParamEffect(tqqeffBCSF, (1 - tqqeffBCSF) * tqqBC + 1)

                    tqqpassbb.setParamEffect(tqqnormSF, 1 * tqqnormSF)
                    tqqpasscc.setParamEffect(tqqnormSF, 1 * tqqnormSF)
                    tqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

        muonCR_model = rl.Model("muonCR_" + year)
        if do_muon_CR:
            templates = {}
            samps = ["QCD", "ttbar", "singlet", "Wjets", "Zjetsc", "Zjetslight", "Zjetsbb"]
            for region in ["pass_bb_", "pass_cc_", "fail_"]:

                ch_name = f"muonCR{region.replace('_', '')}{year}"

                ch = rl.Channel(ch_name)
                muonCR_model.addChannel(ch)
                for sName in samps:
                    templates[sName] = one_bin(infile_path, sName, region, 1, "mucr_", syst="nominal")
                    nominal = templates[sName][0]

                    if nominal < eps:
                        print(f"Sample {sName} is too small, skipping")
                        continue

                    stype = rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + "_" + sName, stype, templates[sName])

                    sample.setParamEffect(
                        sys_lumi_uncor, lumi_err[year[:4]] ** (LUMI[year[:4]] / LUMI["2022-2024"])
                    )
                    if do_systematics:

                        sample.autoMCStats(lnN=True)

                    ch.addSample(sample)

                data_obs = one_bin(infile_path, "Muondata", region, 1, "mucr_", syst="nominal")
                ch.setObservation(data_obs, read_sumw2=True)

            tqqpassbb = muonCR_model["muonCRpassbb" + year + "_ttbar"]
            tqqpasscc = muonCR_model["muonCRpasscc" + year + "_ttbar"]
            tqqfail = muonCR_model["muonCRfail" + year + "_ttbar"]

            sumPass = (
                tqqpassbb.getExpectation(nominal=True).sum()
                + tqqpasscc.getExpectation(nominal=True).sum()
            )
            sumPassbb = tqqpassbb.getExpectation(nominal=True).sum()
            sumPasscc = tqqpasscc.getExpectation(nominal=True).sum()
            sumFail = tqqfail.getExpectation(nominal=True).sum()

            stqqpassbb = muonCR_model["muonCRpassbb" + year + "_singlet"]
            stqqpasscc = muonCR_model["muonCRpasscc" + year + "_singlet"]
            stqqfail = muonCR_model["muonCRfail" + year + "_singlet"]

            sumPass += stqqpassbb.getExpectation(nominal=True).sum()
            sumPass += stqqpasscc.getExpectation(nominal=True).sum()

            sumPassbb += stqqpassbb.getExpectation(nominal=True).sum()
            sumPasscc += stqqpasscc.getExpectation(nominal=True).sum()
            sumFail += stqqfail.getExpectation(nominal=True).sum()

            tqqPF = sumPass / sumFail
            tqqBC = sumPassbb / sumPasscc

            tqqpassbb.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            tqqpasscc.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)

            tqqpassbb.setParamEffect(tqqeffBCSF, 1 * tqqeffBCSF)
            tqqpasscc.setParamEffect(tqqeffBCSF, (1 - tqqeffBCSF) * tqqBC + 1)

            tqqpassbb.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            tqqpasscc.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            tqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

            stqqpassbb.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            stqqpasscc.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            stqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)

            stqqpassbb.setParamEffect(tqqeffBCSF, 1 * tqqeffBCSF)
            stqqpasscc.setParamEffect(tqqeffBCSF, (1 - tqqeffBCSF) * tqqBC + 1)

            stqqpassbb.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            stqqpasscc.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            stqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

        # ---------------------------------------------------------
        # 7. SAVE & RENDER
        # ---------------------------------------------------------
        with (datacard_dir / f"{analysis}_{model_name}Model_{year}.pkl").open("wb") as fout:
            pickle.dump(model_dict[model_name], fout)
        modeldir = datacard_dir / f"{analysis}_{model_name}Model_{year}"
        if do_muon_CR:
            muonCR_model.renderCombine(modeldir)
        model_dict[model_name].renderCombine(modeldir)
        print(f"Datacards saved to {modeldir}")

        out_cards = ""
        for ch in model_dict[model_name]:
            if "/" in ch.name:
                continue
            out_cards += f"{ch.name}={ch.name}.txt "
            with Path(f"{modeldir}/{ch.name}.txt").open("a") as f:
                f.write("\nqcd_norm rateParam * qcd 1.0 [0,20]\n")
        if do_muon_CR:
            for ch in muonCR_model:
                if "/" in ch.name:
                    continue
                out_cards += f"{ch.name}={ch.name}.txt "

        model_cls = pm_config[model_name]["model"]
        maps = " ".join([f"--PO '{m}'" for m in pm_config[model_name]["maps"]])
        t2w_extra = pm_config[model_name].get("t2w_extra", "")
        bpre_in = pm_config[model_name].get("build_preamble", [])
        bpre_out = ""
        if len(bpre_in) > 0:
            bpre_out = "".join([f"{b}\n" for b in bpre_in])
            

        # Construct the text2workspace command dynamically
        t2w_cfg = f"-P {model_cls} --PO verbose {maps} {t2w_extra}"

        # Write the build script
        with (modeldir / "build.sh").open("w") as f:
            f.write("#!/bin/bash\n")
            f.write(bpre_out)
            f.write(f"combineCards.py {out_cards} > model_combined.txt\n")
            f.write(f"text2workspace.py {t2w_cfg} model_combined.txt -o workspace.root\n")
            f.write("echo 'Workspace created: workspace.root'\n")

        (modeldir / "build.sh").chmod(0o755)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--indir", default=None)
    parser.add_argument(
        "--outdir", default="results", help="Output directory for datacards and plots"
    )
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--debug", action="store_true", help="Enter debug mode")
    args = parser.parse_args()
    rhalphabet(args)
