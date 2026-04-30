#!/usr/bin/env python3
"""
Run cutflow for 2022 VBF and ggF using 10 ROOT files each.

Usage (from condor_sf/):
    python run_cutflow_2022.py

Produces:
    - figures_cutflow/2022/VBF/cutflow_VBF.png
    - figures_cutflow/2022/ggF/cutflow_ggF.png
    - printed cutflow tables to stdout
"""
import sys
import os
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from trig_eff_2d import (
    trig_vars_2d,
    trigger_dict_periods,
    TriggerEfficiency2DProcessor,
    plot_cutflow,
)
from coffea import processor
from coffea.nanoevents import NanoAODSchema

YEAR    = "2022"
N_FILES = 10
TRIGGERS = trigger_dict_periods[YEAR]


def main():
    for prod_mode in ["VBF", "ggF"]:
        print(f"\n{'=' * 60}")
        print(f"  {prod_mode} cutflow  ({N_FILES} files, {YEAR})")
        print(f"{'=' * 60}\n")

        json_path = os.path.join("infiles", YEAR, f"{YEAR}_{prod_mode}.json")
        with open(json_path) as f:
            data = json.load(f)

        samples = list(data.values())[0][:N_FILES]
        fileset = {f"{prod_mode}_{YEAR}": samples}
        print(f"Loaded {len(samples)} files from {json_path}")

        runner = processor.Runner(
            executor=processor.FuturesExecutor(compression=None, workers=2),
            schema=NanoAODSchema,
            skipbadfiles=True,
        )

        out = runner(
            fileset,
            treename="Events",
            processor_instance=TriggerEfficiency2DProcessor(
                TRIGGERS, trig_vars_2d, baseline_key=prod_mode
            ),
        )
        output = out[0]

        # Print cutflow table
        print(f"\n{'Step':<55} {'Events':>10} {'Fraction':>10}")
        print("-" * 78)
        total = output["cutflow"]["All events"]
        for step, count in output["cutflow"].items():
            frac = count / total if total > 0 else 0
            print(f"{step:<55} {count:>10,} {frac:>9.1%}")

        # Save plot
        save_dir = os.path.join("figures_cutflow", YEAR, prod_mode)
        plot_cutflow(output["cutflow"], tag=prod_mode, save_dir=save_dir, show_img=False)


if __name__ == "__main__":
    main()
