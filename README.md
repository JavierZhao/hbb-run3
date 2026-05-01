# hbb-run3

CMS Run 3 H→bb trigger efficiency and scale-factor analysis. Aim is to reproduce the plot and
table set from the Run 2 trigger AN ([`docs/TRG_AN_RUN_2.pdf`](docs/TRG_AN_RUN_2.pdf)) for Run 3
data (2022, 2022EE, 2023, 2023BPix), using a muon-triggered tag-and-probe method.

The active code lives in [`condor_sf/`](condor_sf/). Three pipelines:

| Stage                             | Script                              | Inputs                         | Outputs                              |
|-----------------------------------|-------------------------------------|--------------------------------|--------------------------------------|
| (a) Signal MC trigger efficiency  | `trig_eff_2d.py`                    | VBF / ggF MC                   | `trig_eff_2d_<YEAR>/`                |
| (b) Data/MC scale factor          | `scale_factor_2d.py`                | Muon-triggered data + ttbar MC | `sf_2d_<YEAR>/`                      |
| (c) ParkingVBF efficiency study   | `parking_vbf_eff.py`                | ParkingVBF data                | `parking_vbf_<YEAR>/`                |

Closure validation: `apply_scale_factors_2d.py`.


## What is and isn't reproduced from the Run 2 AN

| Run 2 AN element                                         | Run 3 status | Where                                                                     |
|----------------------------------------------------------|:---:|---------------------------------------------------------------------------------------|
| **Tables 2-4** — per-trigger efficiency (450 ≤ pT < 1000 GeV, 40 ≤ mSD < 300 GeV) | ✅ | Two CSVs per period: `per_trigger_efficiency_<year>_*.csv` (period-integrated, with σ columns) and `per_era_efficiency_<year>_*.csv` (AN-style: rows = triggers, columns = data eras + MC). Companion `*_err.csv` carries the per-cell σ. |
| **Figs 3, 6, 10** — 1D efficiency vs leading-jet pT, data + simulation panels  | ✅ | `figures_sf_2d/<year>/Inclusive/<region>/efficiency_vs_pt.png` |
| **Figs 4, 7, 11** — 1D efficiency vs leading-jet mSD (pT ≥ 450 GeV plateau)    | ✅ | `figures_sf_2d/<year>/Inclusive/<region>/efficiency_vs_msd.png` |
| **Figs 5, 8, 12 (left)** — 2D scale-factor map (pT vs mSD)                     | ✅ | `figures_sf_2d/<year>/Inclusive/<region>/scale_factor_pt_vs_msd_2d.png` |
| **Figs 5, 8, 12 (right)** — statistical-uncertainty map for the 2D SF          | ✅ | `figures_sf_2d/<year>/Inclusive/<region>/scale_factor_pt_vs_msd_2d_unc.png`; per-bin `σ(SF) = SF·√((1−ε_d)/k_d + (1−ε_m)/k_m)` (binomial Gaussian, not Clopper-Pearson). Per-bin σ also appears as error bars on the 1D projections and as `*_Err` columns in the per-trigger CSV. |
| **Fig 9** — efficiency vs leading-jet (X<sub>bb</sub>+X<sub>cc</sub>) tagger score | ✅ | `figures_sf_2d/<year>/Inclusive/<region>/efficiency_vs_score.png`; data + MC two-panel, with the kinematic plateau cut applied. |

Plot kinematic axes use the **full range** (`pT 200–1000 GeV`, `mSD 0–300 GeV`) so the
turn-on stays visible. Table efficiency numbers are still computed in the AN window
`450 ≤ pT < 1000 GeV, 40 ≤ mSD < 300 GeV`. Both windows are configurable from the CLI.

Beyond the Run 2 AN, the SF measurement is also produced split by candidate-jet
Xbb+Xcc working point (pass/fail at WP = 0.82) — see `--txbb-region` below.


## How to run

### Option 1 — full production (HTCondor)

From inside `condor_sf/`:

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf

# Per-period jobs. Submit each year separately, or all in one batch.
condor_submit sf_2d-condor.jdl YEAR=2022
condor_submit sf_2d-condor.jdl YEAR=2022EE
condor_submit sf_2d-condor.jdl YEAR=2023
condor_submit sf_2d-condor.jdl YEAR=2023BPix

# Equivalent for signal-MC trigger efficiency
condor_submit trig_eff_2d-condor.jdl YEAR=2023

# ParkingVBF (2023 + 2023BPix queued from one submit)
condor_submit parking_vbf-condor.jdl
```

Outputs land in `condor_sf/sf_2d_<YEAR>/{output, figures_sf_2d}/` (and equivalent for the
other pipelines). Monitor with `condor_q` and `tail -f sf_2d_<YEAR>/*.out`.

The JDLs use the `coffeateam/coffea-dask-almalinux8:2025.10.2-py3.10` apptainer image.
After the job finishes, run the closure step locally:

```bash
python apply_scale_factors_2d.py --year 2023
```

### Option 2 — quick interactive test

Runs the full pipeline on a small subset of files inside the same container the condor
jobs use. Output goes to `condor_sf/my_test_run/` so it doesn't collide with production
runs:

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf
./shell coffeateam/coffea-dask-almalinux8:2025.10.2-py3.10 \
  "python /srv/scale_factor_2d.py \
    --year 2023 \
    --quick-test \
    --quick-files-per-dataset 80 \
    --quick-maxchunks 400 \
    --quick-workers 1 \
    --txbb-region pass \
    --outdir /srv/my_test_run"
```

Useful quick-test knobs (defined in `scale_factor_2d.py`):

| Flag | Default | Meaning |
|---|---|---|
| `--year` | `2022` | Period: `2022`, `2022EE`, `2023`, `2023BPix` |
| `--quick-test` | off | Enable the small-scale debug mode |
| `--quick-files-per-dataset N` | 2 | Files per dataset to read |
| `--quick-maxchunks N` | 2 | Max coffea chunks per dataset |
| `--quick-workers N` | 1 | Futures-executor workers |
| `--txbb-region` | `all` | `pass` (≥0.82), `fail` (<0.82), `inclusive`, or `all` (run pass+fail) |
| `--pt-min/max`, `--msd-min/max` | 450/1000, 40/300 | Kinematic window for the per-trigger CSV |
| `--outdir DIR` | `.` | Base directory for `figures_sf_2d/` and `output/` |
| `--disable-jetid-correction` | off | Use NanoAOD `isTight` instead of correctionlib JetID |


## Repository layout

```
hbb-run3/
├── condor_sf/              ← active 2D pipeline (main working dir)
│   ├── scale_factor_2d.py        scale-factor processor + 1D / 2D plotters + per-trigger CSV
│   ├── trig_eff_2d.py            signal-MC efficiency processor + plotter
│   ├── parking_vbf_eff.py        ParkingVBF dijet efficiency
│   ├── apply_scale_factors_2d.py closure validation
│   ├── *-condor.jdl              HTCondor JDLs (parameterized via $(YEAR))
│   ├── run_*.sh                  in-container runners
│   ├── shell                     apptainer/singularity wrapper
│   ├── infiles/<year>/*.json     dataset → file lists
│   └── README.md (+ topic-specific README_*.md)
├── boostedhiggs/           ← shared processor / corrections (Connor's framework, lightly modified)
└── trig_eff/               ← legacy notebook-driven 1D workflow, kept for reference
```


## References

- Run 2 trigger AN — [`docs/TRG_AN_RUN_2.pdf`](docs/TRG_AN_RUN_2.pdf) — target plot/table set
- [Boosted H(bb) Run 3 status presentation](https://indico.cern.ch/event/1624978/contributions/6875328/) — txbb WP = 0.82 motivation
- [DAZSLE/hbb-run3 corrections](https://github.com/DAZSLE/hbb-run3/blob/main/src/hbb/corrections.py) — JetID correction reference for v14+

Originally based on [@jet96`s](https://github.com/jet96) datasets and Connor's `boostedhiggs` framework.
