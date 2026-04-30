# Trigger Efficiency Analysis - Output Guide

This document describes all outputs produced by the trigger efficiency analysis pipeline and how they correspond to the Analysis Note (AN) requirements.

## Quick Start

To generate all plots and tables for the AN:

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf
./run_all_analyses.sh
```

For a quick test run (10 files only):
```bash
./run_all_analyses.sh --test
```

## Analysis Scripts

### Main Analysis Scripts

1. **scale_factor_2d.py** - 2D Scale Factor Analysis
   - Produces: 2D efficiency maps (pT vs mSD) for Data and MC
   - Produces: 2D scale factor maps (Data/MC ratio)
   - Produces: Per-trigger efficiency tables (CSV format)
   - Kinematic cuts: pT ≥ 450 GeV, 0 ≤ mSD < 300 GeV (for efficiency calculation)
   - Plotting range: pT ≥ 300 GeV, 0 ≤ mSD < 300 GeV

2. **scale_factor.py** - 1D Scale Factor Analysis
   - Produces: 1D efficiency turn-on curves (vs pT, mSD, HT, etc.)
   - Produces: 1D scale factor curves (Data/MC ratio)
   - Variables: pT, mSD, HT, num_ak4, particleNet_XbbVsQCD

3. **trig_eff_2d.py** - 2D Trigger Efficiency (MC or Data only)
   - Produces: 2D trigger efficiency maps
   - Used for: Validation and cross-checks

4. **trig_eff.py** - 1D Trigger Efficiency (MC or Data only)
   - Produces: 1D trigger efficiency curves
   - Used for: Validation and cross-checks

### Wrapper Scripts

- `run_sf_2d.sh` - Runs 2D scale factor analysis for a single period
- `run_sf.sh` - Runs 1D scale factor analysis for a single period
- `run_trig_eff_2d.sh` - Runs 2D trigger efficiency for a single period
- `run_all_analyses.sh` - **Master script** to run all analyses for all periods

## Output Structure

### 2D Scale Factor Plots

**Location:** `figures_sf_2d/{YEAR}/{MODE}/`

**Files:**
- `scale_factor_pt_vs_msd_2d.png` - Three-panel plot showing:
  - Left: Data efficiency map
  - Middle: MC efficiency map
  - Right: Scale factor map (Data/MC)

**Periods:** 2022, 2022EE, 2023, 2023BPix
**Modes:** VBF, ggF

**Example:**
```
figures_sf_2d/2022/VBF/scale_factor_pt_vs_msd_2d.png
figures_sf_2d/2022/ggF/scale_factor_pt_vs_msd_2d.png
figures_sf_2d/2022EE/VBF/scale_factor_pt_vs_msd_2d.png
...
```

### Per-Trigger Efficiency Tables (CSV)

**Location:** `output/{YEAR}/{MODE}/`

**Files:**
- `per_trigger_efficiency_{YEAR}_{MODE}.csv`

**Format:**
| Trigger | Data_Efficiency | Data_Pass | Data_Total | MC_Efficiency | MC_Pass | MC_Total | Scale_Factor |
|---------|----------------|-----------|------------|---------------|---------|----------|--------------|
| AK8PFJet250_... | 0.8542 | 12345 | 14456 | 0.8923 | 9876 | 11067 | 0.9573 |
| QuadPFJet70_... | 0.7234 | 10456 | 14456 | 0.7689 | 8512 | 11067 | 0.9408 |
| Combined (OR) | 0.9123 | 13189 | 14456 | 0.9456 | 10465 | 11067 | 0.9648 |

**Kinematic Region:** 450 ≤ pT < 1000 GeV, 0 ≤ mSD < 300 GeV (as per Run 2 AN)

**Usage in AN:** These tables correspond to Tables 2, 3, 4 in the Run 2 Trigger AN showing individual trigger efficiencies for each data-taking period.

### 1D Scale Factor Plots

**Location:** `sf_{YEAR}_OR/`

**Files:**
- `scale_factor_pt_{YEAR}.png` - Leading jet pT turn-on curve
- `scale_factor_msd_{YEAR}.png` - Leading jet mSD distribution
- `scale_factor_ht_{YEAR}.png` - HT distribution
- `scale_factor_num_ak4_{YEAR}.png` - Number of AK4 jets
- `scale_factor_particleNet_XbbVsQCD_{YEAR}.png` - ParticleNet score

Each plot shows:
- Data efficiency (black)
- MC efficiency (red)
- Scale factor (ratio panel below)

### Coffea Output Files

**Location:** `output/`

**Files:**
- `sf_2d_{YEAR}_{MODE}_data.coffea` - Data histograms for 2D analysis
- `sf_2d_{YEAR}_{MODE}_mc.coffea` - MC histograms for 2D analysis

These can be reloaded for further analysis or replotting without reprocessing.

## Analysis Details

### Trigger Lists by Period

**2022, 2022EE:**
- `HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35`
- `HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65`

**2023, 2023BPix:**
- `HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06`
- `HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70`
- `HLT_VBF_DiPFJet125_45_Mjj720_Detajj3p0`

### Muon Reference Triggers (All periods)
- `HLT_Mu50`
- `HLT_CascadeMu100`
- `HLT_HighPtTkMu100`
- `HLT_IsoMu24`

### Selection Criteria

**Baseline Selection:**
- Leading FatJet: pT > 250 GeV, |η| < 2.5, ParticleNet_XbbVsQCD > 0.4
- Muon tag: pT > 25 GeV, |η| < 2.4, looseId, RelIso < 0.25, ΔR(μ, FatJet) > 0.8

**VBF Selection (additional):**
- Two AK4 jets outside FatJet (ΔR > 0.8)
- Δη(j1, j2) > 3.5
- mjj > 1000 GeV

### Kinematic Regions

**For Efficiency Calculation (CSV tables):**
- 450 ≤ pT < 1000 GeV
- 0 ≤ mSD < 300 GeV

**For Plotting (all figures):**
- 300 ≤ pT < 1000 GeV
- 0 ≤ mSD < 300 GeV

This follows the Run 2 AN methodology where a stricter pT cut is applied for the quoted efficiencies, while plots show the full turn-on curve.

## For the Analysis Note

### Recommended Plots and Tables

1. **Section: Trigger Efficiency Measurement**
   - Include: 2D efficiency maps for one representative period (e.g., 2023)
   - Use: `figures_sf_2d/2023/VBF/scale_factor_pt_vs_msd_2d.png`

2. **Section: Per-Trigger Efficiencies**
   - Include: Per-trigger efficiency tables (similar to Run 2 Tables 2-4)
   - Use: `output/2022/VBF/per_trigger_efficiency_2022_VBF.csv` (convert to LaTeX table)
   - Use: `output/2022EE/VBF/per_trigger_efficiency_2022EE_VBF.csv`
   - Use: `output/2023/VBF/per_trigger_efficiency_2023_VBF.csv`
   - Use: `output/2023BPix/VBF/per_trigger_efficiency_2023BPix_VBF.csv`

3. **Section: 1D Efficiency Turn-on Curves**
   - Include: pT turn-on curves for all periods
   - Use: `sf_{YEAR}_OR/scale_factor_pt_{YEAR}.png`

4. **Section: Scale Factors**
   - Include: Summary of scale factors from all periods
   - Extract from CSV tables: Combined (OR) scale factor values

### Converting CSV to LaTeX

Example Python snippet:
```python
import pandas as pd

df = pd.read_csv('output/2022/VBF/per_trigger_efficiency_2022_VBF.csv')
latex = df.to_latex(index=False, float_format="%.4f")
print(latex)
```

## Running Individual Analyses

### Run single period, single analysis:
```bash
# 2D scale factors for 2022
./run_sf_2d.sh --year 2022

# 1D scale factors for 2023
./run_sf.sh --year 2023

# Test mode (10 files only)
./run_sf_2d.sh --year 2022 --test
```

### Run all periods, all analyses:
```bash
./run_all_analyses.sh
```

## Troubleshooting

### If jobs fail:
1. Check log files in `sf_2d_{YEAR}/` or `sf_{YEAR}_OR/`
2. Run in test mode first: `./run_all_analyses.sh --test`
3. Check that input JSON files exist: `ls infiles/{YEAR}/*.json`

### If plots are missing:
1. Verify output directories were created: `ls -d figures_sf_2d output`
2. Check for Python errors in log files
3. Verify matplotlib backend: `export MPLBACKEND=Agg`

## Notes

- All analyses use TTtoLNu2Q (semileptonic ttbar) samples for MC
- The muon tag-and-probe method is used to measure trigger efficiencies
- Scale factors are computed as Data efficiency / MC efficiency
- Both VBF and ggF production modes are processed separately
- Per-trigger efficiencies are calculated independently before combining with OR

## Contact

For questions about this analysis pipeline, refer to:
- TRG_AN_RUN_2.pdf (Run 2 trigger efficiency AN methodology)
- HIG-21-020 (Run 2 Hbb analysis)
