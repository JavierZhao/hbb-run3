# 2D Scale Factor Analysis

## Overview

The script `scale_factor_2d.py` calculates 2D trigger scale factors by comparing data and MC trigger efficiencies as a function of **leading jet mass** (mSD) and **leading jet pT**. It produces:
- Data efficiency maps (2D)
- MC efficiency maps (2D)
- Scale factor maps (Data/MC ratio in 2D)

## Features

- Processes both **VBF** and **ggF** production modes automatically
- Uses muon-triggered events as reference (tag-and-probe method)
- Applies baseline selection matching the production mode
- Outputs both data and MC results for comparison

## Usage

### Local Execution

```bash
# Basic usage
python scale_factor_2d.py --year 2023

# Test mode (10 files only)
python scale_factor_2d.py --year 2023 --test
```

### Condor Submission

#### Easy Method (Recommended)

```bash
# Submit for 2023
./submit_sf_2d.sh 2023

# Submit for other years
./submit_sf_2d.sh 2022
./submit_sf_2d.sh 2022EE
./submit_sf_2d.sh 2023BPix
```

#### Manual Method

```bash
# Create output directories
mkdir -p sf_2d_2023/{figures_sf_2d,output}

# Submit
condor_submit sf_2d-condor.jdl YEAR=2023
```

## Command Line Arguments

- `--year`: Data-taking period (choices: `2022`, `2022EE`, `2023`, `2023BPix`)
  - Default: `2022`

- `--test`: Test mode flag (no value needed)
  - When enabled, only processes the first 10 ROOT files
  - Useful for quick testing before running at full scale

## Output

### Saved Files

1. **Coffea output files**:
   - `output/sf_2d_{year}_{prod_mode}_data.coffea`
   - `output/sf_2d_{year}_{prod_mode}_mc.coffea`
   - Contains the processed histograms for data and MC

2. **2D SF plots**: `figures_sf_2d/{year}/{prod_mode}/`
   - Three-panel plots showing:
     - Left: Data efficiency
     - Middle: MC efficiency
     - Right: Scale factor (Data/MC)
   - Format: `scale_factor_pt_vs_msd_2d.png`

### Plot Features

Each 3-panel plot shows:
- **X-axis**: Leading jet mass (mSD) [40-200 GeV]
- **Y-axis**: Leading jet pT [200-1000 GeV]
- **Left panel**: Data efficiency [0-1, RdYlGn colormap]
- **Middle panel**: MC efficiency [0-1, RdYlGn colormap]
- **Right panel**: Scale factor [0.8-1.2, RdBu_r colormap]
- **Titles**: Include overall efficiency and SF values
- **CMS labels**: Standard CMS labels with year and COM energy

## Method

### Reference Triggers (Muon)

Used to select unbiased events for efficiency measurement:
- `HLT_Mu50`
- `HLT_CascadeMu100`
- `HLT_HighPtTkMu100`
- `HLT_IsoMu24`

### Signal Triggers by Year

#### 2022 / 2022EE
- `AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35`
- `QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65`

#### 2023 / 2023BPix
- `AK8PFJet250_SoftDropMass40_PNetBB0p06`
- `PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70`
- `VBF_DiPFJet125_45_Mjj720_Detajj3p0`

### Baseline Selection

The script applies different baseline selections depending on the production mode:

#### Common Selection (Both VBF and ggF)
1. At least one AK4 jet with pT > 15 GeV
2. At least one FatJet with:
   - pT > 250 GeV
   - |η| < 2.5
   - ParticleNet XbbVsQCD > 0.4
3. Muon requirements:
   - At least one loose muon (pT > 25 GeV, |η| < 2.4, loose ID, pfRelIso04 < 0.25)
   - ΔR(muon, leading FatJet) > 0.8

#### VBF-Specific Requirements
In addition to common selection:
- At least two AK4 jets with pT > 30 GeV and |η| < 5.0
- ΔR(AK4, leading FatJet) > 0.8
- |Δη(j1, j2)| > 3.5
- mjj > 1000 GeV

#### ggF Selection
Only the common selection (no additional requirements)

## Technical Details

### Binning
- **Leading jet mass (mSD)**: 15 bins from 40 to 200 GeV
- **Leading jet pT**: 20 bins from 200 to 1000 GeV

### Efficiency Calculation
For each (mSD, pT) bin:
- **Denominator**: Events passing baseline + reference muon trigger
- **Numerator**: Events passing baseline + reference muon trigger + signal HLT trigger
- **Efficiency**: Numerator / Denominator

### Scale Factor Calculation
For each (mSD, pT) bin:
- **SF**: Data Efficiency / MC Efficiency

### Data vs MC

The script processes:
- **Data**: Muon-triggered collision data
- **MC**: ttbar simulation

Both datasets are required for SF calculation.

## Input Files

The script expects JSON files in the following structure:
```
infiles/{year}/{year}_MuonData.json  # For data
infiles/{year}/{year}_ttbar.json     # For MC
```

Example:
```
infiles/2023/2023_MuonData.json
infiles/2023/2023_ttbar.json
```

Each JSON file should contain a dictionary mapping dataset names to lists of ROOT file paths.

## Dependencies

Required packages:
- `awkward`
- `coffea`
- `hist`
- `matplotlib`
- `mplhep`
- `numpy`
- `uproot`

## Workflow

1. Script runs twice (once for VBF, once for ggF)
2. For each production mode:
   - Loads both data and MC files
   - Applies reference muon triggers
   - Applies baseline selection
   - Calculates efficiencies in 2D bins
   - Computes scale factors
   - Generates plots
3. Outputs saved to separate directories per production mode

## Notes

- Scale factors are only calculated in bins where both data and MC have sufficient statistics
- Bins with zero denominator are masked (shown as NaN)
- The OR of all signal triggers is used for the numerator
- The script processes both production modes automatically in a single run
