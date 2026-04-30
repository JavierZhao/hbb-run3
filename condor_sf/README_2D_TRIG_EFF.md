# 2D Trigger Efficiency Script

## Overview

The script `trig_eff_2d.py` calculates trigger efficiencies as a function of **leading jet mass** (mSD) and **leading jet pT** and produces 2D efficiency plots for:
- Each individual trigger
- The OR of all triggers

## Usage

### Basic Usage

```bash
python trig_eff_2d.py --year 2022 --prod_mode VBF
```

### Command Line Arguments

- `--year`: Data-taking period (choices: `2022`, `2022EE`, `2023`, `2023BPix`)
  - Default: `2022`

- `--prod_mode`: Production mode for baseline selection (choices: `VBF`, `ggF`)
  - Default: `VBF`

- `--test`: Test mode flag (no value needed)
  - When enabled, only processes the first 10 ROOT files
  - Useful for quick testing before running at full scale

### Examples

```bash
# Test mode: Process only 10 files for quick validation
python trig_eff_2d.py --year 2022 --prod_mode VBF --test

# Process 2022 data with VBF selection (full scale)
python trig_eff_2d.py --year 2022 --prod_mode VBF

# Process 2023 data with ggF selection
python trig_eff_2d.py --year 2023 --prod_mode ggF

# Process 2023BPix data with VBF selection
python trig_eff_2d.py --year 2023BPix --prod_mode VBF
```

## Output

### Saved Files

1. **Coffea output file**: `output/trig_eff_2d_{year}_{prod_mode}.coffea`
   - Contains the processed histograms for later analysis

2. **2D efficiency plots**: `figures_2d/{year}/{prod_mode}/`
   - One plot per trigger showing efficiency as a function of (mSD, pT)
   - One plot for the OR of all triggers
   - Format: `{trigger_name}_pt_vs_msd_2d.png`

### Plot Features

Each 2D plot shows:
- **X-axis**: Leading jet mass (mSD) [40-200 GeV]
- **Y-axis**: Leading jet pT [200-1000 GeV]
- **Color scale**: Trigger efficiency [0-1]
  - Red: Low efficiency
  - Yellow: Medium efficiency
  - Green: High efficiency
- **Title**: Includes trigger name and overall efficiency percentage
- **CMS label**: Standard CMS preliminary label with year and COM energy

## Triggers by Year

### 2022 / 2022EE
- `AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35`
- `QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65`

### 2023 / 2023BPix
- `AK8PFJet250_SoftDropMass40_PNetBB0p06`
- `PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70`
- `VBF_DiPFJet125_45_Mjj720_Detajj3p0`

## Baseline Selection

The script applies different baseline selections depending on the production mode:

### VBF Selection
1. At least one AK4 jet with pT > 15 GeV
2. At least one FatJet with:
   - pT > 250 GeV
   - |η| < 2.5
   - ParticleNet XbbVsQCD > 0.4
3. VBF-specific requirements:
   - At least two AK4 jets with pT > 30 GeV and |η| < 5.0
   - ΔR(AK4, leading FatJet) > 0.8
   - |Δη(j1, j2)| > 3.5
   - mjj > 1000 GeV

### ggF Selection
1. At least one AK4 jet with pT > 15 GeV
2. At least one FatJet with:
   - pT > 250 GeV
   - |η| < 2.5
   - ParticleNet XbbVsQCD > 0.4

## Technical Details

### Binning
- **Leading jet mass (mSD)**: 15 bins from 40 to 200 GeV
- **Leading jet pT**: 20 bins from 200 to 1000 GeV

### Efficiency Calculation
For each bin:
- **Numerator**: Events passing baseline + trigger
- **Denominator**: Events passing baseline only
- **Efficiency**: Numerator / Denominator

### OR Trigger
The OR trigger combines all individual triggers with a logical OR operation, representing the overall trigger efficiency when using any of the configured triggers.

## Dependencies

Required packages:
- `awkward`
- `coffea`
- `hist`
- `matplotlib`
- `mplhep`
- `numpy`
- `uproot`

## Input Files

The script expects JSON files in the following structure:
```
infiles/{year}/{year}_{prod_mode}.json
```

Example: `infiles/2022/2022_VBF.json`

Each JSON file should contain a dictionary mapping dataset names to lists of ROOT file paths.
