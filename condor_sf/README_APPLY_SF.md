# Apply Scale Factors to MC Efficiencies

## Overview

The script `apply_scale_factors_2d.py` takes the calculated 2D scale factors and applies them to MC efficiencies to produce corrected efficiency predictions. This allows you to visualize how well the scale factor corrections work.

## Purpose

After calculating scale factors using `scale_factor_2d.py`, this script:
1. Loads the data and MC efficiency outputs
2. Recalculates efficiencies and scale factors
3. Applies scale factors to MC efficiencies (MC × SF)
4. Compares:
   - Data efficiency
   - MC efficiency (uncorrected)
   - MC efficiency × SF (corrected)
   - Scale factors themselves

## Usage

### Basic Usage

```bash
# Process both VBF and ggF for 2023
python apply_scale_factors_2d.py --year 2023

# Process only VBF
python apply_scale_factors_2d.py --year 2023 --prod_mode VBF

# Process only ggF
python apply_scale_factors_2d.py --year 2023 --prod_mode ggF

# Specify custom input directory
python apply_scale_factors_2d.py --year 2023 --input_dir /path/to/outputs
```

### Command Line Arguments

- `--year`: Data-taking period (choices: `2022`, `2022EE`, `2023`, `2023BPix`)
  - Default: `2023`

- `--prod_mode`: Production mode (choices: `VBF`, `ggF`, `both`)
  - Default: `both`
  - Processes both production modes sequentially

- `--input_dir`: Base directory containing input coffea files
  - Default: current working directory
  - Should contain `output/` subdirectory with `.coffea` files

## Input Requirements

The script expects these files to exist:
```
{input_dir}/output/sf_2d_{year}_{prod_mode}_data.coffea
{input_dir}/output/sf_2d_{year}_{prod_mode}_mc.coffea
```

These files are created by running `scale_factor_2d.py`.

## Output

### Saved Files

All outputs are saved to: `figures_scaled_sf_2d/{year}/{prod_mode}/`

#### 1. 2D Comparison Plot
**File**: `scaled_efficiency_comparison_pt_vs_msd_2d.png`

A 2×2 panel plot showing:
- **Top Left**: Data efficiency (2D heatmap)
- **Top Right**: MC efficiency - Uncorrected (2D heatmap)
- **Bottom Left**: MC efficiency × SF - Corrected (2D heatmap)
- **Bottom Right**: Scale factor map (2D heatmap)

Color scales:
- Efficiencies: RdYlGn (0-1)
- Scale factors: RdBu_r (0.8-1.2)

#### 2. 1D Projection Plots

**Files**:
- `scaled_eff_projection_x_{prod_mode}.png` - Projection onto mSD axis
- `scaled_eff_projection_y_{prod_mode}.png` - Projection onto pT axis

Each plot has:
- **Top panel**: Efficiency comparison
  - Black solid line: Data
  - Blue dashed line: MC (uncorrected)
  - Red dotted line: MC × SF (corrected)
- **Bottom panel**: Ratio to data
  - Blue dashed: MC / Data
  - Red dotted: (MC × SF) / Data
  - Should be close to 1.0 after correction

### Plot Features

- Standard CMS labels with year and COM energy
- Grid lines for easier reading
- Legend identifying each curve
- Overall statistics in titles

## Validation

This script is useful for validating that:
1. Scale factors are correctly calculated
2. MC × SF efficiencies match data better than uncorrected MC
3. No bins have extreme scale factor values
4. Coverage is adequate across the 2D phase space

### Good Validation Signs

- MC × SF curves track data closely
- Ratio plots show MC × SF / Data ≈ 1.0
- Scale factors are mostly in range [0.9, 1.1]
- Smooth variations (no large discontinuities)

### Warning Signs

- Large deviations of MC × SF from data
- Scale factors < 0.8 or > 1.2
- Large statistical fluctuations
- Empty bins in critical regions

## Example Workflow

```bash
# Step 1: Calculate scale factors (must be done first)
python scale_factor_2d.py --year 2023

# Step 2: Apply scale factors and create validation plots
python apply_scale_factors_2d.py --year 2023

# Step 3: Review plots in figures_scaled_sf_2d/2023/
```

## Technical Details

### Binning
Uses the same binning as `scale_factor_2d.py`:
- **Leading jet mass (mSD)**: 15 bins from 40 to 200 GeV
- **Leading jet pT**: 20 bins from 200 to 1000 GeV

### Calculations

For each bin (i, j):
- `eff_data[i,j]` = N_pass_data / N_total_data
- `eff_mc[i,j]` = N_pass_mc / N_total_mc
- `sf[i,j]` = eff_data[i,j] / eff_mc[i,j]
- `eff_mc_scaled[i,j]` = eff_mc[i,j] × sf[i,j]

Ideally: `eff_mc_scaled[i,j]` ≈ `eff_data[i,j]`

### 1D Projections

Projections are calculated by averaging over the other dimension:
- **X-projection** (mSD): Average over all pT bins
- **Y-projection** (pT): Average over all mSD bins

This shows overall trends while reducing statistical fluctuations.

## Dependencies

Required packages:
- `coffea`
- `matplotlib`
- `mplhep`
- `numpy`

## Notes

- This script does NOT require ROOT files, only the coffea outputs
- Can be run quickly after `scale_factor_2d.py` completes
- Useful for presentations and validation
- The corrected MC efficiency should closely match data if SFs are working correctly
