# Trigger Efficiency and Scale Factor Analysis Suite

Complete workflow for calculating 2D trigger efficiencies and scale factors for HH→bb analysis.

## Overview

This suite provides three main analysis scripts:

1. **`trig_eff_2d.py`** - Calculate 2D trigger efficiencies
2. **`scale_factor_2d.py`** - Calculate 2D scale factors (Data/MC)
3. **`apply_scale_factors_2d.py`** - Apply and validate scale factors

All scripts support both **VBF** and **ggF** production modes.

## Quick Start

### Local Execution

```bash
# 1. Calculate trigger efficiencies (VBF + ggF)
cd trig_eff
python trig_eff_2d.py --year 2023

# 2. Calculate scale factors (VBF + ggF)
python scale_factor_2d.py --year 2023

# 3. Apply and validate scale factors
python apply_scale_factors_2d.py --year 2023
```

### Condor Submission

```bash
# 1. Submit trigger efficiency jobs
./condor_sf/submit_trig_eff_2d.sh 2023

# 2. Submit scale factor jobs (after step 1 completes)
./condor_sf/submit_sf_2d.sh 2023

# 3. Apply scale factors (run locally after step 2 completes)
cd trig_eff
python apply_scale_factors_2d.py --year 2023
```

## Analysis Scripts

### 1. Trigger Efficiency (`trig_eff_2d.py`)

**Purpose**: Calculate trigger efficiencies in 2D (pT vs mSD) for signal triggers.

**Output**:
- `output/trig_eff_2d_{year}_{prod_mode}.coffea`
- `figures_2d/{year}/{prod_mode}/*.png` - 2D efficiency maps per trigger

**Details**: See [README_2D_TRIG_EFF.md](README_2D_TRIG_EFF.md)

### 2. Scale Factors (`scale_factor_2d.py`)

**Purpose**: Calculate scale factors by comparing data and MC trigger efficiencies.

**Output**:
- `output/sf_2d_{year}_{prod_mode}_data.coffea`
- `output/sf_2d_{year}_{prod_mode}_mc.coffea`
- `figures_sf_2d/{year}/{prod_mode}/scale_factor_*.png` - SF maps

**Details**: See [README_SF_2D.md](README_SF_2D.md)

### 3. Apply Scale Factors (`apply_scale_factors_2d.py`)

**Purpose**: Apply SFs to MC and validate correction quality.

**Output**:
- `figures_scaled_sf_2d/{year}/{prod_mode}/*.png` - Comparison plots

**Details**: See [README_APPLY_SF.md](README_APPLY_SF.md)

## Complete Workflow

### Standard Analysis Flow

```
Step 1: Calculate Efficiencies
├─ Run: trig_eff_2d.py
├─ Input: Signal MC (VBF/ggF)
└─ Output: Efficiency maps

Step 2: Calculate Scale Factors
├─ Run: scale_factor_2d.py
├─ Input: Data (muon-triggered) + MC (ttbar)
└─ Output: Scale factor maps

Step 3: Validate Scale Factors
├─ Run: apply_scale_factors_2d.py
├─ Input: Outputs from Step 2
└─ Output: Validation plots
```

### Dependency Chain

```
Input Data Files
    ├─ infiles/{year}/{year}_VBF.json
    ├─ infiles/{year}/{year}_ggF.json
    ├─ infiles/{year}/{year}_MuonData.json
    └─ infiles/{year}/{year}_ttbar.json
           ↓
    trig_eff_2d.py (optional, for signal efficiency)
           ↓
    scale_factor_2d.py (required for SFs)
           ↓
    apply_scale_factors_2d.py (validation)
```

## Directory Structure

```
condor_sf/
├── README.md                        # This file
├── README_2D_TRIG_EFF.md           # Trigger efficiency docs
├── README_SF_2D.md                 # Scale factor docs
├── README_APPLY_SF.md              # Apply SF docs
│
├── trig_eff_2d.py                  # Script 1: Efficiencies
├── trig_eff_2d-condor.jdl          # Condor JDL for script 1
├── run_trig_eff_2d.sh              # Wrapper for script 1
├── submit_trig_eff_2d.sh           # Submit script 1
│
├── scale_factor_2d.py              # Script 2: Scale factors
├── sf_2d-condor.jdl                # Condor JDL for script 2
├── run_sf_2d.sh                    # Wrapper for script 2
├── submit_sf_2d.sh                 # Submit script 2
│
├── apply_scale_factors_2d.py       # Script 3: Apply & validate
│
└── infiles/                        # Input JSON files
    ├── 2022/
    ├── 2022EE/
    ├── 2023/
    └── 2023BPix/
```

## Output Organization

After running all scripts:

```
Working Directory/
├── output/
│   ├── trig_eff_2d_{year}_{prod_mode}.coffea
│   ├── sf_2d_{year}_{prod_mode}_data.coffea
│   └── sf_2d_{year}_{prod_mode}_mc.coffea
│
├── figures_2d/                     # From trig_eff_2d.py
│   └── {year}/{prod_mode}/
│       └── {trigger}_pt_vs_msd_2d.png
│
├── figures_sf_2d/                  # From scale_factor_2d.py
│   └── {year}/{prod_mode}/
│       └── scale_factor_pt_vs_msd_2d.png
│
└── figures_scaled_sf_2d/           # From apply_scale_factors_2d.py
    └── {year}/{prod_mode}/
        ├── scaled_efficiency_comparison_pt_vs_msd_2d.png
        ├── scaled_eff_projection_x_{prod_mode}.png
        └── scaled_eff_projection_y_{prod_mode}.png
```

## Supported Years

- `2022` - Run 3 2022 C-D
- `2022EE` - Run 3 2022 E-G
- `2023` - Run 3 2023 C
- `2023BPix` - Run 3 2023 D

## Common Command Line Options

All scripts support:
- `--year {2022,2022EE,2023,2023BPix}` - Data-taking period
- `--test` - Test mode (first 10 files only)

Additional options:
- `trig_eff_2d.py`: (none - automatically processes VBF + ggF)
- `scale_factor_2d.py`: (none - automatically processes VBF + ggF)
- `apply_scale_factors_2d.py`: `--prod_mode {VBF,ggF,both}`, `--input_dir`

## Condor Submission

### Files Required for Condor

- JDL file (`*.jdl`) - Job description
- Shell wrapper (`run_*.sh`) - Environment setup
- Python script (`*.py`) - Analysis code
- Submission script (`submit_*.sh`) - Creates directories and submits

### Condor Submission Process

1. Submission script creates output directories
2. Condor transfers input files and scripts
3. Job runs in Docker container (coffeateam/coffea-dask:latest-py3.9)
4. Output files transferred back on completion

### Monitoring Jobs

```bash
# Check job status
condor_q

# Check detailed status
condor_q -better-analyze <job_id>

# Check logs
tail -f {output_dir}/*.out
tail -f {output_dir}/*.err
```

## Troubleshooting

### Common Issues

1. **Missing input files**
   - Check that JSON files exist in `infiles/{year}/`
   - Verify file paths in JSON are accessible

2. **Empty output**
   - Check that baseline selection isn't too restrictive
   - Verify trigger names match NanoAOD branches

3. **Memory errors on Condor**
   - Increase `request_memory` in JDL file
   - Default is 8 GB, can increase to 16 GB if needed

4. **Import errors (local)**
   - Ensure coffea environment is activated
   - Install missing packages: `pip install coffea matplotlib mplhep`

### Debug Mode

Use `--test` flag to process only 10 files:
```bash
python trig_eff_2d.py --year 2023 --test
python scale_factor_2d.py --year 2023 --test
```

## Analysis Tips

### Choosing Production Modes

- **VBF**: For VBF signal analysis (applies VBF topology cuts)
- **ggF**: For ggF signal analysis (no VBF topology cuts)
- **Both**: Run both modes for complete coverage

### Interpreting Results

#### Good Scale Factors
- Close to 1.0 (ideally 0.95-1.05)
- Smooth variations across phase space
- MC × SF matches data well

#### Warning Signs
- SF < 0.8 or > 1.2
- Large bin-to-bin fluctuations
- Poor statistics (empty bins)
- MC × SF still doesn't match data

### Best Practices

1. Always run test mode first
2. Check logs for warnings/errors
3. Verify input file statistics
4. Review plots before using SFs
5. Document any modifications to selections

## Performance

### Typical Runtime (Full Dataset)

- `trig_eff_2d.py`: ~30-60 minutes per prod_mode
- `scale_factor_2d.py`: ~30-60 minutes per prod_mode
- `apply_scale_factors_2d.py`: ~1 minute (post-processing only)

### Memory Usage

- Typical: 4-8 GB
- Peak: Up to 12 GB for large datasets

## Updates and Modifications

### Adding New Triggers

Edit trigger dictionaries in each script:
```python
trigger_dict_periods = {
    'YEAR': [
        'HLT_TriggerName1',
        'HLT_TriggerName2',
    ]
}
```

### Changing Binning

Edit in respective scripts:
```python
'axis_x': hist.axis.Regular(bins=N, start=MIN, stop=MAX, ...)
'axis_y': hist.axis.Regular(bins=N, start=MIN, stop=MAX, ...)
```

### Modifying Selections

Look for `create_baseline_selection_mask_new()` function and adjust cuts.

## References

- [CMS Trigger Documentation](https://twiki.cern.ch/twiki/bin/view/CMS/Triggers)
- [Coffea Framework](https://coffeateam.github.io/coffea/)
- [Scale Factor Method](https://twiki.cern.ch/twiki/bin/view/CMS/TagAndProbe)

## Support

For issues or questions:
1. Check the individual README files for each script
2. Review logs in output directories
3. Verify input file accessibility
4. Check Condor job status

## Version History

- v1.0 - Initial release with 2D analysis suite
- Supports Run 3 data (2022-2023)
- Automated VBF + ggF processing
