# Analysis Structure: condor_sf vs trig_eff

## Two Different Analyses

### 1. **condor_sf/** - Trigger Scale Factor Measurement (Tag-and-Probe)

**Purpose:** Measure data/MC scale factors for trigger efficiency using muon tag-and-probe method

**Sample:**
- Data: Single-muon data
- MC: TTtoLNu2Q (semileptonic ttbar)

**Method:**
- Muon reference triggers select unbiased sample
- Measure efficiency of HBB triggers
- Calculate SF = ε_data / ε_MC

**Production Mode Separation:** **NO - Inclusive measurement!**
```python
prod_modes = ['Inclusive']  # ✓ Correct
```

**Why Inclusive?**
- Triggers are production-mode agnostic (only depend on jet kinematics)
- ttbar events don't separate into VBF/ggF naturally
- VBF cuts reject ttbar → poor statistics
- **Same SF map applies to both VBF and ggF** based on (pT, mSD)

**Output:**
```
figures_sf_2d/{YEAR}/Inclusive/scale_factor_pt_vs_msd_2d.png
output/{YEAR}/Inclusive/per_trigger_efficiency_{YEAR}_Inclusive.csv
```

**Files:**
- `scale_factor_2d.py` - 2D SF measurement (pT vs mSD)
- `scale_factor.py` - 1D SF measurement (pT, mSD, HT, etc.)

---

### 2. **trig_eff/** - Signal MC Trigger Efficiency Evaluation

**Purpose:** Evaluate trigger efficiency on signal MC samples (for validation/studies)

**Sample:**
- VBF H→bb MC (VBF production)
- ggF H→bb MC (gluon fusion production)

**Method:**
- Direct efficiency: ε = N(pass trigger) / N(total)
- No data comparison
- No scale factors

**Production Mode Separation:** **YES - Keep VBF/ggF separate!**
```python
prod_modes = ['VBF', 'ggF']  # ✓ Correct
```

**Why Separate?**
- Using **different MC samples** (VBF MC vs ggF MC)
- Want to study efficiency differences between production modes
- VBF events have different kinematics than ggF
- Useful for understanding trigger performance on signal

**Output:**
```
figures_trig_eff/{YEAR}/VBF/trig_eff_pt_vs_msd_2d.png
figures_trig_eff/{YEAR}/ggF/trig_eff_pt_vs_msd_2d.png
```

**Files:**
- `scale_factor_2d.py` - 2D trigger efficiency (pT vs mSD)
- `trig_eff.py` - 1D trigger efficiency (pT, mSD, HT, etc.)

---

## Summary Table

| Aspect | condor_sf/ (Tag-and-Probe SF) | trig_eff/ (Signal MC Efficiency) |
|--------|-------------------------------|----------------------------------|
| **Purpose** | Measure trigger SF (Data/MC) | Evaluate trigger eff on signal MC |
| **Sample** | ttbar + Single-muon data | VBF H→bb MC, ggF H→bb MC |
| **Method** | Tag-and-probe | Direct efficiency |
| **VBF/ggF** | **Inclusive** (no separation) | **Separate** (different MC files) |
| **prod_modes** | `['Inclusive']` | `['VBF', 'ggF']` |
| **VBF cuts** | ❌ NO (would reject ttbar) | ✅ YES (define VBF phase space) |
| **Output** | SF maps for analysis | Efficiency plots for validation |

---

## How to Use

### For Trigger SF Measurement (condor_sf/):

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf

# Run tag-and-probe SF measurement (inclusive, no VBF/ggF split)
./submit_all_analyses.sh
```

**Expected output:**
- One SF map per year: `figures_sf_2d/2022/Inclusive/`
- High statistics (no VBF cuts rejecting ttbar)
- SF values close to 1.0

**Use in analysis:**
```python
# Apply same SF to both VBF and ggF events
sf = sf_map.GetBinContent(sf_map.FindBin(jet_msd, jet_pt))
event_weight *= sf
```

### For Signal MC Trigger Efficiency (trig_eff/):

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/trig_eff

# Run trigger efficiency on VBF and ggF signal MC
./run_trig_eff_2d.sh --year 2022
```

**Expected output:**
- Separate efficiency plots for VBF and ggF
- Shows if trigger performance differs by production mode
- Useful for validation and systematic studies

---

## Key Changes Made

### ✅ Updated: condor_sf/scale_factor_2d.py

**Before:**
```python
prod_modes = ['VBF', 'ggF']  # ❌ Wrong - created separate SFs
# Applied VBF topology cuts
```

**After:**
```python
prod_modes = ['Inclusive']  # ✅ Correct - one SF for all
# NO VBF topology cuts in trigger SF measurement
```

**Reason:** Triggers are production-mode agnostic. VBF cuts reject ttbar (poor statistics).

### ✅ Kept Unchanged: trig_eff/scale_factor_2d.py

```python
prod_modes = ['VBF', 'ggF']  # ✅ Correct - different MC samples
# VBF topology cuts ARE used to define VBF phase space
```

**Reason:** Evaluating efficiency on separate signal MC samples (VBF MC vs ggF MC).

---

## References

See `condor_sf/TRIGGER_SF_METHODOLOGY.md` for detailed methodology and supporting evidence from published CMS analyses.
