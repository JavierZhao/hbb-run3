# Trigger Scale Factor Methodology

## Key Principle: Triggers are Production-Mode Agnostic

**Trigger scale factors are NOT production-mode specific!**

Triggers like:
- `HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35`
- `HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65`

These fire based on:
- ✅ Jet kinematics (pT, mSD, η)
- ✅ ParticleNet discriminator scores
- ✅ Number and configuration of jets

They do **NOT** fire based on:
- ❌ Production mechanism (VBF vs ggF)
- ❌ Event topology beyond trigger objects
- ❌ Generator-level information

## Correct Methodology (Implemented)

### 1. Trigger SF Measurement (Tag-and-Probe)

**Sample:** TTtoLNu2Q (semileptonic ttbar)
**Selection:**
- Muon reference triggers (HLT_Mu50, HLT_IsoMu24, etc.)
- ≥1 loose muon (pT > 25 GeV, |η| < 2.4, isolated, ΔR > 0.8 from FatJet)
- Leading FatJet (pT > 250 GeV, |η| < 2.5, ParticleNet_XbbVsQCD > 0.4)
- **NO VBF/ggF separation** - measured inclusively!

**Binning:** 2D in (jet pT, jet mSD)
- pT: 300-1000 GeV (plotting), ≥450 GeV (efficiency calculation)
- mSD: 0-300 GeV

**Output:** One SF map as a function of (pT, mSD)

```
SF(pT, mSD) = ε_data(pT, mSD) / ε_MC(pT, mSD)
```

### 2. Application in Signal Analysis

In the H→bb signal analysis:

**Step 1:** Categorize events into VBF and ggF based on topology
- VBF: Δη(j1,j2) > 3.5, mjj > 1000 GeV
- ggF: All other events

**Step 2:** Apply the **same SF map** to both categories
- For a VBF event with jet (pT=600 GeV, mSD=120 GeV): look up SF(600, 120)
- For a ggF event with jet (pT=600 GeV, mSD=120 GeV): look up SF(600, 120)
- **Same SF value** for both!

## Why VBF Selection Was Wrong

### Previous (Incorrect) Implementation:
```python
prod_modes = ['VBF', 'ggF']  # ❌ Created separate SFs
# Applied VBF topology cuts (Δη > 3.5, mjj > 1000) during SF measurement
```

**Problems:**
1. ❌ VBF cuts **suppress ttbar** (by design!) - only 1,994 data events survived
2. ❌ Low statistics → unreliable SF (~18% data eff vs ~93% MC eff)
3. ❌ Triggers don't care about production mode
4. ❌ Deviates from published CMS methodology

### Current (Correct) Implementation:
```python
prod_modes = ['Inclusive']  # ✅ One SF map for all events
# NO VBF topology cuts during SF measurement
```

**Benefits:**
1. ✅ High statistics (all ttbar events with muon tag)
2. ✅ Reliable efficiency measurement
3. ✅ Follows published CMS methodology
4. ✅ Same SF applies to VBF and ggF based on jet kinematics

## Supporting Evidence

### From Published CMS H→bb Analysis (arXiv:2407.08012)

> "The trigger selection is about 90% efficient with respect to the offline selection for events containing large-radius jets with pt from 450–500 GeV and |η|<2.5, and is fully efficient for pt>500 GeV and |η|<2.5."

**Note:** No distinction between VBF and ggF for trigger efficiency!

### From VBF Topology Studies (arXiv:2206.04965)

> "These selections **suppress the dominant background coming from top-quark production** [ttbar], which is also reduced by requiring a b-veto on the tag jets..."

**This proves:** VBF cuts are designed to **reject** ttbar, making VBF+ttbar incompatible!

### From CMS Tag-and-Probe Documentation

Tag-and-probe is designed for:
- Lepton efficiencies (using Z→ll, J/ψ→ll resonances)
- Trigger efficiencies that depend only on object kinematics
- **NOT** for topology-dependent categorization

## Output Structure

### Files Generated

**Plots:**
```
figures_sf_2d/{YEAR}/Inclusive/scale_factor_pt_vs_msd_2d.png
```

**CSV Tables:**
```
output/{YEAR}/Inclusive/per_trigger_efficiency_{YEAR}_Inclusive.csv
```

**Coffea Files:**
```
output/sf_2d_{YEAR}_Inclusive_data.coffea
output/sf_2d_{YEAR}_Inclusive_mc.coffea
```

### CSV Table Format

```csv
Trigger,Data_Efficiency,Data_Pass,Data_Total,MC_Efficiency,MC_Pass,MC_Total,Scale_Factor
AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35,0.8542,12345,14456,0.8923,9876,11067,0.9573
QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65,0.7234,10456,14456,0.7689,8512,11067,0.9408
Combined (OR),0.9123,13189,14456,0.9456,10465,11067,0.9648
```

Efficiencies calculated in kinematic region: **450 ≤ pT < 1000 GeV, 0 ≤ mSD < 300 GeV**

## How to Use These SFs

### In Your H→bb Analysis:

```python
# 1. Load the SF map (2D histogram or lookup table)
sf_map = load_sf_map("sf_2d_2022_Inclusive.root")

# 2. For each event (VBF or ggF), look up SF based on jet kinematics
for event in events:
    jet_pt = event.FatJet.pt[0]
    jet_msd = event.FatJet.msoftdrop[0]

    # Same SF lookup for both VBF and ggF!
    sf = sf_map.GetBinContent(
        sf_map.FindBin(jet_msd, jet_pt)
    )

    # Apply SF to event weight
    event_weight *= sf
```

### Expected SF Values

From the 2022EE VBF measurement (which was incorrect):
- Data efficiency: 22%
- MC efficiency: 87%
- SF: 0.26 ← **Too low due to poor statistics!**

Expected from inclusive measurement:
- Data efficiency: ~80-90%
- MC efficiency: ~90-95%
- SF: ~0.95-1.00 ← More reasonable!

## Changes Made

### Code Changes:

1. **scale_factor_2d.py** (condor_sf and trig_eff versions)
   - Line ~616: `prod_modes = ['Inclusive']` (was `['VBF', 'ggF']`)
   - Added comments explaining methodology
   - Output directories now use 'Inclusive' naming

2. **Baseline selection function**
   - Added comment: "VBF cuts should NOT be used for trigger SF"
   - VBF cuts only applied when explicitly requested (not for 'Inclusive')

### Documentation Updates:

1. **README_ANALYSIS_OUTPUTS.md** - Updated to reflect inclusive measurement
2. **TRIGGER_SF_METHODOLOGY.md** (this file) - Complete methodology documentation

## Running the Updated Analysis

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf

# Test run (recommended first)
./submit_all_analyses.sh --year 2022 --test

# Full production
./submit_all_analyses.sh
```

**Expected improvements:**
- ✅ Much higher statistics (10-20x more events)
- ✅ More reliable efficiencies (both data and MC should be ~90%)
- ✅ Scale factors closer to 1.0
- ✅ Smaller uncertainties
- ✅ Methodology matches published analyses

## References

1. [CMS boosted H→bb measurement (arXiv:2407.08012)](https://arxiv.org/html/2407.08012)
2. [CMS PAS HIG-23-012](https://cds.cern.ch/record/2904879/files/HIG-23-012-pas.pdf)
3. [VBF topology and ttbar suppression (arXiv:2206.04965)](https://arxiv.org/abs/2206.04965)
4. [CMS Tag-and-Probe Documentation](https://cms-opendata-workshop.github.io/workshop-lesson-tagandprobe/aio/index.html)
5. Run 2 Trigger Efficiency AN (TRG_AN_RUN_2.pdf) - methodology slide

## Summary

**Key Takeaway:** Trigger SFs are measured **inclusively** in bins of (pT, mSD), with no VBF/ggF separation. The same SF map applies to all H→bb events regardless of production mechanism. VBF topology cuts should never be used during trigger SF measurement with muon tag-and-probe.
