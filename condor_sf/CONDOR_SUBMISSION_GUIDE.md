# Condor Job Submission Guide

This guide covers how to submit trigger efficiency analysis jobs to HTCondor.

## Quick Start

### Recommended: Parallel Submission (Fastest)

Submit 4 jobs (one per period) that run in parallel:

```bash
cd /uscms/home/zzhao2/nobackup/hbb-run3/condor_sf
./submit_all_analyses.sh
```

This submits 4 condor jobs:
- Job 1: 2022 (2D SF + 1D SF)
- Job 2: 2022EE (2D SF + 1D SF)
- Job 3: 2023 (2D SF + 1D SF)
- Job 4: 2023BPix (2D SF + 1D SF)

### Alternative: Sequential Submission

Submit 1 job that processes all periods sequentially (slower):

```bash
./submit_all_analyses.sh --sequential
```

### Test Mode (Recommended First)

Always test with a small subset first (10 files per dataset):

```bash
./submit_all_analyses.sh --test
```

## Submission Options

### 1. Parallel Submission (Default, Recommended)

**Command:**
```bash
./submit_all_analyses.sh [--parallel] [--test]
```

**What it does:**
- Submits 4 separate condor jobs (one per data-taking period)
- Each job runs both 2D SF and 1D SF analyses
- Jobs run in parallel → **Much faster**
- Uses: `run_year-condor.jdl`

**Output locations:**
```
year_2022_<cluster_id>/
year_2022EE_<cluster_id>/
year_2023_<cluster_id>/
year_2023BPix_<cluster_id>/
```

**Use case:** Production runs, fastest turnaround

---

### 2. Sequential Submission

**Command:**
```bash
./submit_all_analyses.sh --sequential [--test]
```

**What it does:**
- Submits 1 condor job that processes all periods sequentially
- Runs 2D SF + 1D SF for each period, one after another
- Uses: `run_all-condor.jdl`

**Output locations:**
```
all_analyses_<cluster_id>/
```

**Use case:** When you want everything in one place, don't mind waiting

---

### 3. Single Year Submission

**Command:**
```bash
./submit_all_analyses.sh --year 2022 [--test]
```

**What it does:**
- Submits 1 condor job for a specific period only
- Runs both 2D SF and 1D SF for that period
- Uses: `run_year-condor.jdl` with YEAR parameter

**Output locations:**
```
year_2022_<cluster_id>/
```

**Use case:** Testing, rerunning a single period, debugging

---

## Direct Condor Submission (Advanced)

If you prefer to submit directly without the wrapper script:

### Parallel (4 jobs):
```bash
condor_submit run_year-condor.jdl
```

### Sequential (1 job):
```bash
condor_submit run_all-condor.jdl
```

### Single year:
```bash
condor_submit run_year-condor.jdl YEAR=2022
```

### Test mode:
```bash
condor_submit run_year-condor.jdl TEST_FLAG="--test"
condor_submit run_all-condor.jdl TEST_FLAG="--test"
```

## Monitoring Jobs

### Check job status:
```bash
condor_q                    # Summary view
condor_q -nobatch           # Detailed view
condor_q -analyze <job_id>  # Why job is idle/held
```

### Watch live output:
```bash
# Parallel jobs
tail -f year_2022_<cluster_id>/<cluster_id>_0.out
tail -f year_2022EE_<cluster_id>/<cluster_id>_0.out
tail -f year_2023_<cluster_id>/<cluster_id>_0.out
tail -f year_2023BPix_<cluster_id>/<cluster_id>_0.out

# Sequential job
tail -f all_analyses_<cluster_id>/<cluster_id>_0.out
```

### Check for errors:
```bash
# Look at error logs
cat year_2022_<cluster_id>/<cluster_id>_0.err

# Check condor log
cat year_2022_<cluster_id>/<cluster_id>_0.log
```

### Remove jobs:
```bash
condor_rm <cluster_id>      # Remove specific job
condor_rm <username>        # Remove all your jobs
```

## Output Structure

After jobs complete successfully, outputs will be transferred back:

### Parallel Jobs Output:
```
year_2022_<cluster_id>/
├── figures_sf_2d/          # 2D SF plots
│   └── 2022/
│       ├── VBF/
│       │   └── scale_factor_pt_vs_msd_2d.png
│       └── ggF/
│           └── scale_factor_pt_vs_msd_2d.png
├── output/                 # CSV tables and coffea files
│   └── 2022/
│       ├── VBF/
│       │   └── per_trigger_efficiency_2022_VBF.csv
│       └── ggF/
│           └── per_trigger_efficiency_2022_ggF.csv
├── sf_2022_OR/            # 1D SF plots
│   ├── scale_factor_pt_2022.png
│   ├── scale_factor_msd_2022.png
│   └── ...
├── sf_2d_2022/            # 2D SF logs
├── <cluster_id>_0.out     # stdout
├── <cluster_id>_0.err     # stderr
└── <cluster_id>_0.log     # condor log
```

Repeat for `year_2022EE_<cluster_id>/`, `year_2023_<cluster_id>/`, etc.

### Sequential Job Output:
```
all_analyses_<cluster_id>/
├── figures_sf_2d/          # 2D SF plots for ALL years
│   ├── 2022/
│   ├── 2022EE/
│   ├── 2023/
│   └── 2023BPix/
├── output/                 # CSV tables for ALL years
│   ├── 2022/
│   ├── 2022EE/
│   ├── 2023/
│   └── 2023BPix/
├── sf_2022_OR/            # 1D SF plots
├── sf_2022EE_OR/
├── sf_2023_OR/
├── sf_2023BPix_OR/
├── sf_2d_2022/            # Logs
├── sf_2d_2022EE/
├── sf_2d_2023/
├── sf_2d_2023BPix/
├── <cluster_id>_0.out
├── <cluster_id>_0.err
└── <cluster_id>_0.log
```

## Resource Requirements

Each job requests:
- **Memory**: 8 GB
- **CPUs**: 2
- **Time**: ~1-3 hours per year (full dataset)
- **Time**: ~5-10 minutes per year (test mode)

## Troubleshooting

### Job goes on hold:
```bash
condor_q -analyze <job_id>
cat <output_dir>/<cluster_id>_0.err
```

Common issues:
- **Memory exceeded**: Increase `request_memory` in JDL
- **Missing files**: Check `transfer_input_files` includes all needed scripts
- **Singularity error**: Check container name is correct

### Job completes but outputs missing:
```bash
# Check what was transferred back
cat <output_dir>/<cluster_id>_0.log | grep -i "transfer"

# Check for Python errors
cat <output_dir>/<cluster_id>_0.out | grep -i "error\|exception\|traceback"
```

Common issues:
- Output directories not created → Already fixed with `os.makedirs(..., exist_ok=True)`
- Empty datasets → Check JSON files in `infiles/`

### Job runs forever:
- Check if it's actually processing: `tail -f <output_dir>/<cluster_id>_0.out`
- If hung, remove and resubmit: `condor_rm <job_id>`
- Consider running in test mode first

## Recommended Workflow

### First Time:
1. **Test with small sample:**
   ```bash
   ./submit_all_analyses.sh --year 2022 --test
   ```
2. **Check outputs:**
   - Verify plots are generated
   - Check CSV tables look reasonable
   - Review logs for errors

3. **Run full production:**
   ```bash
   ./submit_all_analyses.sh
   ```

### Rerunning a Single Period:
```bash
./submit_all_analyses.sh --year 2023
```

### Quick Validation:
```bash
./submit_all_analyses.sh --test
```

## File Transfer Details

### What gets transferred TO the job:
- `shell` - Singularity wrapper
- `run_sf_2d.sh`, `run_sf.sh` - Execution scripts
- `scale_factor_2d.py`, `scale_factor.py` - Analysis code
- `infiles/` - JSON file lists (NOT the ROOT files)

### What gets transferred BACK:
- `figures_sf_2d/` - All PNG plots (2D)
- `output/` - CSV tables and coffea files
- `sf_<YEAR>_OR/` - All PNG plots (1D)
- Log directories

### NOT transferred (stays on worker node):
- ROOT files (read via XRootD)
- Temporary matplotlib cache
- Python bytecode (*.pyc)

## Comparison: Local vs Condor

### Running Locally:
```bash
./run_all_analyses.sh
```
- **Pros**: Immediate feedback, easier debugging
- **Cons**: Ties up your terminal, slower (sequential)
- **Time**: ~2-4 hours per year

### Running on Condor (Parallel):
```bash
./submit_all_analyses.sh
```
- **Pros**: Runs in background, parallel processing, faster
- **Cons**: Need to wait for scheduling, harder to debug
- **Time**: ~1-3 hours total (all 4 years in parallel)

### Running on Condor (Sequential):
```bash
./submit_all_analyses.sh --sequential
```
- **Pros**: All outputs in one place
- **Cons**: Slower than parallel, ties up worker node longer
- **Time**: ~4-8 hours total (all 4 years sequential)

## Best Practices

1. **Always test first** with `--test` flag
2. **Use parallel mode** for production runs
3. **Monitor jobs** - check logs periodically
4. **Clean up old outputs** before rerunning
5. **Save cluster IDs** for tracking jobs later
6. **Check resource limits** if jobs fail

## Summary

| Method | Command | Jobs | Time | Use Case |
|--------|---------|------|------|----------|
| Parallel | `./submit_all_analyses.sh` | 4 | ~1-3h | **Production (recommended)** |
| Sequential | `./submit_all_analyses.sh --sequential` | 1 | ~4-8h | All in one place |
| Single year | `./submit_all_analyses.sh --year 2022` | 1 | ~1-2h | Testing/rerun |
| Test mode | `./submit_all_analyses.sh --test` | 4 | ~10m | **Validation (recommended first)** |

For the Analysis Note, use **parallel mode** to get all results quickly!
