# AF3 Protein-Protein Co-Folding Workflow Guide

## Overview

This pipeline runs AlphaFold3 (AF3) protein-protein interaction (PPI) structure predictions on the MSU ICER HPC cluster. Given a CSV of binder and target protein sequences, it generates all pairwise co-folded structures (M binders × N targets) as a SLURM array job, then compiles the interface TM (ipTM) scores into a summary CSV.

---

## Prerequisites

The following resources must already exist before running the pipeline. Contact the lab if any are missing.

| Resource | Path |
|----------|------|
| AF3 Singularity container | `/mnt/research/woldring_lab/AlphaFold3/image/alphafold3.sif` |
| AF3 source code | `/mnt/research/woldring_lab/AlphaFold3/code/` |
| AF3 model weights | `/mnt/research/woldring_lab/AlphaFold3/weights/` |
| Sequence databases | `/mnt/research/common-data/alphafold/database_3` |
| This repo (cloned) | `/mnt/home/woldring/AF3-HPC/` |

---

## Workflow at a Glance

```
Step 1  Prepare your input CSV
           ↓
Step 2  python3 generate_ppi_jsons.py --csv <your.csv> --outdir <input_dir> --num-seeds 1
           ↓  (note the printed array range)
Step 3  Edit cofold_ppi_af3.sb  →  sbatch cofold_ppi_af3.sb
           ↓  (wait for jobs to finish)
Step 4  bash collect_iptm.sh iptm_scores.csv
```

---

## Step 1: Prepare Your Input CSV

Create a CSV file with two columns: `binders` and `targets`. Each row holds one protein sequence. The pipeline generates every binder–target combination (M × N total jobs).

**Example (`af3_inputs.csv`):**

```
binders,targets
AEAEYAKEPEYAVYEIDGLP...,GTDSGSEVLPDSFPSAPA...
MGHGGVSVYDQFDVFRGGL...,PKYVKQNTLKLAT...
ACDEFGHIKLMNPQRSTVWY...,
```

- Blank cells are skipped (columns do not need equal length).
- Duplicate sequences trigger a warning but are not removed.
- Save the file to `/mnt/home/woldring/AF3-HPC/` or any accessible path.

---

## Step 2: Generate JSON Input Files

Run `generate_ppi_jsons.py` to create one AF3 input JSON per binder-target pair.

```bash
cd ~/AF3-HPC

python3 generate_ppi_jsons.py \
    --csv /mnt/home/woldring/AF3-HPC/af3_inputs.csv \
    --outdir /mnt/scratch/woldring/af3/inputs/my_run_label \
    --num-seeds 1
```

**Arguments:**

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--csv` | Yes | — | Path to your input CSV |
| `--outdir` | Yes | — | Directory to write JSON files into. **Use a per-run subdirectory** (see [Running Multiple Batches](#running-multiple-batches-without-overwriting)) |
| `--num-seeds` | No | `1` | Seeds per job. Must match `AF3_NUM_SEEDS` in `cofold_ppi_af3.sb` |
| `--template` | No | `AF3_PPI.json` | Path to the JSON template (do not change unless you have a custom template) |
| `--skip-existing` | No | off | Skip JSON files that already exist in `--outdir` |

**Expected output:**
```
Generated 240 JSON files (80 binders x 3 targets).
Total combinations: 240
Set SLURM array to: 1-240%10
```

Note the array range printed (`1-240%10`) — you will need it in the next step.

---

## Step 3: Submit the SLURM Job

Open `cofold_ppi_af3.sb` and update the following variables before submitting.

### Variables to Edit

| Variable | Location in file | What to set |
|----------|-----------------|-------------|
| `RUN_LABEL` | line ~24 | A short unique name for this run, e.g. `"ppi_my_run"`. Controls both input and output subdirectories. |
| `CSV_PATH` | line ~33 | Full path to your input CSV |
| `#SBATCH --array` | line ~11 | The range printed by `generate_ppi_jsons.py`, e.g. `1-240%10` |

**Variables you should not need to change** (already configured for ICER):

| Variable | Value | Description |
|----------|-------|-------------|
| `AF3_PERSONAL_DIR` | `/mnt/home/woldring/AF3-HPC` | Repo directory |
| `AF3_RESOURCES_DIR` | `/mnt/research/woldring_lab/AlphaFold3` | Shared AF3 resources |
| `AF3_SCRATCH_DIR` | `/mnt/scratch/woldring/af3` | Root scratch directory |
| `AF3_NUM_SEEDS` | `1` | Seeds per job (match `--num-seeds` from Step 2) |
| `AF3_NUM_SAMPLES_PER_SEED` | `5` | Poses generated per seed |

> **Total structures per binder-target pair** = `AF3_NUM_SEEDS` × `AF3_NUM_SAMPLES_PER_SEED`
> (default: 1 × 5 = **5 structures per pair**)

### Submit

```bash
sbatch cofold_ppi_af3.sb
```

Monitor jobs:
```bash
squeue -u woldring
```

Cancel all jobs if needed:
```bash
scancel -u woldring
```

---

## Step 4: Collect ipTM Scores

Once jobs are complete, open `collect_iptm.sh` and confirm `BASE_DIR` matches your run's output directory:

```bash
BASE_DIR="/mnt/scratch/woldring/af3/outputs/my_run_label"
```

Then run:

```bash
bash collect_iptm.sh iptm_scores.csv
```

**Output CSV format:**

```
job_tag,iptm
binder1_target1,0.85
binder1_target2,0.31
binder2_target1,0.72
...
```

- Binder-target pairs where AF3 did not complete are silently skipped.
- `iptm` ranges from 0 to 1; higher values indicate better predicted binding.

---

## Directory Structure

| Directory | Purpose |
|-----------|---------|
| `/mnt/home/woldring/AF3-HPC/` | Repo root: scripts, template JSON, input CSVs |
| `/mnt/research/woldring_lab/AlphaFold3/image/` | `alphafold3.sif` Singularity container (~1–2 GB) |
| `/mnt/research/woldring_lab/AlphaFold3/code/` | AF3 Python source code |
| `/mnt/research/woldring_lab/AlphaFold3/weights/` | AF3 model parameter weights |
| `/mnt/research/common-data/alphafold/database_3` | Sequence databases (uniref90, uniprot, mgy, bfd, pdb_seqres) |
| `/mnt/scratch/woldring/af3/inputs/{RUN_LABEL}/` | Pre-generated JSON input files (one per binder-target pair) |
| `/mnt/scratch/woldring/af3/outputs/{RUN_LABEL}/` | AF3 structure outputs and confidence scores |
| `/mnt/scratch/woldring/af3/logs/` | SLURM `.out` / `.err` logs per array task |

> **Note:** `/mnt/scratch/` is not backed up and is periodically purged. Copy important results to `/mnt/research/` after runs complete.

---

## Key Variables Reference

A consolidated view of every variable a user needs to be aware of across all three files.

| Variable | File | Current Value | Description |
|----------|------|--------------|-------------|
| `RUN_LABEL` | `cofold_ppi_af3.sb` | `ppi_rfd3_run` | Unique run name — sets both input and output subdirs. **Change for every new run.** |
| `CSV_PATH` | `cofold_ppi_af3.sb` | `af3_inputs_from_rfd3.csv` | Full path to input CSV. **Change for every new run.** |
| `#SBATCH --array` | `cofold_ppi_af3.sb` | `1-240%10` | Array range from `generate_ppi_jsons.py` output. **Change for every new run.** |
| `AF3_NUM_SEEDS` | `cofold_ppi_af3.sb` | `1` | Seeds per job. Must match `--num-seeds` passed to `generate_ppi_jsons.py`. |
| `AF3_NUM_SAMPLES_PER_SEED` | `cofold_ppi_af3.sb` | `5` | Diffusion poses per seed. |
| `AF3_SCRATCH_DIR` | `cofold_ppi_af3.sb` | `/mnt/scratch/woldring/af3` | Root scratch directory. |
| `AF3_PERSONAL_DIR` | `cofold_ppi_af3.sb` | `/mnt/home/woldring/AF3-HPC` | Location of this repo. |
| `AF3_RESOURCES_DIR` | `cofold_ppi_af3.sb` | `/mnt/research/woldring_lab/AlphaFold3` | Shared AF3 resources directory. |
| `--csv` | `generate_ppi_jsons.py` | (required arg) | Input CSV path. |
| `--outdir` | `generate_ppi_jsons.py` | (required arg) | JSON output directory. Use `inputs/{RUN_LABEL}`. |
| `--num-seeds` | `generate_ppi_jsons.py` | `1` | Seeds embedded in each JSON. |
| `BASE_DIR` | `collect_iptm.sh` | `outputs/ppi_rfd3_run` | Must match `outputs/{RUN_LABEL}`. **Change for every new run.** |

---

## Running Multiple Batches Without Overwriting

Each run is isolated by `RUN_LABEL`. For every new batch, change **three things**:

| # | What to change | Where |
|---|---------------|-------|
| 1 | `RUN_LABEL` | `cofold_ppi_af3.sb` line ~24 |
| 2 | `CSV_PATH` | `cofold_ppi_af3.sb` line ~33 |
| 3 | `BASE_DIR` | `collect_iptm.sh` line 6 |

Also pass the matching `--outdir` and `--csv` to `generate_ppi_jsons.py`.

**Example for a third run:**

```bash
# generate_ppi_jsons.py
python3 generate_ppi_jsons.py \
    --csv /mnt/home/woldring/AF3-HPC/af3_inputs_batch3.csv \
    --outdir /mnt/scratch/woldring/af3/inputs/ppi_batch3 \
    --num-seeds 1

# cofold_ppi_af3.sb  (edit these lines)
RUN_LABEL="ppi_batch3"
CSV_PATH='/mnt/home/woldring/AF3-HPC/af3_inputs_batch3.csv'
#SBATCH --array=1-150%10   # from generate_ppi_jsons.py output

# collect_iptm.sh  (edit this line)
BASE_DIR="/mnt/scratch/woldring/af3/outputs/ppi_batch3"
```

---

## Expected Outputs

After a successful run, each completed binder-target pair produces the following under `outputs/{RUN_LABEL}/{job_tag}/`:

```
binder1_target1/
└── binder1_target1/
    ├── seed-0_sample-0/
    │   └── binder1_target1_model.cif       ← 3D structure (CIF format)
    ├── seed-0_sample-1/
    │   └── binder1_target1_model.cif
    ├── seed-0_sample-2/
    │   └── binder1_target1_model.cif
    ├── seed-0_sample-3/
    │   └── binder1_target1_model.cif
    ├── seed-0_sample-4/
    │   └── binder1_target1_model.cif
    ├── binder1_target1_summary_confidences.json
    ├── af3.stdout.log
    └── af3.stderr.log
```

**`_summary_confidences.json` fields:**

| Field | Description |
|-------|-------------|
| `iptm` | Interface TM-score — key metric for binding quality (0–1, higher is better) |
| `ptm` | Predicted TM-score for the full complex |
| `ranking_score` | AF3's internal ranking score |
| `has_clash` | Whether steric clashes were detected |
| `fraction_disordered` | Fraction of disordered residues |
| `chain_iptm` | Per-chain ipTM breakdown |
| `chain_ptm` | Per-chain pTM breakdown |

**SLURM logs** are written to `/mnt/scratch/woldring/af3/logs/`:
- `af3_ppi_{jobid}_{taskid}.out` — standard output
- `af3_ppi_{jobid}_{taskid}.err` — standard error

---

## Troubleshooting

### `ERROR: CSV not found`
The `CSV_PATH` in `cofold_ppi_af3.sb` is wrong or the file does not exist at that path. Verify the full path:
```bash
ls /mnt/home/woldring/AF3-HPC/your_file.csv
```

### `ERROR: JSON not found for binderX_targetY`
`generate_ppi_jsons.py` was not run before submitting, or `--outdir` does not match `AF3_INPUT_DIR` in the SLURM script. Re-run `generate_ppi_jsons.py` with `--outdir` set to `$AF3_SCRATCH_DIR/inputs/$RUN_LABEL`.

### `Unable to initialize backend 'rocm'` / `Unable to initialize backend 'tpu'`
These are harmless warnings. AF3 falls back to CUDA (NVIDIA GPU) automatically.

### `ValueError: GPU compute capability 7.x`
The job landed on a V100 GPU. AF3 requires compute capability 8.0+. The script requests `--gpus=a100:1` to avoid this. If you see this error, confirm the directive is present and resubmit.

### `KeyError: '...'` in `templates.py`
AF3 encounters a malformed template database entry. The script automatically detects this and retries with a patched `templates.py` that skips bad entries. If the retry also fails, check `af3.stderr.log` for a different error following the KeyError.

### Jobs finish immediately with exit code 0 but no output
The skip-if-complete check found existing CIF files and skipped the job. This is normal for resubmissions. To force a rerun, delete the output directory for that pair:
```bash
rm -rf /mnt/scratch/woldring/af3/outputs/{RUN_LABEL}/binderX_targetY/
```
