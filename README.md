# AlphaFold3 Batch Runner (SLURM)

This repository contains SLURM batch pipelines for running **AlphaFold3 structure predictions** on an HPC cluster — supporting both:

- **Protein-ligand** co-folding (one protein + one small molecule per job)
- **Protein-protein** co-folding (all-vs-all binder × target matrix)

---

# Protein–Ligand Co-folding

## Overview

Each SLURM array task:

1. Reads one row from the CSV
2. Generates a JSON input file
3. Runs AlphaFold3
4. Saves predictions and logs

Total structures per complex:
```
structures = AF3_NUM_SEEDS * AF3_NUM_SAMPLES_PER_SEED
```

Default:
```
20 seeds * 5 samples = 100 structures per complex
```

## Input CSV Format

The CSV **must contain these columns**:

```
pdb_id,smiles,sequence,ligand
```

## Running the Pipeline

Submit the job with:

```bash
sbatch cofold_af3_general.sb
```

The script uses a **SLURM job array**:

```
#SBATCH --array=41-134%10
```

Meaning:

- tasks **41–134** from the CSV
- **10 jobs run simultaneously**

Modify this range depending on your dataset.

## Output

Each complex produces a directory:

```
outputs/<RUN_LABEL>/<job_tag>/
```

---

# Protein–Protein Co-folding

## Overview

Runs AlphaFold3 on every combination of binders and targets from a CSV — an **M × N all-vs-all matrix**.

Two-step workflow:

1. **Pre-generate all JSON input files** with `generate_ppi_jsons.py`
2. **Submit the SLURM array job** with `cofold_ppi_af3.sb`

Total structures per pair:
```
structures = AF3_NUM_SEEDS * AF3_NUM_SAMPLES_PER_SEED
```

Default:
```
20 seeds * 5 samples = 100 structures per pair
```

## Input CSV Format

The CSV **must contain these columns**:

```
binders,targets
```

- Each non-blank entry in the `binders` column is one binder sequence.
- Each non-blank entry in the `targets` column is one target sequence.
- The columns do not need to be the same length.

Example:

```csv
binders,targets
MSEQBINDER1,MSEQTARGET1
MSEQBINDER2,MSEQTARGET2
MSEQBINDER3,
```

This produces a 3 binders × 2 targets = **6 total jobs**.

## Step 1 — Pre-generate JSON Files

Run this **once** before submitting the SLURM array:

```bash
python3 generate_ppi_jsons.py \
    --csv ppi.csv \
    --outdir /path/to/af3/inputs \
    [--template AF3_PPI.json] \
    [--num-seeds 20] \
    [--skip-existing]
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `--csv` | Yes | Input CSV with `binders` and `targets` columns |
| `--outdir` | Yes | Directory to write JSON files into |
| `--template` | No | Path to `AF3_PPI.json` (default: next to the script) |
| `--num-seeds` | No | Number of model seeds per job (default: 20) |
| `--skip-existing` | No | Skip writing a JSON if the file already exists |

The script prints the correct `--array` range to use in the next step:

```
Generated 6 JSON files (3 binders x 2 targets).
Total combinations: 6

Set SLURM array to: 1-6%10
```

Output files are named `binder{i}_target{j}.json` (1-based indices).

## Step 2 — Submit the SLURM Array

Edit the configurable variables at the top of `cofold_ppi_af3.sb`:

```bash
AF3_PERSONAL_DIR=/path/to/your/af3/dir   # contains AF3_PPI.json
AF3_RESOURCES_DIR=/path/to/AlphaFold3    # image, code, weights
AF3_SCRATCH_DIR=/path/to/scratch/AF3     # inputs, outputs, tmp

RUN_LABEL="my_ppi_run"                   # descriptive label for this run
CSV_PATH='/path/to/your/ppi.csv'
```

Update the array range to match the output of `generate_ppi_jsons.py`:

```bash
#SBATCH --array=1-6%10
```

Then submit:

```bash
sbatch cofold_ppi_af3.sb
```

## Output

Each binder-target pair produces a directory:

```
outputs/<RUN_LABEL>/binder{i}_target{j}/
```

Jobs that already have the expected number of output CIF files are automatically skipped (safe to re-submit after failures).

---

# Requirements

This pipeline assumes access to:
- **SLURM cluster**
- **GPU nodes**
- **AlphaFold3 container image**
- **AlphaFold3 weights**
- **AlphaFold3 databases**

Required software:

```
python3
singularity
bash
slurm
```

---

# Required AlphaFold3 Resources

You must provide paths to:

```
AF3_IMAGE            AlphaFold3 container image
AF3_CODE_DIR         AlphaFold3 code directory
AF3_MODEL_PARAMETERS model weights
AF3_DATABASES_DIR    AF3 sequence + structure databases
```

Example from the scripts:

```
AF3_RESOURCES_DIR=/mnt/research/woldring_lab/AlphaFold3
```

---

# Configurable Variables

### Input CSV

```
CSV_PATH=/path/to/input.csv
```

### Output label

```
RUN_LABEL="experiment_name"
```

### Resource directories

```
AF3_PERSONAL_DIR
AF3_RESOURCES_DIR
AF3_SCRATCH_DIR
```

### Number of samples

```
AF3_NUM_SEEDS=20
AF3_NUM_SAMPLES_PER_SEED=5
```

---

Good luck!
