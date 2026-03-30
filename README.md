# AlphaFold3 Batch Runner (SLURM)

This repository contains a SLURM batch pipeline for running **AlphaFold3 protein-ligand structure predictions** across many complexes using a CSV input file.

The script automatically:

- Reads **protein sequences and ligand SMILES** from a CSV file.
- Generates **AlphaFold3 JSON input files**
- Runs AlphaFold3 inside a **singularity container**
- Launches **parallel SLURM array jobs**
- Produces multiple **diffusion samples per seed**
- Automatically skips completed jobs

---

# Overview

Each SLURM array task:

1. Reads one row from the CSV (pdb_id,smiles,sequence,ligand)
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
You can change this number based on your preference.

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

---

# Input CSV Format

The CSV **must contain these columns**:

```
pdb_id,smiles,sequence,ligand
```

---

# Configurable Variables

You will likely need to modify the following in the script.

### Input CSV

```
CSV_PATH=/path/to/input.csv
```

---

### Output label

```
RUN_LABEL="experiment_name"
```

---

### Resource directories

```
AF3_PERSONAL_DIR
AF3_RESOURCES_DIR
AF3_SCRATCH_DIR
```

---

### Number of samples

```
AF3_NUM_SEEDS=20
AF3_NUM_SAMPLES_PER_SEED=5
```

---

# Running the Pipeline

Submit the job with:

```
sbatch cofold_af3_general.sh
```

---

The script uses a **SLURM job array**:

```
#SBATCH --array=1-100%10
```

Meaning:

- tasks **1-100** from the CSV
- **10 jobs run simultaneously**

Modify this range depending on your dataset.

---

# Output

Each complex produces a directory:

```
outputs/<RUN_LABEL>/<job_tag>/
```

---

Good luck!
