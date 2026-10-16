# Run AlphaFold 3 MSA-Only Pipeline

This subdirectory runs AlphaFold 3 (AF3) to generate multiple sequence alignments (MSAs) only. No structure inference is run. Intermediate files and output directories are created and saved within this subdirectory to keep the main repository tidy.

## Required input

A CSV file with one row per protein and two columns:

| Column       | Description                |
|--------------|----------------------------|
| `protein_id` | Unique ID for each protein |
| `sequence`   | Amino acid sequence        |

Exact column names can vary slightly. If yours differ, pass them as arguments to `make_af3_msa_precompute.py`.

## Workflow

```bash
cd AF3-HPC/Run_MSA_Only_AF3/scripts/

# 1. Generate one AF3 input JSON per protein sequence
python make_af3_msa_precompute.py --seq-list-csv {input_csv_file}

# 2. Run the MSA-only AF3 job on the cluster
sbatch af3_apo_MSAonly.sh

# 3. Convert the MSA JSONs to .a3m files
python Get_MSA_from_json.py --out-dir ../outputs
```

> **Before running step 2:** edit `af3_apo_MSAonly.sh` to
> - change all 3 instances of the `SCRATCH_DIR` definition (replace with your username), and
> - set the Slurm array size to match the number of proteins in your CSV.

## Outputs

**1. `make_af3_msa_precompute.py`** creates an intermediate directory containing:
- one JSON file per protein sequence
- a manifest text file
- a text file listing all intermediate JSON files

**2. `af3_apo_MSAonly.sh`** creates an output directory (`outputs/`) with:
- one subdirectory per protein sequence
- within each, another subdirectory containing the `.json` file with the MSA data

**3. `Get_MSA_from_json.py`** writes paired and unpaired MSAs in `.a3m` format to the same directory as each protein's JSON file.

## Resulting directory structure

```
outputs/
├── protein_A/
│   └── protein_a/
│       ├── <MSA data>.json
│       ├── <paired MSA>.a3m
│       └── <unpaired MSA>.a3m
└── protein_B/
    └── ...
```