#!/bin/bash --login
#SBATCH --job-name=MSAonly_AF3
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --output=/mnt/scratch/needsang/AF3-HPC/Run_MSA_Only_AF3/logs/%x_%A_%a.out             # CHANGE ME
#SBATCH --error=/mnt/scratch/needsang/AF3-HPC/Run_MSA_Only_AF3/logs/%x_%A_%a.err              # CHANGE ME
#SBATCH --array=1-1                                                                           # CHANGE ME

# ======================================================================
# CHANGE THE ARRAY, error log output, and other directories 
# (in particualr the input dir)
# ======================================================================

set -euo pipefail

module purge

export AF3_RESOURCES_DIR=/mnt/research/woldring_lab/AlphaFold3        
export AF3_SCRATCH_DIR=/mnt/scratch/needsang/AF3-HPC/Run_MSA_Only_AF3        # CHANGE ME

export AF3_IMAGE="${AF3_RESOURCES_DIR}/image/alphafold3.sif"                
export AF3_CODE_DIR="${AF3_RESOURCES_DIR}/code"
export AF3_INPUT_DIR="${AF3_SCRATCH_DIR}/json_dir"
export AF3_OUTPUT_DIR="${AF3_SCRATCH_DIR}/outputs" 
export AF3_MODEL_PARAMETERS="${AF3_RESOURCES_DIR}/weights"
export AF3_DATABASES_DIR=/mnt/research/common-data/alphafold/database_3


# 0 = CPU-only singularity exec, recommended for data-pipeline-only mode
# 1 = inlcude --nv, useful if container requires it for some reason
AF3_USE_NV="${AF3_USE_NV:-0}"   

# Extra run_alphafold.py arguments to supply at submission:
#   e.g., sbatch --export=ALL,AF3_EXTRA_ARGS='--jackhammer_n_cpu=8' ...
AF3_EXTRA_ARGS="${AF3_EXTRA_ARGS:-}"

AF3_TMP_HOST="${AF3_SCRATCH_DIR}/tmp"


mkdir -p "${AF3_OUTPUT_DIR}" "${AF3_SCRATCH_DIR}/logs" "${AF3_TMP_HOST}"

[[ -d "${AF3_INPUT_DIR}" ]] || { echo "ERROR: AF3_INPUT_DIR not found: ${AF3_INPUT_DIR}" >&2; exit 1; }
[[ -f "${AF3_IMAGE}" ]] || { echo "ERROR: AF3 image not found: ${AF3_IMAGE}" >&2; exit 1; }
[[ -d "${AF3_CODE_DIR}" ]] || { echo "ERROR: AF3 code dir not found: ${AF3_CODE_DIR}" >&2; exit 1; }
[[ -d "${AF3_MODEL_PARAMETERS}" ]] || { echo "ERROR: AF3 model parameters dir not found: ${AF3_MODEL_PARAMETERS}" >&2; exit 1; }
[[ -d "${AF3_DATABASES_DIR}" ]] || { echo "ERROR: AF3 databases dir not found: ${AF3_DATABASES_DIR}" >&2; exit 1; }

mapfile -t JSONS < <(find "${AF3_INPUT_DIR}" -maxdepth 1 -type f -name '*.json' | sort)
N_JSONS="${#JSONS[@]}"

if [[ "${N_JSONS}" -lt 1 ]]; then
    echo "ERROR: no JSON files found in ${AF3_INPUT_DIR}" >&2
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
if (( TASK_ID < 1 || TASK_ID > N_JSONS )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${TASK_ID} outside valid range 1-${N_JSONS}" >&2
    exit 1
fi

JSON_HOST="${JSONS[$((TASK_ID - 1))]}"
job_tag="$(basename "${JSON_HOST}" .json)"

AF3_STDOUT="${AF3_SCRATCH_DIR}/logs/${job_tag}.msa_only.stdout.log"
AF3_STDERR="${AF3_SCRATCH_DIR}/logs/${job_tag}.msa_only.stderr.log"

echo "AF3 MSA-only/data-pipeline-only task"
echo "Task ID: ${TASK_ID}/${N_JSONS}"
echo "Job tag: ${job_tag}"
echo "Input JSON: ${JSON_HOST}"
echo "Output root: ${AF3_OUTPUT_DIR}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-unset}"
echo "TMPDIR: ${AF3_TMP_HOST}"

SINGULARITY_ARGS=()
if [[ "${AF3_USE_NV}" == "1" ]]; then
    SINGULARITY_ARGS+=(--nv)
fi 


set +e
singularity exec "${SINGULARITY_ARGS[@]}" \
    --bind "${AF3_INPUT_DIR}:/root/af_input" \
    --bind "${AF3_OUTPUT_DIR}:/root/af_output" \
    --bind "${AF3_MODEL_PARAMETERS}:/root/models" \
    --bind "${AF3_DATABASES_DIR}:/root/public_databases" \
    --bind "${AF3_CODE_DIR}:/root/af_code" \
    --bind "${AF3_TMP_HOST}:/tmp/local" \
    --env TMPDIR=/tmp/local \
    "${AF3_IMAGE}" \
    python3 /root/af_code/alphafold3/run_alphafold.py \
        --json_path "/root/af_input/${job_tag}.json" \
        --model_dir /root/models \
        --db_dir /root/public_databases \
        --output_dir "/root/af_output/${job_tag}" \
        --run_data_pipeline=true \
        --run_inference=false \
        ${AF3_EXTRA_ARGS} \
        > >(tee -a "${AF3_STDOUT}") \
        2> >(tee -a "${AF3_STDERR}" >&2)
status=$?
set -e

if [[ "${status}" -ne 0 ]]; then
    echo "ERROR: AF3 MSA-only run failed for ${job_tag} with status ${status}" >&2
    exit "${status}"
fi

DATA_JSON_COUNT=$(find "${AF3_OUTPUT_DIR}/${job_tag}" -type f -name '*_data.json' 2>/dev/null | wc -l | tr -d ' ')
if [[ "${DATA_JSON_COUNT}" -lt 1 ]]; then
    echo "ERROR: AF3 completed but no *_data.json was found under ${AF3_OUTPUT_DIR}/${job_tag}" >&2
    exit 2
fi

echo "Completed MSA-only AF3 data pipeline for ${job_tag}."
find "${AF3_OUTPUT_DIR}/${job_tag}" -type f -name '*_data.json' | sort
