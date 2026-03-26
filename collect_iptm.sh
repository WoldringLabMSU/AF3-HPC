#!/usr/bin/env bash
# collect_iptm.sh
# Extract ipTM scores from AF3 PPI summary confidence files into a CSV.
#
# Usage:
#   bash collect_iptm.sh [output.csv]
#
# If no output file is given, writes to stdout.

BASE_DIR="/mnt/scratch/woldring/af3/outputs/ppi_ABY_run"
OUTPUT="${1:-/dev/stdout}"

echo "job_tag,iptm" > "${OUTPUT}"

for top_dir in "${BASE_DIR}"/binder*_target*/; do
    job_tag="$(basename "${top_dir}")"
    json="${top_dir}${job_tag}/${job_tag}_summary_confidences.json"

    if [[ ! -f "${json}" ]]; then
        continue
    fi

    iptm="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['iptm'])" "${json}")"
    echo "${job_tag},${iptm}" >> "${OUTPUT}"
done
