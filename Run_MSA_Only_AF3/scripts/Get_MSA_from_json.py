import json 
from pathlib import Path
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--out-dir", type=Path, required=True, help="Path to output directory, e.g. ../outputs")
args = parser.parse_args()

input_dir = args.out_dir

if not input_dir.is_dir():
    parser.error(f"{input_dir} is not a directory")

for prot_dir in input_dir.iterdir():
    if not prot_dir.is_dir():
        continue
    
    json_files = list(prot_dir.rglob("*_data.json"))
    print(json_files)
    for json_file in json_files:
        with open (json_file) as f:
            data = json.load(f)

        try:
            protein = data["sequences"][0]["protein"]
            unpaired_msa = protein["unpairedMsa"]
            paired_msa = protein["pairedMsa"]
        except (KeyError, IndexError):
            print(f"Skipping {json_file}: not a valid AF3 JSON.")
            continue    

        unpaired_msa = data['sequences'][0]['protein']['unpairedMsa']
        paired_msa = data['sequences'][0]['protein']['pairedMsa']

        unpaired_msa = unpaired_msa.replace("\\n", "\n")
        paired_msa = paired_msa.replace("\\n", "\n")

        if not unpaired_msa.strip():
            print(f"WARNING: {json_file} has an empty unpaired MSA.")
            continue

        if not paired_msa.strip():
            print(f"WARNING: {json_file} has an empty paired MSA.")

        unpaired_MSA_file = json_file.with_name(f"{json_file.stem}_unpaired.a3m")
        paired_MSA_file = json_file.with_name(f"{json_file.stem}_paired.a3m")

        with open(unpaired_MSA_file, "w") as f:
            f.write(unpaired_msa)

        if paired_msa.strip():
            with open(paired_MSA_file, "w") as f:
                f.write(paired_msa)

        print(f"{json_file}: MSA files extracted and saved successfully.")
        print(f"Unpaired MSA saved to: {unpaired_MSA_file}")
        print(f"Paired MSA saved to: {paired_MSA_file}")
        print("--------------------------------------------------")
        