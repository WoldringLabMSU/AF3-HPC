"""
Date: March 2026
Description: Organization of outputs from an AF3 SLURM job run. Extracts the structural files and json confidence scores.
Then stores them in a fresh directory that includes the name of the ligand and a number. Then stores them in another directory
with the ligand name. Repeats until all outputs have been processed.
"""
import shutil
from pathlib import Path


# CHANGE ME!!!
root = Path(config["data_path"]) # the path to your desired directory where outputs are located 
out_root = Path(config["out_path"]) # the path you want your NEW outputs to go in

# TRUE: if you want to see the EXPECTED output printed in the terminal but not officially run the code
# FALSE: if you want to run the script fully
DRY_RUN = False

# CHANGE ME!!! Only change the first line in the path to the directory that comes right after your current path in 'root'
# e.g. 1b3_completely_overlapping
paths = root.glob("<CHANGE_ME>/*/*/seed-*")

# to count number of poses: 20 seeds x 5 runs for each = 100 poses for each ligand
global_num = 1

for i, p in enumerate(sorted(paths)):
    if DRY_RUN and i > 2:
        break
    if not p.is_dir():
        continue
    if DRY_RUN and global_num > 10:
        break
    # extract necessary names
    ligand_folder = p.parent.parent.name
    lig = ligand_folder.split('_', 1)[0]
    lig = lig.lower()

    oatp_raw = p.parent.name
    oatp = oatp_raw.upper()

    oatp_dir = out_root / oatp
    oatp_dir.mkdir(parents=True, exist_ok=True)

    lig_dir = oatp_dir / lig
    lig_dir.mkdir(parents=True, exist_ok=True)

    new_base = f"{oatp}_{lig}_AF3_{global_num}"
    out_dir = lig_dir / new_base

    if DRY_RUN: 
        print(f"[MKDIR] {out_dir.name}, parent: {out_dir.parent.name}")
    else:
        out_dir.mkdir(exist_ok=True)
    
    structure = p / "model.cif"
    confidence_scores = p / "summary_confidences.json"

    new_struct_path = out_dir / f"{new_base}.cif"
    new_confidence_path = out_dir / f"{new_base}.json"

    if not structure.exists() or not confidence_scores.exists():
        print(f"[SKIP] Missing files in {p}")
        continue

    if DRY_RUN:
        print(f"[COPY] {structure.name} -> {new_struct_path.name}")
        print(f"[COPY] {confidence_scores.name} -> {new_confidence_path.name}")
    else:
        shutil.copy(structure, new_struct_path)
        shutil.copy(confidence_scores, new_confidence_path)

    print(f"[INDEX {global_num}] {lig}")
    global_num += 1
    if global_num > 100:
        global_num = 1