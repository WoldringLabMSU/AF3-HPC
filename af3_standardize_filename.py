import shutil
from pathlib import Path


root = Path("/mnt/scratch/patimara/AF3/outputs")
out_root = Path("/mnt/research/woldring_lab/Members/Patimar/standardized_names/AF3")

DRY_RUN = False

paths = root.glob("1b3_partially_overlapping_af3/*/*/seed-*")

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