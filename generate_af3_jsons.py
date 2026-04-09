#!/usr/bin/env python3
"""
AF3 JSON generator for the CPU/GPU split pipeline.

Run this script TWICE:

  Stage 1 -- before the CPU data-pipeline jobs:
    python generate_af3_jsons.py --stage protein --csv_paths 1b1/*.csv 1b3/*.csv 2b1/*.csv

    Writes one protein-only JSON per unique OATP and oatp_json_list.txt.
    Then submit af3_data_pipeline.sb.

  Stage 2 -- after the CPU data-pipeline jobs finish:
    python generate_af3_jsons.py --stage pairs --csv_paths 1b1/*.csv 1b3/*.csv 2b1/*.csv

    Reads the precomputed JSONs (MSA + templates filled in), writes one
    protein+ligand JSON per OATP-ligand pair, and pair_json_list.txt.
    Then submit af3_inference.sb.

Default paths assume the layout:
  /mnt/research/woldring_lab/af3/   -- JSON staging area and list files
  /mnt/scratch/woldring/af3/        -- AF3 outputs (data pipeline + inference)
"""

import argparse
import csv
import json
import sys
import zlib
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_seeds(name: str, n_seeds: int) -> list:
    """Deterministic CRC32-based seeds -- matches cofold_af3_general.sb."""
    max_seed = 2147483647
    base = (zlib.crc32(name.encode("utf-8")) % (max_seed - 1)) + 1
    seeds = []
    for i in range(n_seeds):
        s = (base + i) % max_seed
        if s == 0:
            s = 1
        seeds.append(int(s))
    return seeds


def sanitize(name: str) -> str:
    """Strip characters that would break filesystem paths."""
    return "".join(c if (c.isalnum() or c in "-_.") else "_" for c in name)


def read_csvs(csv_paths: list) -> dict:
    """
    Parse one or more input CSVs (columns: pdb_id, smiles, sequence, ligand).

    Returns:
        {pdb_id: {"sequence": str, "ligands": [{"smiles": str, "name": str}]}}

    The protein sequence is taken from the first row seen for each pdb_id;
    all ligand rows for that pdb_id are collected from every CSV.
    """
    proteins: dict = {}
    needed = {"pdb_id", "smiles", "sequence", "ligand"}

    for csv_path in csv_paths:
        csv_path = Path(csv_path)
        if not csv_path.exists():
            print(f"WARNING: CSV not found, skipping: {csv_path}", file=sys.stderr)
            continue

        with open(csv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            missing = needed - set(reader.fieldnames or [])
            if missing:
                print(
                    f"WARNING: {csv_path} is missing columns {sorted(missing)}, skipping.",
                    file=sys.stderr,
                )
                continue

            for row in reader:
                pdb_id   = row["pdb_id"].strip()
                sequence = "".join(row["sequence"].split())   # collapse whitespace
                smiles   = row["smiles"].strip()
                ligand   = row["ligand"].strip()

                if not pdb_id or not sequence or not smiles:
                    continue

                if pdb_id not in proteins:
                    proteins[pdb_id] = {"sequence": sequence, "ligands": []}

                proteins[pdb_id]["ligands"].append({"smiles": smiles, "name": ligand})

    return proteins


def find_precomputed_json(precomputed_dir: Path, pdb_id: str) -> Path:
    """
    Locate the augmented JSON written by AF3 after --run_inference=false.

    AF3 writes it to {output_dir}/{name}/{name}.json.  We also try a glob
    fallback in case the filename differs slightly.
    """
    job_dir = precomputed_dir / pdb_id
    if not job_dir.is_dir():
        return None

    # Primary candidate
    candidate = job_dir / f"{pdb_id}.json"
    if candidate.exists():
        return candidate

    # Fallback: any JSON in the directory that isn't a confidence/summary file
    for p in sorted(job_dir.glob("*.json")):
        if "confidence" not in p.stem and "summary" not in p.stem:
            return p

    return None


# ---------------------------------------------------------------------------
# Stage 1: protein-only JSONs
# ---------------------------------------------------------------------------

def write_protein_json(pdb_id: str, sequence: str, out_path: Path, n_seeds: int) -> None:
    job = {
        "name": pdb_id,
        "modelSeeds": make_seeds(pdb_id, n_seeds),
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence,
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(job, fh, indent=2)


def stage_protein(args) -> None:
    csv_paths    = args.csv_paths
    out_dir      = Path(args.protein_json_dir)
    list_file    = Path(args.protein_list)
    n_seeds      = args.n_seeds

    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_csvs(csv_paths)
    if not proteins:
        print("ERROR: no proteins found in any CSV.", file=sys.stderr)
        sys.exit(1)

    json_paths = []
    for pdb_id in sorted(proteins):
        sequence = proteins[pdb_id]["sequence"]
        out_path = out_dir / f"{pdb_id}.json"
        write_protein_json(pdb_id, sequence, out_path, n_seeds)
        json_paths.append(str(out_path.resolve()))
        print(f"  wrote {out_path}  ({len(proteins[pdb_id]['ligands'])} ligands will be paired later)")

    list_file.parent.mkdir(parents=True, exist_ok=True)
    list_file.write_text("\n".join(json_paths) + "\n")

    n = len(json_paths)
    print(f"\nWrote {n} protein JSON(s) to {out_dir}")
    print(f"Job list: {list_file}")
    print(f"\nNext step:")
    print(f"  Edit af3_data_pipeline.sb  -->  set  #SBATCH --array=1-{n}")
    print(f"  Then:  sbatch af3_data_pipeline.sb")


# ---------------------------------------------------------------------------
# Stage 2: protein+ligand pair JSONs
# ---------------------------------------------------------------------------

def write_pair_json(
    protein_entry: dict,
    pair_name: str,
    smiles: str,
    out_path: Path,
    n_seeds: int,
) -> None:
    job = {
        "name": pair_name,
        "modelSeeds": make_seeds(pair_name, n_seeds),
        "sequences": [
            protein_entry,                        # protein with MSA/templates from data pipeline
            {
                "ligand": {
                    "id": "B",
                    "smiles": smiles,
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(job, fh, indent=2)


def stage_pairs(args) -> None:
    csv_paths       = args.csv_paths
    precomputed_dir = Path(args.precomputed_dir)
    out_dir         = Path(args.pair_json_dir)
    list_file       = Path(args.pair_list)
    n_seeds         = args.n_seeds

    out_dir.mkdir(parents=True, exist_ok=True)

    proteins = read_csvs(csv_paths)
    if not proteins:
        print("ERROR: no proteins found in any CSV.", file=sys.stderr)
        sys.exit(1)

    json_paths: list = []
    skipped_proteins: list = []

    for pdb_id in sorted(proteins):
        precomputed_json = find_precomputed_json(precomputed_dir, pdb_id)
        if precomputed_json is None:
            print(
                f"WARNING: no precomputed data-pipeline JSON for {pdb_id} "
                f"in {precomputed_dir}. "
                f"Has af3_data_pipeline.sb finished for this protein?",
                file=sys.stderr,
            )
            skipped_proteins.append(pdb_id)
            continue

        with open(precomputed_json) as fh:
            precomputed = json.load(fh)

        # Pull out the protein entry -- it now contains unpairedMsa, pairedMsa, templates
        protein_entry = None
        for entry in precomputed.get("sequences", []):
            if "protein" in entry:
                protein_entry = entry
                break

        if protein_entry is None:
            print(
                f"WARNING: precomputed JSON for {pdb_id} has no protein entry, skipping.",
                file=sys.stderr,
            )
            skipped_proteins.append(pdb_id)
            continue

        n_ligands = len(proteins[pdb_id]["ligands"])
        print(f"  {pdb_id}: generating {n_ligands} pair JSON(s)  (source: {precomputed_json.name})")

        for lig in proteins[pdb_id]["ligands"]:
            ligand_name = sanitize(lig["name"])
            smiles      = lig["smiles"]
            # Match naming convention from cofold_af3_general.sb: ${ligand}_OATP${pdb_id}_af3
            pair_name   = f"{ligand_name}_OATP{pdb_id}_af3"
            out_path    = out_dir / f"{pair_name}.json"
            write_pair_json(protein_entry, pair_name, smiles, out_path, n_seeds)
            json_paths.append(str(out_path.resolve()))

    if skipped_proteins:
        print(
            f"\nWARNING: {len(skipped_proteins)} protein(s) had no precomputed data "
            f"and were skipped: {skipped_proteins}",
            file=sys.stderr,
        )

    list_file.parent.mkdir(parents=True, exist_ok=True)
    list_file.write_text("\n".join(json_paths) + "\n")

    n = len(json_paths)
    print(f"\nWrote {n} pair JSON(s) to {out_dir}")
    print(f"Job list: {list_file}")
    print(f"\nNext step:")
    print(f"  Edit af3_inference.sb  -->  set  #SBATCH --array=1-{n}%25")
    print(f"  Then:  sbatch af3_inference.sb")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--stage",
        required=True,
        choices=["protein", "pairs"],
        help=(
            "protein: write one protein-only JSON per unique OATP (run before CPU jobs). "
            "pairs: write one protein+ligand JSON per OATP-ligand pair (run after CPU jobs)."
        ),
    )
    parser.add_argument(
        "--csv_paths",
        nargs="+",
        required=True,
        metavar="CSV",
        help="One or more CSV files with columns: pdb_id, smiles, sequence, ligand.",
    )
    parser.add_argument(
        "--n_seeds",
        type=int,
        default=20,
        help="Model seeds per complex (default: 20, matching cofold_af3_general.sb).",
    )

    # Stage 1 outputs
    parser.add_argument(
        "--protein_json_dir",
        default="/mnt/research/woldring_lab/af3/oatp_protein_jsons",
        help="Output dir for protein-only JSONs [stage=protein].",
    )
    parser.add_argument(
        "--protein_list",
        default="/mnt/research/woldring_lab/af3/oatp_json_list.txt",
        help="Output path for the protein JSON list file [stage=protein].",
    )

    # Stage 2 inputs/outputs
    parser.add_argument(
        "--precomputed_dir",
        default="/mnt/scratch/woldring/af3/precomputed_oatps",
        help=(
            "Directory where af3_data_pipeline.sb wrote its output JSONs "
            "(each protein gets a subdirectory here) [stage=pairs]."
        ),
    )
    parser.add_argument(
        "--pair_json_dir",
        default="/mnt/research/woldring_lab/af3/pair_jsons",
        help="Output dir for protein+ligand pair JSONs [stage=pairs].",
    )
    parser.add_argument(
        "--pair_list",
        default="/mnt/research/woldring_lab/af3/pair_json_list.txt",
        help="Output path for the pair JSON list file [stage=pairs].",
    )

    args = parser.parse_args()

    if args.stage == "protein":
        stage_protein(args)
    else:
        stage_pairs(args)


if __name__ == "__main__":
    main()
