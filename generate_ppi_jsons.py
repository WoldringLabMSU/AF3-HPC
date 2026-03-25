#!/usr/bin/env python3
"""
generate_ppi_jsons.py

Pre-generate all M x N AlphaFold3 JSON input files for protein-protein
co-folding from a CSV with 'binders' and 'targets' columns.

Run this script before submitting the SLURM array job.

Usage:
    python3 generate_ppi_jsons.py \\
        --csv ppi.csv \\
        --outdir /path/to/af3/inputs \\
        [--template AF3_PPI.json] \\
        [--num-seeds 20] \\
        [--skip-existing]

Output:
    One JSON file per binder-target pair, named binder{i}_target{j}.json
    (1-based indices matching the SLURM array task IDs used by cofold_ppi_af3.sb).
"""

import argparse
import copy
import csv
import json
import os
import sys
import zlib


MAX_SEED = 2_147_483_647


def generate_seeds(job_tag: str, n_seeds: int) -> "list[int]":
    """Generate deterministic model seeds from a job tag via CRC32."""
    base_seed = (zlib.crc32(job_tag.encode("utf-8")) % (MAX_SEED - 1)) + 1
    seeds = []
    for k in range(n_seeds):
        s = (base_seed + k) % MAX_SEED
        if s == 0:
            s = 1
        seeds.append(int(s))
    return seeds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate AF3 JSON files for all binder-target combinations."
    )
    parser.add_argument("--csv", required=True, help="Input CSV with 'binders' and 'targets' columns")
    parser.add_argument(
        "--template",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "AF3_PPI.json"),
        help="Path to AF3_PPI.json template (default: AF3_PPI.json next to this script)",
    )
    parser.add_argument("--outdir", required=True, help="Directory to write JSON files into")
    parser.add_argument("--num-seeds", type=int, default=1, help="Number of model seeds per job (default: 1)")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip writing a JSON if the output file already exists",
    )
    return parser.parse_args()


def load_sequences(csv_path: str) -> "tuple[list[str], list[str]]":
    """Return (binders, targets) lists from the CSV, skipping blank entries."""
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            sys.exit(f"ERROR: CSV file is empty: {csv_path}")
        missing = {"binders", "targets"} - set(reader.fieldnames)
        if missing:
            sys.exit(
                f"ERROR: CSV must have 'binders' and 'targets' columns; "
                f"missing: {sorted(missing)}; found: {reader.fieldnames}"
            )
        binders = []
        targets = []
        for row in reader:
            b = row["binders"].strip()
            t = row["targets"].strip()
            if b:
                binders.append(b)
            if t:
                targets.append(t)

    if not binders:
        sys.exit("ERROR: No binder sequences found in CSV.")
    if not targets:
        sys.exit("ERROR: No target sequences found in CSV.")

    # Warn on duplicates within each column
    if len(set(binders)) < len(binders):
        print("WARNING: Duplicate binder sequences detected in CSV.", file=sys.stderr)
    if len(set(targets)) < len(targets):
        print("WARNING: Duplicate target sequences detected in CSV.", file=sys.stderr)

    return binders, targets


def load_template(template_path: str) -> dict:
    if not os.path.isfile(template_path):
        sys.exit(f"ERROR: Template JSON not found: {template_path}")
    with open(template_path) as f:
        data = json.load(f)

    # Validate that the template has exactly two protein entries
    protein_entries = [item["protein"] for item in data.get("sequences", []) if "protein" in item]
    if len(protein_entries) != 2:
        sys.exit(
            f"ERROR: Template must have exactly 2 protein entries; found {len(protein_entries)}. "
            f"Check {template_path}"
        )
    return data


def main() -> None:
    args = parse_args()

    binders, targets = load_sequences(args.csv)
    template = load_template(args.template)

    M, N = len(binders), len(targets)
    total = M * N

    os.makedirs(args.outdir, exist_ok=True)

    generated = 0
    skipped = 0

    for i, binder_seq in enumerate(binders):
        for j, target_seq in enumerate(targets):
            job_tag = f"binder{i + 1}_target{j + 1}"
            out_path = os.path.join(args.outdir, f"{job_tag}.json")

            if args.skip_existing and os.path.isfile(out_path):
                skipped += 1
                continue

            data = copy.deepcopy(template)
            data["name"] = job_tag

            protein_entries = [item["protein"] for item in data["sequences"] if "protein" in item]
            protein_entries[0]["sequence"] = binder_seq   # Chain A = binder
            protein_entries[1]["sequence"] = target_seq   # Chain B = target

            data["modelSeeds"] = generate_seeds(job_tag, args.num_seeds)

            with open(out_path, "w") as f:
                json.dump(data, f, indent=2)

            generated += 1

    print(f"Generated {generated} JSON files ({M} binders x {N} targets).")
    if skipped:
        print(f"Skipped {skipped} already-existing files.")
    print(f"Total combinations: {total}")
    print(f"\nSet SLURM array to: 1-{total}%10")


if __name__ == "__main__":
    main()
