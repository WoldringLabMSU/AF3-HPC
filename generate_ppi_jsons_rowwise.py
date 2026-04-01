#!/usr/bin/env python3
"""
generate_ppi_jsons_rowwise_canonical_ids.py

Generate AlphaFold3 JSON input files from a CSV with 'binders' and 'targets'
columns, preserving ROW-BY-ROW pairing while assigning stable numeric IDs to
unique binder and target sequences based on first appearance.

Example:
    Row 1: binder A, target X -> binder1_target1.json
    Row 2: binder A, target Y -> binder1_target2.json
    Row 3: binder A, target Z -> binder1_target3.json
    Row 4: binder B, target X -> binder2_target1.json
    Row 5: binder B, target Z -> binder2_target3.json

This avoids all-vs-all expansion while giving repeated binders/targets stable
IDs across the whole CSV.

Usage:
    python3 generate_ppi_jsons_rowwise_canonical_ids.py \
        --csv ppi.csv \
        --outdir /path/to/af3/inputs \
        [--template AF3_PPI.json] \
        [--num-seeds 20] \
        [--skip-existing] \
        [--mapping-csv pair_id_mapping.csv]
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
        description="Generate AF3 JSON files from CSV rows using stable binder/target IDs based on first appearance."
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
    parser.add_argument(
        "--mapping-csv",
        default=None,
        help="Optional CSV path to write row-to-binder/target ID mappings",
    )
    return parser.parse_args()


def load_pairs(csv_path: str) -> "list[tuple[int, str, str]]":
    """
    Return ordered list of (source_row_number, binder_seq, target_seq),
    skipping rows with blanks.
    """
    pairs = []

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

        for row_num, row in enumerate(reader, start=2):  # header is row 1
            b = row["binders"].strip()
            t = row["targets"].strip()

            if not b and not t:
                continue
            if not b or not t:
                print(
                    f"WARNING: Skipping row {row_num} because one sequence is blank "
                    f"(binder present={bool(b)}, target present={bool(t)}).",
                    file=sys.stderr,
                )
                continue

            pairs.append((row_num, b, t))

    if not pairs:
        sys.exit("ERROR: No complete binder-target pairs found in CSV.")

    return pairs


def load_template(template_path: str) -> dict:
    if not os.path.isfile(template_path):
        sys.exit(f"ERROR: Template JSON not found: {template_path}")

    with open(template_path) as f:
        data = json.load(f)

    protein_entries = [item["protein"] for item in data.get("sequences", []) if "protein" in item]
    if len(protein_entries) != 2:
        sys.exit(
            f"ERROR: Template must have exactly 2 protein entries; found {len(protein_entries)}. "
            f"Check {template_path}"
        )

    return data


def assign_stable_ids(pairs: "list[tuple[int, str, str]]"):
    """
    Assign binder IDs and target IDs based on first appearance of each unique sequence.
    Returns:
        annotated_pairs: list of dicts with row number, ids, sequences, job_tag
        binder_ids: dict sequence -> binder_id
        target_ids: dict sequence -> target_id
    """
    binder_ids = {}
    target_ids = {}
    annotated_pairs = []

    next_binder_id = 1
    next_target_id = 1

    for row_num, binder_seq, target_seq in pairs:
        if binder_seq not in binder_ids:
            binder_ids[binder_seq] = next_binder_id
            next_binder_id += 1

        if target_seq not in target_ids:
            target_ids[target_seq] = next_target_id
            next_target_id += 1

        binder_id = binder_ids[binder_seq]
        target_id = target_ids[target_seq]
        job_tag = f"binder{binder_id}_target{target_id}"

        annotated_pairs.append({
            "source_row": row_num,
            "binder_id": binder_id,
            "target_id": target_id,
            "binder_sequence": binder_seq,
            "target_sequence": target_seq,
            "job_tag": job_tag,
        })

    return annotated_pairs, binder_ids, target_ids


def main() -> None:
    args = parse_args()

    pairs = load_pairs(args.csv)
    annotated_pairs, binder_ids, target_ids = assign_stable_ids(pairs)
    template = load_template(args.template)

    os.makedirs(args.outdir, exist_ok=True)

    generated = 0
    skipped = 0

    for item in annotated_pairs:
        job_tag = item["job_tag"]
        out_path = os.path.join(args.outdir, f"{job_tag}.json")

        if args.skip_existing and os.path.isfile(out_path):
            skipped += 1
            continue

        data = copy.deepcopy(template)
        data["name"] = job_tag

        protein_entries = [entry["protein"] for entry in data["sequences"] if "protein" in entry]
        protein_entries[0]["sequence"] = item["binder_sequence"]   # Chain A = binder
        protein_entries[1]["sequence"] = item["target_sequence"]   # Chain B = target

        data["modelSeeds"] = generate_seeds(job_tag, args.num_seeds)

        with open(out_path, "w") as f:
            json.dump(data, f, indent=2)

        generated += 1

    if args.mapping_csv:
        with open(args.mapping_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "source_row",
                    "job_tag",
                    "binder_id",
                    "target_id",
                    "binder_sequence",
                    "target_sequence",
                ],
            )
            writer.writeheader()
            for item in annotated_pairs:
                writer.writerow(item)

    print(f"Generated {generated} JSON files from {len(annotated_pairs)} row-wise binder-target pairs.")
    print(f"Unique binders: {len(binder_ids)}")
    print(f"Unique targets: {len(target_ids)}")
    if skipped:
        print(f"Skipped {skipped} already-existing files.")
    if args.mapping_csv:
        print(f"Mapping table written to: {args.mapping_csv}")
    print(f"Total jobs: {len(annotated_pairs)}")
    print(f"\nSet SLURM array to: 1-{len(annotated_pairs)}%10")


if __name__ == "__main__":
    main()
