#!/usr/bin/env python3
"""
pdb_to_csv.py

Fetch protein sequences for a list of PDB IDs from RCSB and write a
binders,targets CSV compatible with generate_ppi_jsons.py.

Each PDB ID is expected to contain exactly two protein chains. PDB IDs
with a different number of chains are skipped with a warning.

Usage:
    python3 pdb_to_csv.py --pdb-list pdb_ids.txt --output-csv pdb_pairs.csv

Input file format (pdb_ids.txt) — one PDB ID per line, blank lines and
lines starting with '#' are ignored:
    1ABC
    2XYZ
    # this is a comment
    4DEF
"""

import argparse
import csv
import sys
import time
import urllib.error
import urllib.request


def fetch_fasta(pdb_id: str, retries: int = 3) -> str:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.URLError as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                raise RuntimeError(f"Failed to fetch FASTA for {pdb_id}: {e}")


def parse_fasta(fasta_text: str) -> "list[tuple[str, str]]":
    """Return list of (header, sequence) tuples from a multi-FASTA string."""
    entries = []
    current_header = None
    current_seq = []
    for line in fasta_text.strip().splitlines():
        if line.startswith(">"):
            if current_header is not None:
                entries.append((current_header, "".join(current_seq)))
                current_seq = []
            current_header = line[1:].strip()
        else:
            current_seq.append(line.strip())
    if current_header is not None and current_seq:
        entries.append((current_header, "".join(current_seq)))
    return entries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch PDB sequences from RCSB and write a binders/targets CSV."
    )
    parser.add_argument(
        "--pdb-list", required=True,
        help="Text file with PDB IDs, one per line"
    )
    parser.add_argument(
        "--output-csv", required=True,
        help="Output CSV path (binders, targets columns + pdb_id for reference)"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with open(args.pdb_list) as f:
        pdb_ids = [
            line.strip().upper()
            for line in f
            if line.strip() and not line.startswith("#")
        ]

    if not pdb_ids:
        sys.exit("ERROR: No PDB IDs found in input file.")

    print(f"Fetching sequences for {len(pdb_ids)} PDB ID(s)...")

    rows = []
    skipped = []

    for pdb_id in pdb_ids:
        try:
            fasta = fetch_fasta(pdb_id)
            entries = parse_fasta(fasta)
        except RuntimeError as e:
            print(f"  WARNING: {e}", file=sys.stderr)
            skipped.append(pdb_id)
            continue

        if len(entries) != 2:
            print(
                f"  WARNING: {pdb_id} returned {len(entries)} chain(s), expected 2. Skipping.",
                file=sys.stderr,
            )
            skipped.append(pdb_id)
            continue

        binder_seq = entries[0][1]
        target_seq = entries[1][1]
        rows.append({"pdb_id": pdb_id, "binders": binder_seq, "targets": target_seq})
        print(f"  {pdb_id}: chain 1 ({len(binder_seq)} aa)  chain 2 ({len(target_seq)} aa)")

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pdb_id", "binders", "targets"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} pair(s) to {args.output_csv}.")
    if skipped:
        print(f"Skipped {len(skipped)} PDB ID(s): {', '.join(skipped)}")


if __name__ == "__main__":
    main()
