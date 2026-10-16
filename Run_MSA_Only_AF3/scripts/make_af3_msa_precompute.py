#!/usr/bin/env python3
"""
AUTHOR: BEN DOLGIKH
Edited by: Annie Needs

USAGE: python make_af3_msa_precompute_jsons.py --seq-list-csv ../inputs/{input csv file}  [optional: ---out-dir, -name-col, --seq-col, --chain-id, --model-seed, --template-mode, --dedupe-by-sequence, --error-on-name-conflict, --skip-existing, --clean-out-dir]


Generate AlphaFold 3 monomer JSON files for a data-pipeline-only / MSA-only
precompute run from a sequence-list CSV.

Intended workflow:
  1. Provide a *_seq_list.csv containing one row per protein sequence.
  2. Generate one AF3 JSON per unique sequence, with unpairedMsa=null,
     pairedMsa=null, and templates=[].
  3. Run AF3 with --run_data_pipeline=true --run_inference=false.
  4. Use the resulting *_data.json files as an MSA cache for downstream AF3 predictions

Default expected CSV columns:
  name,sequence

Useful alternatives:
  --name-col <col>
  --seq-col <col>

Outputs:
  <out-dir>/*.json
  <out-dir>/msa_precompute_manifest.tsv
  <out-dir>/msa_precompute_run_list.txt
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUOJ*-]+$", re.IGNORECASE)


def sniff_delimiter(path: Path, default: str = ",") -> str:
    sample = path.read_text(errors="replace")[:8192]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return default


def sanitize_name(value: Any) -> str:
    s = str(value).strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    s = s.strip("._-")
    return s or "unnamed"


def clean_sequence(seq: Any) -> str:
    s = str(seq).strip().replace(" ", "").replace("\n", "").replace("\r", "").upper()
    if not s:
        raise ValueError("empty sequence")
    if not AA_RE.match(s):
        raise ValueError(f"sequence contains unexpected characters: {s[:80]}...")
    return s.replace("*", "")


def seq_hash(seq: str, n: int = 12) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:n]


def infer_columns(fieldnames: List[str], name_col: Optional[str], seq_col: Optional[str]) -> Tuple[str, str]:
    fields = set(fieldnames)

    if name_col:
        if name_col not in fields:
            raise SystemExit(f"ERROR: --name-col {name_col!r} not found. Observed columns: {fieldnames}")
        chosen_name = name_col
    else:
        name_candidates = [
            "name", "protein_name", "seq_name", "sequence_name", "id", "protein_id",
            "query_header", "partner_gene", "gene", "label", "description",
        ]
        chosen_name = next((c for c in name_candidates if c in fields), None)
        if chosen_name is None:
            raise SystemExit(
                "ERROR: could not infer name column. Provide --name-col. "
                f"Observed columns: {fieldnames}"
            )

    if seq_col:
        if seq_col not in fields:
            raise SystemExit(f"ERROR: --seq-col {seq_col!r} not found. Observed columns: {fieldnames}")
        chosen_seq = seq_col
    else:
        seq_candidates = [
            "sequence", "protein_sequence", "seq", "aa_sequence", "protein_seq",
            "binders", "targets", "partner_sequence",
        ]
        chosen_seq = next((c for c in seq_candidates if c in fields), None)
        if chosen_seq is None:
            raise SystemExit(
                "ERROR: could not infer sequence column. Provide --seq-col. "
                f"Observed columns: {fieldnames}"
            )

    return chosen_name, chosen_seq


def read_sequence_rows(path: Path, name_col: Optional[str], seq_col: Optional[str]) -> List[Dict[str, str]]:
    delim = sniff_delimiter(path)
    rows: List[Dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: input CSV has no header: {path}")
        chosen_name, chosen_seq = infer_columns(reader.fieldnames, name_col, seq_col)
        for i, row in enumerate(reader, start=1):
            raw_name = (row.get(chosen_name) or "").strip()
            raw_seq = row.get(chosen_seq)
            if not raw_name and not raw_seq:
                continue
            if not raw_name:
                raw_name = f"row{i}"
            try:
                seq = clean_sequence(raw_seq)
            except Exception as exc:
                raise SystemExit(f"ERROR: bad sequence in row {i}: {exc}") from exc
            rows.append({
                "row_index": str(i),
                "input_name": raw_name,
                "sequence": seq,
                "length": str(len(seq)),
                "sequence_hash": seq_hash(seq),
            })
    if not rows:
        raise SystemExit(f"ERROR: no usable rows found in {path}")
    return rows


def uniquify_job_name(base: str, used: set[str]) -> str:
    base = sanitize_name(base)
    if base not in used:
        used.add(base)
        return base
    k = 2
    while f"{base}_{k}" in used:
        k += 1
    out = f"{base}_{k}"
    used.add(out)
    return out


def build_jobs(rows: List[Dict[str, str]], dedupe_by_sequence: bool, error_on_name_conflict: bool) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Return unique jobs and manifest rows mapping every input row to a job."""
    used_names: set[str] = set()
    name_to_seq: Dict[str, str] = {}
    seq_to_job: OrderedDict[str, Dict[str, str]] = OrderedDict()
    jobs: List[Dict[str, str]] = []
    manifest_rows: List[Dict[str, str]] = []

    for row in rows:
        clean_name = sanitize_name(row["input_name"])
        sequence = row["sequence"]
        sh = row["sequence_hash"]

        if clean_name in name_to_seq and name_to_seq[clean_name] != sequence:
            msg = f"name {clean_name!r} is associated with multiple different sequences"
            if error_on_name_conflict:
                raise SystemExit(f"ERROR: {msg}")
            clean_name = f"{clean_name}_{sh}"
        name_to_seq[clean_name] = sequence

        if dedupe_by_sequence and sequence in seq_to_job:
            job = seq_to_job[sequence]
        else:
            job_name = uniquify_job_name(clean_name, used_names)
            job = {
                "job_name": job_name,
                "canonical_input_name": row["input_name"],
                "sequence": sequence,
                "length": str(len(sequence)),
                "sequence_hash": sh,
            }
            jobs.append(job)
            seq_to_job[sequence] = job

        manifest_rows.append({
            "input_row_index": row["row_index"],
            "input_name": row["input_name"],
            "job_name": job["job_name"],
            "json_file": f"{job['job_name']}.json",
            "sequence_length": str(len(sequence)),
            "sequence_hash": sh,
            "deduped_to_existing_sequence": "1" if job["canonical_input_name"] != row["input_name"] or job["sequence"] == sequence and job["job_name"] != sanitize_name(row["input_name"]) else "0",
        })

    return jobs, manifest_rows


def make_af3_json(job: Dict[str, str], chain_id: str, model_seed: int, template_mode: str) -> Dict[str, Any]:
    protein = {
        "id": chain_id,
        "sequence": job["sequence"],
        "unpairedMsa": None,
        "pairedMsa": None,
    }

    if template_mode == "none":
        # Explicit template-free data pipeline: build MSA, do not search/use templates.
        protein["templates"] = []
    elif template_mode == "search":
        # Omit templates so AF3 may run template search during the data pipeline.
        pass
    else:
        raise ValueError(f"unknown template_mode: {template_mode}")

    return {
        "name": job["job_name"],
        "modelSeeds": [int(model_seed)],
        "sequences": [{"protein": protein}],
        "dialect": "alphafold3",
        "version": 1,
    }


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate AF3 monomer JSONs for data-pipeline-only MSA precomputation."
    )
    parser.add_argument("--seq-list-csv", required=True, help="Input *_seq_list.csv.")
    parser.add_argument("--out-dir", defualt="../json_dir", help="Directory to write AF3 JSON files.")
    parser.add_argument("--name-col", default=None, help="Name column. Default: auto-detect.")
    parser.add_argument("--seq-col", default=None, help="Sequence column. Default: auto-detect.")
    parser.add_argument("--chain-id", default="A", help="Protein chain ID in monomer JSON. Default: A.")
    parser.add_argument("--model-seed", type=int, default=1, help="Single model seed stored in JSON. Inference is not run.")
    parser.add_argument(
        "--template-mode",
        choices=["none", "search"],
        default="none",
        help="none: set templates=[] and only build MSA. search: omit templates and allow AF3 template search.",
    )
    parser.add_argument(
        "--dedupe-by-sequence",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write one JSON per unique sequence. Default: true.",
    )
    parser.add_argument("--error-on-name-conflict", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--clean-out-dir", action="store_true", help="Delete existing *.json in --out-dir before writing.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seq_list = Path(args.seq_list_csv).resolve()
    out_dir = Path(args.out_dir).resolve()

    if not seq_list.is_file():
        raise SystemExit(f"ERROR: input CSV not found: {seq_list}")

    out_dir.mkdir(parents=True, exist_ok=True)
    if args.clean_out_dir:
        for p in out_dir.glob("*.json"):
            p.unlink()

    rows = read_sequence_rows(seq_list, args.name_col, args.seq_col)
    jobs, manifest_rows = build_jobs(
        rows,
        dedupe_by_sequence=args.dedupe_by_sequence,
        error_on_name_conflict=args.error_on_name_conflict,
    )

    written = 0
    skipped = 0
    run_list_paths: List[str] = []

    for job in jobs:
        out_path = out_dir / f"{job['job_name']}.json"
        run_list_paths.append(str(out_path))
        if args.skip_existing and out_path.exists():
            skipped += 1
            continue
        data = make_af3_json(job, args.chain_id, args.model_seed, args.template_mode)
        out_path.write_text(json.dumps(data, indent=2) + "\n")
        written += 1

    # Add absolute JSON path to manifest after writing.
    for row in manifest_rows:
        row["json_path"] = str(out_dir / row["json_file"])

    manifest_path = out_dir / "msa_precompute_manifest.tsv"
    write_tsv(
        manifest_path,
        manifest_rows,
        [
            "input_row_index", "input_name", "job_name", "json_file", "json_path",
            "sequence_length", "sequence_hash", "deduped_to_existing_sequence",
        ],
    )

    run_list = out_dir / "msa_precompute_run_list.txt"
    run_list.write_text("\n".join(sorted(run_list_paths)) + "\n")

    print("MSA precompute JSON generation complete.")
    print(f"Input rows: {len(rows)}")
    print(f"Unique jobs: {len(jobs)}")
    print(f"JSONs written: {written}")
    print(f"JSONs skipped: {skipped}")
    print(f"JSON dir: {out_dir}")
    print(f"Manifest: {manifest_path}")
    print(f"Run list: {run_list}")
    print("\nSubmit array as:")
    print(f"  sbatch --array=1-{len(jobs)}%<CONCURRENCY> cofold_ppi_af3_msa_only.sb")


if __name__ == "__main__":
    main()