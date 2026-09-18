#!/usr/bin/env python3
"""Batch SpliceAI scoring for SPRED1 NM_152594.3 variants."""

import argparse
import csv
import os
import time
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import genebe as gnb
import requests

TARGET_GENE = "SPRED1"
TARGET_TRANSCRIPT = "ENST00000299084.9"
TARGET_REFSEQ = "NM_152594.3"

GENOME = "hg38"
HG = "38"
GENCODE_SET = "basic"
DISTANCE = 500
MASK = 0

API = os.environ.get("SPLICEAI_API", "http://localhost:8080/spliceai/")

FIELDS = [
    "HGVS",
    "resolved_hg38",
    "transcript",
    "transcript_priority",
    "gene",
    "DS_AG",
    "DS_AL",
    "DS_DG",
    "DS_DL",
    "SpliceAI_max",
    "max_type",
    "status",
]

DS_KEYS = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")


def format_2dp(value):
    """Format with conventional half-up rounding: 0.005 -> 0.01."""
    d = Decimal(str(value))
    return format(d.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP), ".2f")


def empty_row(hgvs, resolved="", status="ERROR"):
    return {
        "HGVS": hgvs,
        "resolved_hg38": resolved,
        "transcript": "",
        "transcript_priority": "",
        "gene": "",
        "DS_AG": "",
        "DS_AL": "",
        "DS_DG": "",
        "DS_DL": "",
        "SpliceAI_max": "",
        "max_type": "",
        "status": status,
    }


def read_variants(path):
    """Read unique non-empty HGVS records while preserving input order."""
    variants = []
    seen = set()
    with open(path) as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#") or value in seen:
                continue
            variants.append(value)
            seen.add(value)
    return variants


def read_completed(output_path, input_variants):
    """Load prior OK rows for --resume."""
    path = Path(output_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}

    completed = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise RuntimeError(
                "Existing output has unexpected columns; move or rename it "
                "before using --resume."
            )
        for row in reader:
            hgvs = (row.get("HGVS") or "").strip()
            if hgvs and row.get("status") == "OK":
                completed[hgvs] = row

    return {v: completed[v] for v in input_variants if v in completed}


def resolve_hgvs(variants):
    """Convert transcript HGVS variants to GRCh38 chrom-pos-ref-alt."""
    if not variants:
        return {}

    print(f"Resolving {len(variants)} variants with GeneBe...", flush=True)
    resolved = gnb.parse_variants(variants, genome=GENOME)

    if len(resolved) != len(variants):
        raise RuntimeError(
            f"GeneBe returned {len(resolved)} results for {len(variants)} inputs."
        )

    return dict(zip(variants, resolved))


def query_spliceai(session, genomic_variant, retries=4):
    params = {
        "hg": HG,
        "bc": GENCODE_SET,
        "distance": DISTANCE,
        "mask": MASK,
        "variant": genomic_variant,
    }

    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = session.get(API, params=params, timeout=900)
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(5 * attempt)

    raise last_error


def target_record(scores):
    """Return only the SPRED1 MANE Select transcript used for this project."""
    for record in scores or []:
        if (
            record.get("g_name") == TARGET_GENE
            and record.get("t_id") == TARGET_TRANSCRIPT
        ):
            return record
    return None


def score_variant(session, hgvs, resolved):
    if not resolved:
        return empty_row(hgvs, status="ERROR: HGVS_TO_HG38_FAILED")

    data = query_spliceai(session, resolved)
    target = target_record(data.get("scores", []))

    if target is None:
        return empty_row(
            hgvs,
            resolved,
            f"ERROR: TARGET_TRANSCRIPT_NOT_FOUND ({TARGET_TRANSCRIPT})",
        )

    raw = {
        key: Decimal(str(target.get(key, "0") or "0"))
        for key in DS_KEYS
    }

    max_type = max(DS_KEYS, key=lambda key: raw[key])
    max_score = raw[max_type]

    return {
        "HGVS": hgvs,
        "resolved_hg38": resolved,
        "transcript": target.get("t_id", ""),
        "transcript_priority": target.get("t_priority", ""),
        "gene": target.get("g_name", ""),
        "DS_AG": format_2dp(raw["DS_AG"]),
        "DS_AL": format_2dp(raw["DS_AL"]),
        "DS_DG": format_2dp(raw["DS_DG"]),
        "DS_DL": format_2dp(raw["DS_DL"]),
        "SpliceAI_max": format_2dp(max_score),
        "max_type": max_type,
        "status": "OK",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Batch-score SPRED1 variants using a local SpliceAI API."
    )
    parser.add_argument("input", help="Text file with one HGVS variant per line")
    parser.add_argument("output", help="Output CSV")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Keep prior OK rows and retry missing/error variants",
    )
    args = parser.parse_args()

    variants = read_variants(args.input)
    if not variants:
        raise SystemExit("No variants found in input file.")

    previous = read_completed(args.output, variants) if args.resume else {}
    remaining = [v for v in variants if v not in previous]

    print(f"Loaded {len(variants)} HGVS variants", flush=True)
    print(
        f"Settings: hg38, GENCODE={GENCODE_SET}, distance={DISTANCE}, "
        f"mask={MASK}, transcript={TARGET_TRANSCRIPT}",
        flush=True,
    )
    print(f"Variants to score: {len(remaining)}", flush=True)

    resolved = resolve_hgvs(remaining)
    output = Path(args.output)

    # Rewrite a clean file containing prior successful rows first.
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for variant in variants:
            if variant in previous:
                writer.writerow(previous[variant])
        handle.flush()

    session = requests.Session()

    # Append each newly completed row immediately.
    with output.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)

        for i, hgvs in enumerate(remaining, start=1):
            try:
                row = score_variant(session, hgvs, resolved.get(hgvs))
            except Exception as exc:
                row = empty_row(
                    hgvs,
                    resolved.get(hgvs, ""),
                    f"ERROR: {type(exc).__name__}: {exc}",
                )

            writer.writerow(row)
            handle.flush()

            if row["status"] == "OK":
                print(
                    f"[{i}/{len(remaining)}] {hgvs} -> "
                    f"{row['SpliceAI_max']} {row['max_type']}",
                    flush=True,
                )
            else:
                print(
                    f"[{i}/{len(remaining)}] {hgvs} -> {row['status']}",
                    flush=True,
                )

    # Restore original input order and eliminate duplicates.
    latest = {}
    with output.open(newline="") as handle:
        for row in csv.DictReader(handle):
            latest[row["HGVS"]] = row

    temp = output.with_suffix(output.suffix + ".tmp")
    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for hgvs in variants:
            writer.writerow(
                latest.get(hgvs, empty_row(hgvs, status="ERROR: NO_RESULT_ROW"))
            )

    temp.replace(output)

    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    ok = sum(row["status"] == "OK" for row in rows)
    errors = len(rows) - ok

    print(f"Finished: {ok} OK, {errors} errors", flush=True)


if __name__ == "__main__":
    main()
