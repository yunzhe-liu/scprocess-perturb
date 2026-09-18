#!/usr/bin/env python3
import csv
import gzip
import json
import os
import sys
from pathlib import Path


def atomic_replace(writer, destination):
    temporary = destination.with_name(destination.name + f".tmp-{os.getpid()}")
    try:
        writer(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def unique_columns(columns):
    """Preserve first occurrence of each metadata field for a stable TSV contract."""
    result = []
    seen = set()
    for column in columns:
        if column not in seen:
            result.append(column)
            seen.add(column)
    return result


def main(root_text, metadata_text):
    root = Path(root_text)
    metadata_path = Path(metadata_text)
    output_dir = root / "merged"
    output_dir.mkdir(parents=True, exist_ok=True)
    shard_paths = sorted((root / "outputs").glob("shard_*/mixscape_metadata.tsv.gz"))
    if not shard_paths:
        raise RuntimeError("No shard metadata files were found")
    missing_files = [str(path) for path in shard_paths if not path.is_file()]
    if missing_files:
        raise RuntimeError("Missing shard metadata: " + ", ".join(missing_files))

    result_by_cell = {}
    ntc_reference = {}
    target_duplicates = []
    result_columns = None
    for shard_index, path in enumerate(shard_paths, start=1):
        with gzip.open(path, "rt", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames is None or reader.fieldnames[0] != "":
                raise RuntimeError(f"Unexpected shard header: {path}")
            columns = [name for name in reader.fieldnames if name.startswith("mixscape_")]
            if not columns:
                raise RuntimeError(f"No Mixscape result columns: {path}")
            if result_columns is None:
                result_columns = columns
            elif columns != result_columns:
                raise RuntimeError(f"Inconsistent result columns: {path}")
            for row in reader:
                cell_id = row[""]
                values = tuple(row[column] for column in result_columns)
                is_ntc = row["is_ntc"].strip().lower() == "true"
                if is_ntc:
                    old = ntc_reference.setdefault(cell_id, values)
                    if old != values:
                        raise RuntimeError(f"Inconsistent NTC result across shards: {cell_id}")
                    if shard_index == 1:
                        result_by_cell[cell_id] = values
                elif cell_id in result_by_cell:
                    target_duplicates.append(cell_id)
                else:
                    result_by_cell[cell_id] = values

    expected_cells = []
    expected_ntc = 0
    with gzip.open(metadata_path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames is None or "cell_id" not in reader.fieldnames or "is_ntc" not in reader.fieldnames:
            raise RuntimeError("Canonical metadata lacks cell_id or is_ntc")
        metadata_columns = unique_columns(reader.fieldnames)
        for row in reader:
            expected_cells.append(row["cell_id"])
            expected_ntc += row["is_ntc"].strip().lower() == "true"

    expected_set = set(expected_cells)
    result_set = set(result_by_cell)
    missing_cells = expected_set - result_set
    unexpected_cells = result_set - expected_set
    if target_duplicates or missing_cells or unexpected_cells or len(expected_cells) != len(expected_set):
        raise RuntimeError("Merge coverage validation failed")

    output_path = output_dir / "mixscape_metadata.tsv.gz"
    def write_output(temporary):
        with gzip.open(temporary, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=metadata_columns + result_columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            with gzip.open(metadata_path, "rt", newline="") as source:
                for row in csv.DictReader(source, delimiter="\t"):
                    values = result_by_cell[row["cell_id"]]
                    row.update(dict(zip(result_columns, values)))
                    writer.writerow(row)
    atomic_replace(write_output, output_path)

    report = {
        "status": "complete",
        "n_shards": len(shard_paths),
        "expected_cells": len(expected_cells),
        "expected_ntc": expected_ntc,
        "merged_cells": len(result_by_cell),
        "result_columns": result_columns,
        "target_duplicate_count": len(target_duplicates),
        "missing_cell_count": len(missing_cells),
        "unexpected_cell_count": len(unexpected_cells),
        "ntc_agreement": "verified",
        "output": str(output_path),
    }
    report_path = output_dir / "merge_validation.json"
    atomic_replace(lambda temporary: temporary.write_text(json.dumps(report, indent=2) + "\n"), report_path)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: merge_pooled_shards.py <shard_root> <canonical_metadata.tsv.gz>")
    main(sys.argv[1], sys.argv[2])
