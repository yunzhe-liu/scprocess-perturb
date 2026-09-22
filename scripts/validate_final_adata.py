#!/usr/bin/env python3
"""Independent validation of a finalized perturbation-screen H5AD."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


ELIGIBLE = {"single_guide", "concordant_construct"}
STATUS_COLUMNS = [
    "perturbation_status_method", "perturbation_status_score",
    "perturbation_status_label", "perturbation_status_scorable",
    "perturbation_status_reason",
]


def digest_dataset(dataset: h5py.Dataset, chunk_values: int) -> str:
    digest = hashlib.sha256()
    digest.update(str(dataset.shape).encode())
    digest.update(dataset.dtype.str.encode())
    if dataset.ndim == 0:
        digest.update(np.asarray(dataset[()]).tobytes())
        return digest.hexdigest()
    step = max(1, chunk_values // max(1, int(np.prod(dataset.shape[1:]))))
    for start in range(0, dataset.shape[0], step):
        digest.update(np.ascontiguousarray(dataset[start : start + step]).tobytes())
    return digest.hexdigest()


def matrix_digests(handle: h5py.File, path: str, chunk_values: int) -> dict:
    item = handle[path]
    if isinstance(item, h5py.Dataset):
        return {"shape": list(item.shape), "dense": digest_dataset(item, chunk_values)}
    return {
        "shape": list(map(int, item.attrs["shape"])),
        **{name: digest_dataset(item[name], chunk_values) for name in ("data", "indices", "indptr")},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--status-method", choices=("none", "mixscape", "ps"), required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--chunk-values", type=int, default=10_000_000)
    parser.add_argument("--merge-existing", action="store_true")
    args = parser.parse_args()
    with h5py.File(args.source, "r") as source, h5py.File(args.final, "r") as final:
        source_x = matrix_digests(source, "X", args.chunk_values)
        final_x = matrix_digests(final, "X", args.chunk_values)
        source_counts = matrix_digests(source, "layers/counts", args.chunk_values)
        final_counts = matrix_digests(final, "layers/counts", args.chunk_values)
    backed = ad.read_h5ad(args.final, backed="r")
    try:
        obs = backed.obs.copy()
        missing = set(STATUS_COLUMNS).difference(obs.columns)
        if missing:
            raise ValueError(f"Final H5AD lacks status columns: {sorted(missing)}")
        methods = set(obs["perturbation_status_method"].astype(str))
        if methods != {args.status_method}:
            raise ValueError(f"Unexpected final status methods: {sorted(methods)}")
        eligible = obs["assignment_structure"].astype(str).isin(ELIGIBLE)
        scorable = obs["perturbation_status_scorable"].astype(bool)
        reason = obs["perturbation_status_reason"].astype("object")
        if args.status_method == "none":
            semantic_ok = bool((~scorable).all() and (reason == "not_run").all())
        else:
            noneligible_ok = bool((~scorable.loc[~eligible]).all() and (reason.loc[~eligible] == "not_eligible").all())
            eligible_reason_ok = bool(reason.loc[eligible & scorable].isna().all())
            semantic_ok = noneligible_ok and eligible_reason_ok
        if not semantic_ok:
            raise ValueError("Final status field combinations violate the contract")
        if "expression_contract" not in backed.uns or "finalization" not in backed.uns:
            raise ValueError("Finalization provenance is missing from .uns")
        summary = {
            "status": "PASS", "source": str(args.source.resolve()), "final": str(args.final.resolve()),
            "status_method": args.status_method, "shape": list(backed.shape),
            "cells": int(backed.n_obs), "genes": int(backed.n_vars),
            "eligible_cells": int(eligible.sum()), "noneligible_cells": int((~eligible).sum()),
            "scorable_cells": int(scorable.sum()),
            "x_exact": source_x == final_x, "counts_exact": source_counts == final_counts,
            "x_digest": source_x, "counts_digest": source_counts,
            "status_semantics": "PASS", "h5ad_reopen": "PASS",
        }
        if not summary["x_exact"] or not summary["counts_exact"]:
            raise ValueError("X or layers['counts'] changed during finalization")
    finally:
        backed.file.close()
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.report.with_name(args.report.name + f".partial-{os.getpid()}")
        if args.merge_existing and args.report.exists():
            combined = json.loads(args.report.read_text(encoding="utf-8"))
            combined["independent_validation"] = summary
            combined["status"] = "PASS"
        else:
            combined = summary
        temporary.write_text(json.dumps(combined, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, args.report)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
