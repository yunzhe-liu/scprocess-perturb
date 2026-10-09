#!/usr/bin/env python3
"""Write one compact provenance YAML beside the final perturbation H5AD."""

from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import yaml
from anndata.io import read_elem


SCHEMA_VERSION = "1.0"
KEY_PACKAGES = ("anndata", "h5py", "numpy", "pandas", "scipy", "PyYAML")
KEY_ENV_PACKAGES = {
    "python", "snakemake", "simpleaf", "piscem", "alevin-fry", "salmon",
    "r-base", "r-seurat", "r-mixtools", "scmageck", "presto",
}


def _builtin(value: Any) -> Any:
    """Convert AnnData/numpy values into YAML-safe Python values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_builtin(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): _builtin(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_builtin(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _decode_payload(payload: str) -> dict:
    try:
        raw = base64.b64decode(payload.encode("ascii"), validate=True)
        value = json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid base64-encoded JSON payload") from exc
    if not isinstance(value, dict):
        raise ValueError("provenance payload must decode to a JSON object")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path_value: str | Path | None) -> dict | None:
    if not path_value:
        return None
    path = Path(path_value).expanduser()
    record = {"path": str(path.resolve(strict=False)), "exists": path.is_file()}
    if path.is_file():
        stat = path.stat()
        record.update({
            "size_bytes": int(stat.st_size),
            "modified_at": datetime.fromtimestamp(
                stat.st_mtime, tz=timezone.utc
            ).isoformat(),
        })
    return record


def _git_record(project_dir: Path) -> dict:
    def run(*args: str) -> str | None:
        try:
            result = subprocess.run(
                ["git", "-C", str(project_dir), *args], check=True,
                capture_output=True, text=True, timeout=15,
            )
        except (FileNotFoundError, subprocess.SubprocessError):
            return None
        return result.stdout.strip()

    commit = run("rev-parse", "HEAD")
    # Ignore unrelated untracked files; tracked workflow modifications still
    # mark the run as dirty and therefore non-identical to the recorded commit.
    status = run("status", "--porcelain", "--untracked-files=no")
    return {
        "repository": run("config", "--get", "remote.origin.url"),
        "commit": commit,
        "commit_short": commit[:12] if commit else None,
        "dirty": bool(status) if status is not None else None,
    }


def _package_versions() -> dict:
    versions = {"python": platform.python_version()}
    for package in KEY_PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _environment_records(project_dir: Path) -> list[dict]:
    records = []
    env_dir = project_dir / "envs"
    if not env_dir.is_dir():
        return records
    paths = sorted(set(env_dir.glob("*.yaml")) | set(env_dir.glob("*.post-deploy.sh")))
    for path in paths:
        specification = {}
        if path.suffix == ".yaml":
            try:
                specification = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            except (OSError, yaml.YAMLError):
                specification = {}
        key_dependencies = []
        for dependency in specification.get("dependencies", []):
            if not isinstance(dependency, str):
                continue
            name = dependency.split("=", 1)[0].lower()
            if name in KEY_ENV_PACKAGES:
                key_dependencies.append(dependency)
        records.append({
            "file": str(path.relative_to(project_dir)),
            "sha256": _sha256(path),
            "key_dependencies": key_dependencies,
        })
    return records


def _input_records(config: dict, groups: dict) -> dict:
    group_records = []
    for name in sorted(groups):
        group = groups[name]
        group_records.append({
            "group": name,
            "expression": _file_record(group.get("gex_h5")),
            "guide_r1": [
                _file_record(path) for path in group.get("sgRNA_r1", [])
            ],
            "guide_r2": [
                _file_record(path) for path in group.get("sgRNA_r2", [])
            ],
        })
    integration = config.get("integration", {})
    return {
        "groups": group_records,
        "guide_library": _file_record(
            config.get("assignment", {}).get("guide_csv")
            or config.get("guide_csv")
        ),
        "counts_source": _file_record(integration.get("counts_source")),
        "normalized_source": _file_record(integration.get("normalized_source")),
    }


def _read_json(path: Path | None) -> dict | None:
    if path is None:
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read provenance input {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"provenance input {path} must contain a JSON object")
    return value


def _decode_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _matrix_shape(node: h5py.Dataset | h5py.Group) -> list[int]:
    shape = node.attrs.get("shape") if isinstance(node, h5py.Group) else node.shape
    return [int(value) for value in shape]


def _h5ad_summary(path: Path) -> tuple[list[int], list[str], list[str], list[str], dict]:
    wanted_uns = (
        "assignment_score", "integration", "manifest", "standardization",
        "expression_contract", "finalization",
    )
    with h5py.File(path, "r") as handle:
        shape = _matrix_shape(handle["X"])
        obs = handle["obs"]
        order = obs.attrs.get("column-order", [])
        obs_fields = sorted(_decode_text(value) for value in order)
        if not obs_fields:
            index_name = _decode_text(obs.attrs.get("_index", "_index"))
            obs_fields = sorted(key for key in obs.keys() if key != index_name)
        layers = sorted(map(str, handle.get("layers", {}).keys()))
        obsm = sorted(map(str, handle.get("obsm", {}).keys()))
        uns = handle.get("uns")
        metadata = {
            key: _builtin(read_elem(uns[key]))
            for key in wanted_uns
            if uns is not None and key in uns
        }
    return shape, obs_fields, layers, obsm, metadata


def _operation_records(config: dict, status_manifest: dict | None) -> list[dict]:
    extraction = config.get("guide_extraction", {}).get("method", "simpleaf")
    methods = config.get("assignment", {}).get("methods", [])
    assignment = methods[0] if len(methods) == 1 else methods
    status_method = str(
        config.get("perturbation_status", {}).get("method", "none")
    ).lower()
    status_state = "skipped"
    if status_method != "none":
        status_state = "skipped" if (status_manifest or {}).get("skipped") else "success"
    return [
        {"stage": "guide_reference_generation", "status": "success"},
        {"stage": "barcode_whitelist_extraction", "status": "success"},
        {"stage": "guide_umi_quantification", "method": extraction, "status": "success"},
        {"stage": "lane_merge", "status": "success"},
        {"stage": "guide_assignment", "method": assignment, "status": "success"},
        {"stage": "perturbation_integration", "status": "success"},
        {"stage": "perturbation_status_estimation", "method": status_method, "status": status_state},
        {"stage": "data_standardization_and_validation", "status": "success"},
    ]


def build_provenance(
    final_h5ad: Path,
    project_dir: Path,
    config: dict,
    groups: dict,
    status_method: str,
    snakemake_version: str,
    status_manifest: dict | None,
    status_validation: dict | None,
) -> dict:
    shape, obs_fields, layers, obsm, metadata = _h5ad_summary(final_h5ad)

    stat = final_h5ad.stat()
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "workflow": {
            "name": "scprocess-perturb",
            "engine": "Snakemake",
            "engine_version": snakemake_version,
            **_git_record(project_dir),
            "status": "success",
        },
        "software": {
            "provenance_runtime": _package_versions(),
            "declared_environments": _environment_records(project_dir),
        },
        "inputs": _input_records(config, groups),
        "parameters": _builtin(config),
        "operations": _operation_records(config, status_manifest),
        "matrix_and_annotations": {
            "shape": shape,
            "x": metadata.get("expression_contract", {}).get("representation"),
            "layers": layers,
            "obsm": obsm,
            "obs_fields": obs_fields,
            "metadata": metadata,
        },
        "validation": {
            "final_h5ad": "PASS",
            "status_method": status_method,
            "perturbation_status_manifest": _builtin(status_manifest),
            "perturbation_status_validation": _builtin(status_validation),
        },
        "output": {
            "path": str(final_h5ad.resolve()),
            "size_bytes": int(stat.st_size),
            "modified_at": datetime.fromtimestamp(
                stat.st_mtime, tz=timezone.utc
            ).isoformat(),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-h5ad", type=Path, required=True)
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--config-b64", required=True)
    parser.add_argument("--groups-b64", required=True)
    parser.add_argument("--status-method", required=True)
    parser.add_argument("--snakemake-version", default="unknown")
    parser.add_argument("--status-manifest", type=Path)
    parser.add_argument("--status-validation", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not args.final_h5ad.is_file():
        raise FileNotFoundError(f"final H5AD does not exist: {args.final_h5ad}")
    config = _decode_payload(args.config_b64)
    groups = _decode_payload(args.groups_b64)
    document = build_provenance(
        args.final_h5ad,
        args.project_dir.resolve(),
        config,
        groups,
        args.status_method,
        args.snakemake_version,
        _read_json(args.status_manifest),
        _read_json(args.status_validation),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + f".partial-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(document, handle, sort_keys=False, allow_unicode=True)
        # Verify that the complete temporary document is parseable before publish.
        verified = yaml.safe_load(temporary.read_text(encoding="utf-8"))
        if not isinstance(verified, dict) or verified.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("generated provenance YAML failed validation")
        os.replace(temporary, args.output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"[done] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
