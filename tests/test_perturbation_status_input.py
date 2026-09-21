from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import io, sparse


ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT / "scripts" / "perturbation_status" / "adapt_input.py"


def make_input(path: Path, counts: np.ndarray, features: list[str]) -> None:
    matrix = sparse.csr_matrix(counts, dtype=np.int32)
    obs = pd.DataFrame(
        {
            "guide_id": ["g1", "g2"],
            "target_label": ["A", "non-targeting"],
            "perturbation_group": ["A", "non-targeting"],
            "is_ntc": [False, True],
            "batch_id": ["b1", "b1"],
            "assignment_structure": ["single_guide", "single_guide"],
        },
        index=pd.Index(["c1", "c2"], name="cell_id"),
    )
    adata = ad.AnnData(
        X=matrix.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(features, name="feature_id")),
    )
    adata.layers["counts"] = matrix
    adata.write_h5ad(path)


def run_adapter(directory: Path, source: Path, mode: str = "auto") -> subprocess.CompletedProcess:
    output = directory / "output"
    config = {
        "dataset": "test",
        "source_h5ad": str(source),
        "input_dir": str(output),
        "row_chunk_size": 1,
        "feature_mode": mode,
        "selection": {
            "mode": "eligible",
            "eligible_assignment_structures": ["single_guide", "concordant_construct"],
            "max_cells": 0,
        },
    }
    config_path = directory / "adapter.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(ADAPTER), "--config", str(config_path)],
        capture_output=True,
        text=True,
    )


class PerturbationStatusInputTests(unittest.TestCase):
    def test_auto_aggregates_strict_usa_as_s_plus_a(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input.h5ad"
            make_input(
                source,
                np.asarray([[1, 2, 10, 20, 3, 0], [0, 1, 30, 40, 2, 4]]),
                ["g1_S", "g2_S", "g1_U", "g2_U", "g1_A", "g2_A"],
            )
            process = run_adapter(directory, source)
            self.assertEqual(process.returncode, 0, process.stderr)
            with gzip.open(directory / "output" / "counts.mtx.gz", "rb") as stream:
                result = io.mmread(stream).toarray()
            np.testing.assert_array_equal(result, np.asarray([[4, 2], [2, 5]]))
            features = pd.read_csv(
                directory / "output" / "features.tsv.gz", sep="\t", header=None
            )
            self.assertEqual(features[0].tolist(), ["g1", "g2"])
            manifest = json.loads(
                (directory / "output" / "adapter_manifest.json").read_text()
            )
            self.assertEqual(manifest["feature_mode_resolved"], "usa_sa")
            self.assertEqual(manifest["output_shape"], [2, 2])

    def test_auto_preserves_gene_level_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input.h5ad"
            counts = np.asarray([[1, 2], [3, 4]])
            make_input(source, counts, ["g1", "g2"])
            process = run_adapter(directory, source)
            self.assertEqual(process.returncode, 0, process.stderr)
            with gzip.open(directory / "output" / "counts.mtx.gz", "rb") as stream:
                result = io.mmread(stream).toarray()
            np.testing.assert_array_equal(result, counts.T)
            manifest = json.loads(
                (directory / "output" / "adapter_manifest.json").read_text()
            )
            self.assertEqual(manifest["feature_mode_resolved"], "gene")

    def test_auto_rejects_malformed_usa_like_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input.h5ad"
            make_input(
                source,
                np.asarray([[1, 2, 3, 4], [5, 6, 7, 8]]),
                ["g1_S", "g2_S", "g1_U", "g2_A"],
            )
            process = run_adapter(directory, source)
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("fails the strict S/U/A block contract", process.stderr)


if __name__ == "__main__":
    unittest.main()
