from __future__ import annotations

import base64
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
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "write_provenance.py"


def encode(value: dict) -> str:
    return base64.b64encode(json.dumps(value).encode("utf-8")).decode("ascii")


class ProvenanceTests(unittest.TestCase):
    def test_compact_provenance_is_complete_and_parseable(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            expression = directory / "expression.h5ad"
            final = directory / "perturbation_adata.h5ad"
            output = directory / "provenance.yaml"

            counts = sparse.csr_matrix(np.asarray([[1, 0], [0, 2]], dtype=np.int32))
            adata = ad.AnnData(
                X=counts.astype(np.float32),
                obs=pd.DataFrame(index=["c1", "c2"]),
                var=pd.DataFrame(index=["g1", "g2"]),
            )
            adata.layers["counts"] = counts
            adata.uns["integration"] = {
                "cell_id_matching_selected": "exact",
                "alignment": {"matched_cells": 2},
            }
            adata.uns["standardization"] = {
                "normalization_method": "log1p_cp10k",
            }
            adata.uns["expression_contract"] = {
                "representation": "normalized log-transformed expression",
                "counts_layer": "counts",
            }
            adata.uns["finalization"] = {
                "status_method": "none",
                "cell_filtering": "not_performed",
            }
            adata.write_h5ad(expression)
            adata.write_h5ad(final)

            config = {
                "proj_dir": str(ROOT),
                "guide_csv": str(directory / "library.csv"),
                "guide_extraction": {"method": "simpleaf"},
                "assignment": {"guide_design": "dual", "methods": ["pgmm_em"]},
                "perturbation_status": {"method": "none"},
            }
            groups = {
                "lane_01": {
                    "gex_h5": str(expression),
                    "sgRNA_r1": [],
                    "sgRNA_r2": [],
                }
            }
            subprocess.run([
                sys.executable, str(SCRIPT),
                "--final-h5ad", str(final),
                "--project-dir", str(ROOT),
                "--config-b64", encode(config),
                "--groups-b64", encode(groups),
                "--status-method", "none",
                "--snakemake-version", "7.32.4",
                "--output", str(output),
            ], check=True)

            document = yaml.safe_load(output.read_text(encoding="utf-8"))
            self.assertEqual(document["schema_version"], "1.0")
            self.assertEqual(document["workflow"]["status"], "success")
            self.assertEqual(document["workflow"]["engine_version"], "7.32.4")
            self.assertEqual(document["parameters"]["assignment"]["methods"], ["pgmm_em"])
            self.assertEqual(document["matrix_and_annotations"]["shape"], [2, 2])
            self.assertEqual(document["validation"]["final_h5ad"], "PASS")
            self.assertEqual(document["output"]["size_bytes"], final.stat().st_size)
            self.assertFalse(list(directory.glob("provenance.yaml.partial-*")))


if __name__ == "__main__":
    unittest.main()
