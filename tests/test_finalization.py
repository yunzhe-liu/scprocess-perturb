from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
FINALIZE = ROOT / "scripts/finalize_perturbation_adata.py"
VALIDATE = ROOT / "scripts/validate_final_adata.py"


def make_input(path: Path) -> None:
    counts = sparse.csr_matrix(np.asarray([[1, 0, 3], [0, 2, 0], [4, 1, 0]], dtype=np.int32))
    totals = np.asarray(counts.sum(axis=1)).ravel()
    normalized = counts.multiply(10000 / totals[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data).astype(np.float32)
    obs = pd.DataFrame(
        {
            "guide_id": ["g1", "g2", "g1;g2"],
            "target_label": ["A", "non-targeting", "A+B"],
            "perturbation_group": ["A", "non-targeting", "A+B"],
            "is_ntc": [False, True, False],
            "batch_id": ["b1", "b1", "b2"],
            "assignment_structure": ["single_guide", "single_guide", "mixed_construct"],
        },
        index=pd.Index(["c1", "c2", "c3"], name="cell_id"),
    )
    adata = ad.AnnData(X=normalized, obs=obs, var=pd.DataFrame(index=["a", "b", "c"]))
    adata.layers["counts"] = counts
    adata.write_h5ad(path, compression="gzip")


class FinalizationTests(unittest.TestCase):
    def run_finalize(self, directory: Path, method: str, status: Path | None = None):
        source = directory / "input.h5ad"
        output = directory / f"final_{method}.h5ad"
        report = directory / f"qc_{method}.json"
        command = [sys.executable, str(FINALIZE), "--input", str(source), "--status-method", method,
                   "--output", str(output), "--report", str(report), "--scan-chunk-values", "2",
                   "--min-free-disk-gb", "0"]
        if status:
            command.extend(["--status-table", str(status)])
        subprocess.run(command, check=True)
        return source, output, report

    def test_none_has_fixed_status_interface_and_preserves_matrices(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            make_input(directory / "input.h5ad")
            source, output, report = self.run_finalize(directory, "none")
            validation = directory / "validation.json"
            subprocess.run([sys.executable, str(VALIDATE), "--source", str(source), "--final", str(output),
                            "--status-method", "none", "--report", str(validation), "--chunk-values", "2"], check=True)
            result = ad.read_h5ad(output)
            self.assertTrue((result.obs["perturbation_status_method"] == "none").all())
            self.assertTrue(result.obs["perturbation_status_score"].isna().all())
            self.assertTrue((result.obs["perturbation_status_reason"] == "not_run").all())
            self.assertFalse(result.obs["perturbation_status_scorable"].any())
            self.assertEqual(json.loads(report.read_text())["status"], "PASS")
            self.assertEqual(json.loads(validation.read_text())["status"], "PASS")

    def test_status_covers_eligible_cells_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            make_input(directory / "input.h5ad")
            status = directory / "status.tsv"
            pd.DataFrame(
                {
                    "cell_id": ["c1", "c2"], "method": ["mixscape", "mixscape"],
                    "target_label": ["A", "non-targeting"],
                    "is_ntc": [False, True],
                    "assignment_structure": ["single_guide", "single_guide"],
                    "score": [0.8, 0.0], "native_status": ["KO", "non-targeting"],
                    "scorable": [True, True], "unscorable_reason": [pd.NA, pd.NA],
                }
            ).to_csv(status, sep="\t", index=False)
            source, output, _ = self.run_finalize(directory, "mixscape", status)
            validation = directory / "validation.json"
            subprocess.run([sys.executable, str(VALIDATE), "--source", str(source), "--final", str(output),
                            "--status-method", "mixscape", "--report", str(validation)], check=True)
            result = ad.read_h5ad(output)
            self.assertEqual(result.obs.loc["c3", "perturbation_status_reason"], "not_eligible")
            self.assertFalse(result.obs.loc["c3", "perturbation_status_scorable"])
            self.assertTrue(result.obs.loc["c1", "perturbation_status_scorable"])

    def test_missing_eligible_status_fails_without_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            make_input(directory / "input.h5ad")
            status = directory / "status.tsv"
            pd.DataFrame(
                {"cell_id": ["c1"], "target_label": ["A"], "is_ntc": [False],
                 "assignment_structure": ["single_guide"], "method": ["ps"],
                 "score": [0.3], "native_status": ["modeled"],
                 "scorable": [True], "unscorable_reason": [pd.NA]}
            ).to_csv(status, sep="\t", index=False)
            output = directory / "final.h5ad"
            process = subprocess.run(
                [sys.executable, str(FINALIZE), "--input", str(directory / "input.h5ad"),
                 "--status-method", "ps", "--status-table", str(status), "--output", str(output),
                 "--report", str(directory / "report.json"), "--min-free-disk-gb", "0"],
                capture_output=True, text=True,
            )
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("Status coverage mismatch", process.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
