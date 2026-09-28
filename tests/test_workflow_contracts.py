from __future__ import annotations

import gzip
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
sys.path.insert(0, str(ROOT / "scripts"))

from integrate_multimodal import _h5ad_safe_frame


class WorkflowContractTests(unittest.TestCase):
    def test_h5ad_safe_frame_removes_string_backed_categories(self):
        frame = pd.DataFrame({
            "assigned_guide": pd.Series(["g1", "g2"], dtype="string")
        })
        safe = _h5ad_safe_frame(frame)
        self.assertIsInstance(safe["assigned_guide"].dtype, pd.CategoricalDtype)
        self.assertEqual(safe["assigned_guide"].cat.categories.dtype, object)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "safe.h5ad"
            ad.AnnData(
                X=sparse.eye(2, format="csr"), obs=safe, var=safe.copy()
            ).write_h5ad(output)
            self.assertEqual(ad.read_h5ad(output).n_obs, 2)

    def test_plain_barcode_translation_uses_distinct_backup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            barcodes = root / "barcodes.txt"
            table = root / "translation.txt"
            barcodes.write_text("AAAAAAAAAAAAAAAA\n", encoding="utf-8")
            table.write_text("AAAAAAAAAAAAAAAA\tCCCCCCCCCCCCCCCC\n", encoding="utf-8")
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "translate_barcodes.py"),
                str(barcodes), "--trans-table", str(table),
            ], check=True)
            self.assertEqual(barcodes.read_text(), "CCCCCCCCCCCCCCCC\n")
            self.assertEqual(
                (root / "barcodes.txt_from_backup").read_text(),
                "AAAAAAAAAAAAAAAA\n",
            )

    def test_long_guide_library_builds_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            library = root / "library.csv"
            fasta = root / "guides.fasta"
            t2g = root / "t2g.tsv"
            pd.DataFrame({
                "construct_id": ["c1", "c1"],
                "guide_id": ["g1", "g2"],
                "guide_sequence": ["ACGT", "TGCA"],
                "target_label": ["GENE1", "GENE1"],
                "guide_position": [1, 2],
            }).to_csv(library, index=False)
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "feature_reference_adapter.py"),
                "--csv", str(library), "--out-fasta", str(fasta),
                "--out-t2g", str(t2g),
            ], check=True)
            self.assertEqual(fasta.read_text(), ">g1\nACGT\n>g2\nTGCA\n")
            self.assertEqual(t2g.read_text(), "g1\tg1\ng2\tg2\n")

    def test_fishash_score_is_larger_for_stronger_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            native = root / "native.csv"
            output = root / "assignments.csv"
            pd.DataFrame({
                "cell": ["c1", "c1"], "gRNA": ["g1", "g2"],
                "UMI_counts": [4, 3], "log_pval": [-12.0, -5.0],
            }).to_csv(native, index=False)
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "standardize_assignment.py"),
                "--input", str(native), "--output", str(output),
                "--method", "fishash",
            ], check=True)
            result = pd.read_csv(output)
            self.assertEqual(result["guide_id"].tolist(), ["g1", "g2"])
            self.assertEqual(result["score"].tolist(), [12.0, 5.0])
            self.assertEqual(result["native_score"].tolist(), [-12.0, -5.0])

    def test_perturbation_obs_matches_sanitized_guide_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            library = root / "library.csv"
            assignment = root / "assignments.csv"
            output = root / "perturbation_obs.csv"
            raw_guide = "GENE_+_123.45-P1P2"
            pd.DataFrame({
                "construct_id": ["construct_1"],
                "guide_id": [raw_guide],
                "target_label": ["GENE"],
            }).to_csv(library, index=False)
            pd.DataFrame({
                "cell_barcode": ["cell_1"],
                "guide_id": ["GENE___123_45-P1P2"],
                "rank": [1], "score": [0.99],
            }).to_csv(assignment, index=False)
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "make_perturbation_obs.py"),
                "--input", str(assignment), "--output", str(output),
                "--guide-design", "dual", "--guide-csv", str(library),
                "--method", "pgmm_em",
            ], check=True)
            result = pd.read_csv(output, keep_default_na=False)
            self.assertEqual(result.loc[0, "perturbation"], "construct_1")
            self.assertEqual(result.loc[0, "n_guides_assigned"], 1)

    def test_integration_excludes_expression_cells_without_assignment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expression = root / "expression.h5ad"
            assignment = root / "assignments.csv"
            output = root / "integrated.h5ad"
            cells = ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1", "GGGGGGGGGGGGGGGG-1"]
            counts = sparse.csr_matrix(np.asarray([[1, 0], [0, 2], [3, 1]], dtype=np.int32))
            ad.AnnData(
                X=counts,
                obs=pd.DataFrame(index=cells),
                var=pd.DataFrame(index=["A", "B"]),
            ).write_h5ad(expression)
            pd.DataFrame({
                "cell_barcode": cells[:2], "guide_id": ["g1", "g2"],
                "umi_count": [3, 4], "rank": [1, 1], "score": [0.9, 0.8],
                "score_type": ["prob_gaussian", "prob_gaussian"],
                "method": ["pgmm_em", "pgmm_em"],
            }).to_csv(assignment, index=False)
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "integrate_multimodal.py"),
                "--gex", f"lane_1={expression}", "--assign", str(assignment),
                "--guide-design", "single", "--input-kind", "auto",
                "--min-free-disk-gb", "0", "--max-process-memory-gb", "192",
                "--out", str(output),
            ], check=True)
            result = ad.read_h5ad(output)
            self.assertEqual(result.obs_names.tolist(), cells[:2])
            self.assertFalse(result.obs["guide_assignment_missing"].any())
            self.assertEqual(result.uns["manifest"]["excluded_unassigned_gex_cells"], 1)

    def test_integration_selects_dataset_level_cell_matching(self):
        def run_case(root, obs, assignment_cells):
            expression = root / "expression.h5ad"
            assignment = root / "assignments.csv"
            barcodes = root / "merged_barcodes.tsv"
            output = root / "integrated.h5ad"
            counts = sparse.csr_matrix(
                np.arange(1, len(obs) + 1, dtype=np.int32).reshape(-1, 1)
            )
            ad.AnnData(
                X=counts, obs=obs, var=pd.DataFrame(index=["A"])
            ).write_h5ad(expression)
            pd.DataFrame({
                "cell_barcode": assignment_cells,
                "guide_id": [f"g{i}" for i in range(len(assignment_cells))],
                "umi_count": 3, "rank": 1, "score": 0.9,
                "score_type": "prob_gaussian", "method": "pgmm_em",
            }).to_csv(assignment, index=False)
            barcodes.write_text("\n".join(assignment_cells) + "\n", encoding="utf-8")
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "integrate_multimodal.py"),
                "--gex", f"dataset={expression}", "--assign", str(assignment),
                "--barcodes", str(barcodes),
                "--guide-design", "single", "--input-kind", "auto",
                "--cell-id-matching", "auto", "--min-free-disk-gb", "0",
                "--max-process-memory-gb", "192", "--out", str(output),
            ], check=True)
            return ad.read_h5ad(output)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            papalexi = root / "papalexi"
            papalexi.mkdir()
            barcode_a = "AAAAAAAAAAAAAAAA"
            barcode_c = "CCCCCCCCCCCCCCCC"
            result = run_case(
                papalexi,
                pd.DataFrame(index=[
                    f"l1_{barcode_a}", f"l2_{barcode_a}", f"l2_{barcode_c}"
                ]),
                [f"{barcode_a}-L01", f"{barcode_c}-L01"],
            )
            self.assertEqual(
                result.uns["manifest"]["cell_id_matching_selected"], "pure_16mer"
            )
            self.assertEqual(result.uns["manifest"]["duplicate_expression_rows_removed"], 1)
            self.assertEqual(result.uns["manifest"]["merged_barcode_overlap"], 2)
            self.assertEqual(
                result.obs["source_cell_id"].tolist(),
                [f"l1_{barcode_a}", f"l2_{barcode_c}"],
            )

            replogle = root / "replogle"
            replogle.mkdir()
            obs = pd.DataFrame(
                {"barcode_16mer": [barcode_a, barcode_a], "lane": [1, 2]},
                index=["source_1", "source_2"],
            )
            result = run_case(
                replogle, obs, [f"{barcode_a}-L01", f"{barcode_a}-L02"]
            )
            self.assertEqual(
                result.uns["manifest"]["cell_id_matching_selected"], "barcode_lane"
            )
            self.assertEqual(result.n_obs, 2)

    def test_status_stage_skips_cleanly_without_eligible_cells(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.h5ad"
            output = root / "status"
            counts = sparse.csr_matrix(np.asarray([[1, 2]], dtype=np.int32))
            obs = pd.DataFrame({
                "guide_id": ["g1;g2"], "target_label": ["A+B"],
                "perturbation_group": ["A+B"], "is_ntc": [False],
                "batch_id": ["b1"], "assignment_structure": ["mixed_construct"],
            }, index=pd.Index(["c1"], name="cell_id"))
            adata = ad.AnnData(X=counts.astype(np.float32), obs=obs,
                               var=pd.DataFrame(index=["A", "B"]))
            adata.layers["counts"] = counts
            adata.write_h5ad(source)
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "run_perturbation_status.py"),
                "--method", "ps", "--input", str(source),
                "--output-dir", str(output),
            ], check=True)
            manifest = pd.read_json(output / "run_manifest.json", typ="series")
            self.assertTrue(bool(manifest["skipped"]))
            self.assertEqual(manifest["skip_reason"], "no_status_eligible_cells")
            self.assertFalse((output / "work").exists())


if __name__ == "__main__":
    unittest.main()
