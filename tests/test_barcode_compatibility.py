from __future__ import annotations

import csv
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
CHECKER = ROOT / "scripts" / "check_barcode_compatibility.py"


class BarcodeCompatibilityTests(unittest.TestCase):
    def run_case(self, expression_ids: list[str], assignment_ids: list[str]) -> dict:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expression = root / "expression.h5ad"
            assignment = root / "assignment.csv"
            ad.AnnData(
                X=sparse.csr_matrix(np.zeros((len(expression_ids), 1), dtype=np.int8)),
                obs=pd.DataFrame(index=expression_ids),
                var=pd.DataFrame(index=["gene"]),
            ).write_h5ad(expression)
            with assignment.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=["cell"])
                writer.writeheader()
                writer.writerows({"cell": cell} for cell in assignment_ids)
            result = subprocess.run(
                [sys.executable, str(CHECKER), "--expression", str(expression),
                 "--assignment", str(assignment)],
                text=True, capture_output=True, check=False,
            )
            return json.loads(result.stdout)

    def test_exact_ids(self):
        report = self.run_case(["AAAAAAAAAAAAAAAA-1"], ["AAAAAAAAAAAAAAAA-1"])
        self.assertEqual(report["verdict"], "PASS")
        self.assertEqual(report["exact_id_matches"], 1)

    def test_prefix_and_suffix_lane_formats(self):
        barcode = "CCCCCCCCCCCCCCCC"
        report = self.run_case(
            [f"l01_{barcode}", f"lane02_{barcode}"],
            [f"{barcode}-L01", f"{barcode}-2"],
        )
        self.assertEqual(report["verdict"], "PASS")
        self.assertEqual(report["explicit_or_structured_group_matches"], 2)

    def test_unique_barcode_fallback(self):
        barcode = "GGGGGGGGGGGGGGGG"
        report = self.run_case([f"sample_{barcode}"], [barcode])
        self.assertEqual(report["verdict"], "PASS")
        self.assertEqual(report["unique_16mer_matches"], 1)

    def test_repeated_barcode_without_group_is_rejected(self):
        barcode = "TTTTTTTTTTTTTTTT"
        report = self.run_case(
            [f"sampleA_{barcode}", f"sampleB_{barcode}"], [barcode]
        )
        self.assertEqual(report["verdict"], "NEEDS_GROUP_INPUT", report)
        self.assertEqual(report["unresolved_with_expression_candidates"], 1)


if __name__ == "__main__":
    unittest.main()
