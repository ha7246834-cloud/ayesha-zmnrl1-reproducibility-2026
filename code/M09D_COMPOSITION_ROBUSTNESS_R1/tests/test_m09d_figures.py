from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_figures as figures  # noqa: E402


class FigureTests(unittest.TestCase):
    def test_render_and_verify_three_traceable_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tree = root / "tree.nwk"
            tree.write_text(
                "((ZmNRG2_001:0.1,Sb_GeneID_1:0.1)99/100:0.2,"
                "(ZmNRG2_002:0.1,Sb_GeneID_2:0.1)85/96:0.2);\n",
                encoding="utf-8",
            )
            rows = [
                {"tree_label": "ZmNRG2_001", "species": "Zea mays", "classification": "ROBUST_COMPUTATIONAL_CANDIDATE"},
                {"tree_label": "ZmNRG2_002", "species": "Zea mays", "classification": "UNRESOLVED"},
                {"tree_label": "Sb_GeneID_1", "species": "Sorghum bicolor", "classification": "CONTEXT_ONLY"},
                {"tree_label": "Sb_GeneID_2", "species": "Sorghum bicolor", "classification": "CONTEXT_ONLY"},
            ]
            artifacts = figures.render_phylogram(
                tree, rows, root / "figure", "Synthetic unrooted robustness tree"
            )
            got = figures.verify_figure_artifacts(tree, rows, artifacts)
            self.assertEqual(got["leaf_count"], 4)
            self.assertEqual({path.suffix for path in artifacts}, {".svg", ".pdf", ".png"})
            self.assertTrue(got["all_labels_present_in_svg"])

    def test_duplicate_source_label_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tree = root / "tree.nwk"
            tree.write_text("(A:0.1,B:0.1,C:0.1,D:0.1);\n", encoding="utf-8")
            rows = [
                {"tree_label": "A", "species": "Zea mays"},
                {"tree_label": "A", "species": "Zea mays"},
                {"tree_label": "C", "species": "Sorghum bicolor"},
                {"tree_label": "D", "species": "Sorghum bicolor"},
            ]
            with self.assertRaisesRegex(ValueError, "labels"):
                figures.render_phylogram(tree, rows, root / "bad", "Bad")

    def test_render_robustness_matrix_accounts_for_23_unique_maize_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for index in range(1, 24):
                rows.append(
                    {
                        "maize_family_alias": f"ZmNRG2_{index:03d}",
                        "sorghum_candidate_ids": f"Sb_GeneID_{index}",
                        "reciprocal_best_hit_status": "YES" if index <= 16 else "NO",
                        "m09c_final_classification": "SUPPORTED_COMPUTATIONAL_ORTHOLOGY_CANDIDATE",
                        "r1_state": "STRONG",
                        "r2_state": "WEAK_OR_ABSENT",
                        "r3_state": "STRONG",
                        "r4_state": "STRONG",
                        "final_computational_classification": (
                            "ROBUST_COMPUTATIONAL_CANDIDATE"
                            if index <= 11
                            else "POSSIBLE_CO_ORTHOLOGY"
                        ),
                    }
                )
            artifacts = figures.render_robustness_matrix(
                rows, root / "matrix", "Synthetic 23-gene robustness matrix"
            )
            got = figures.verify_matrix_artifacts(rows, artifacts)
            self.assertEqual(got["maize_gene_count"], 23)
            self.assertEqual({path.suffix for path in artifacts}, {".svg", ".pdf", ".png"})

    def test_robustness_matrix_rejects_missing_maize_gene(self):
        rows = [
            {
                "maize_family_alias": f"ZmNRG2_{index:03d}",
                "sorghum_candidate_ids": "Sb",
                "reciprocal_best_hit_status": "YES",
                "m09c_final_classification": "X",
                "r1_state": "STRONG",
                "r2_state": "STRONG",
                "r3_state": "STRONG",
                "r4_state": "STRONG",
                "final_computational_classification": "ROBUST_COMPUTATIONAL_CANDIDATE",
            }
            for index in range(1, 23)
        ]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "23 unique"):
                figures.render_robustness_matrix(rows, Path(tmp) / "bad", "Bad")


if __name__ == "__main__":
    unittest.main()
