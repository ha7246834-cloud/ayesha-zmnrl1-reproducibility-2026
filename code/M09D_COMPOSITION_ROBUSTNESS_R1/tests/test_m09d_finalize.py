from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_finalize as finalize  # noqa: E402


class FinalizationPrimitivesTests(unittest.TestCase):
    def test_support_lookup_is_root_invariant(self):
        universe = {"A", "B", "C", "D", "E"}
        supports = {("A", "B"): (90.0, 99.0)}
        self.assertEqual(finalize.support_for({"A", "B"}, universe, supports), "90/99")
        self.assertEqual(
            finalize.support_for({"C", "D", "E"}, universe, supports), "90/99"
        )
        self.assertEqual(finalize.support_for({"A", "C"}, universe, supports), "NOT_PRESENT")

    def test_no_rbh_rollup_uses_least_conclusive_outcome(self):
        self.assertEqual(
            finalize.rollup_relationships(["POSSIBLE_CO_ORTHOLOGY", "UNRESOLVED"]),
            "UNRESOLVED",
        )
        self.assertEqual(
            finalize.rollup_relationships(["POSSIBLE_CO_ORTHOLOGY", "CONTRADICTED"]),
            "CONTRADICTED",
        )
        self.assertEqual(
            finalize.rollup_relationships(["POSSIBLE_CO_ORTHOLOGY"]),
            "POSSIBLE_CO_ORTHOLOGY",
        )

    def test_complete_manifest_covers_all_files_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "nested").mkdir()
            (root / "a.txt").write_text("a\n", encoding="utf-8")
            (root / "nested/b.txt").write_text("b\n", encoding="utf-8")
            self.assertEqual(finalize._write_complete_manifest(root), 2)
            manifest = (root / "SHA256SUMS.txt").read_text(encoding="utf-8")
            self.assertIn("a.txt", manifest)
            self.assertIn("nested/b.txt", manifest)
            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                finalize._write_complete_manifest(root)


if __name__ == "__main__":
    unittest.main()
