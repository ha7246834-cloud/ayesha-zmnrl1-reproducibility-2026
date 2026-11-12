from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_core as core  # noqa: E402


class SourceAndDiagnosticTests(unittest.TestCase):
    def test_parse_fasta_rejects_duplicate_ids(self):
        self.assertTrue(hasattr(core, "parse_fasta"), "parse_fasta is required")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "duplicate.faa"
            path.write_text(">A\nACD\n>A\nEFG\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                core.parse_fasta(path)

    def test_parse_fasta_preserves_aligned_missing_symbols(self):
        self.assertTrue(hasattr(core, "parse_fasta"), "parse_fasta is required")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aligned.faa"
            path.write_text(">A\nACD-X?\n>B\nA-DZX?\n", encoding="utf-8")
            got = core.parse_fasta(path, aligned=True)
        self.assertEqual(got, {"A": "ACD-X?", "B": "A-DZX?"})

    def test_validate_taxa_requires_exact_23_maize_and_17_sorghum(self):
        self.assertTrue(hasattr(core, "validate_taxa"), "validate_taxa is required")
        maize = {f"ZmNRG2_{i:03d}" for i in range(1, 24)}
        sorghum = {f"Sb_GeneID_{8000000 + i}" for i in range(17)}
        got = core.validate_taxa(maize | sorghum)
        self.assertEqual(got, {"maize": 23, "sorghum": 17, "total": 40})
        with self.assertRaisesRegex(ValueError, "23 maize"):
            core.validate_taxa((maize - {"ZmNRG2_023"}) | sorghum)

    def test_alignment_metrics_are_gap_aware_and_hand_derived(self):
        self.assertTrue(hasattr(core, "alignment_metrics"), "alignment_metrics is required")
        records = {"A": "ACD-", "B": "A-DX"}
        got = core.alignment_metrics(records)
        self.assertEqual(got["taxa"], 2)
        self.assertEqual(got["columns"], 4)
        self.assertEqual(got["standard_residues"], 5)
        self.assertEqual(got["gap_or_ambiguity"], 3)
        self.assertAlmostEqual(got["gap_or_ambiguity_fraction"], 3 / 8)
        self.assertEqual(got["per_taxon"]["A"]["nongap_standard"], 3)
        self.assertAlmostEqual(got["per_taxon"]["B"]["gap_or_ambiguity_fraction"], 0.5)

    def test_composition_diagnostics_zero_for_identical_composition(self):
        self.assertTrue(
            hasattr(core, "composition_diagnostics"),
            "composition_diagnostics is required",
        )
        rows = core.composition_diagnostics({"A": "ACDE", "B": "ACDE", "C": "ACDE"})
        self.assertEqual([row["taxon"] for row in rows], ["A", "B", "C"])
        for row in rows:
            self.assertEqual(row["observed_standard_residues"], 4)
            self.assertAlmostEqual(row["jensen_shannon_distance"], 0.0, places=12)

    def test_composition_diagnostics_ignores_gap_and_ambiguity(self):
        self.assertTrue(
            hasattr(core, "composition_diagnostics"),
            "composition_diagnostics is required",
        )
        rows = core.composition_diagnostics({"A": "AAAA--XX", "B": "CCCC----"})
        by_taxon = {row["taxon"]: row for row in rows}
        self.assertEqual(by_taxon["A"]["observed_standard_residues"], 4)
        self.assertEqual(by_taxon["B"]["observed_standard_residues"], 4)
        self.assertAlmostEqual(by_taxon["A"]["gap_or_ambiguity_fraction"], 0.5)
        self.assertGreater(by_taxon["A"]["jensen_shannon_distance"], 0.0)

    def test_benjamini_hochberg_matches_hand_calculation(self):
        self.assertTrue(
            hasattr(core, "benjamini_hochberg"),
            "benjamini_hochberg is required",
        )
        got = core.benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
        want = [0.02, 0.04, 0.04, 0.008]
        self.assertEqual(len(got), len(want))
        for observed, expected in zip(got, want):
            self.assertAlmostEqual(observed, expected, places=12)

    def test_benjamini_hochberg_preserves_missing_values(self):
        self.assertTrue(
            hasattr(core, "benjamini_hochberg"),
            "benjamini_hochberg is required",
        )
        self.assertEqual(core.benjamini_hochberg([0.01, None, 0.02]), [0.02, None, 0.02])

    def test_terminal_branch_diagnostics_flags_only_mad_outlier(self):
        self.assertTrue(
            hasattr(core, "terminal_branch_diagnostics"),
            "terminal_branch_diagnostics is required",
        )
        rows = core.terminal_branch_diagnostics(
            {"A": 0.10, "B": 0.11, "C": 0.09, "D": 1.00}
        )
        by_taxon = {row["taxon"]: row for row in rows}
        self.assertFalse(by_taxon["A"]["long_branch_flag"])
        self.assertFalse(by_taxon["B"]["long_branch_flag"])
        self.assertFalse(by_taxon["C"]["long_branch_flag"])
        self.assertTrue(by_taxon["D"]["long_branch_flag"])
        self.assertEqual(
            {row["threshold_method"] for row in rows},
            {"MEDIAN_PLUS_3_SCALED_MAD"},
        )

    def test_terminal_branch_diagnostics_uses_tukey_when_mad_zero(self):
        self.assertTrue(
            hasattr(core, "terminal_branch_diagnostics"),
            "terminal_branch_diagnostics is required",
        )
        rows = core.terminal_branch_diagnostics(
            {"A": 0.1, "B": 0.1, "C": 0.1, "D": 0.1}
        )
        self.assertEqual({row["threshold_method"] for row in rows}, {"TUKEY_UPPER_FENCE"})
        self.assertFalse(any(row["long_branch_flag"] for row in rows))


if __name__ == "__main__":
    unittest.main()

