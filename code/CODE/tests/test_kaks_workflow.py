"""Fail-closed orchestration tests; synthetic fixtures are not research data."""
import importlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


class WorkflowTests(unittest.TestCase):
    @property
    def workflow(self):
        self.assertIsNotNone(importlib.util.find_spec("kaks_workflow"), "Workflow is not implemented")
        return importlib.import_module("kaks_workflow")

    def test_incomplete_m09d_cannot_create_downstream_outputs(self):
        w = self.workflow
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(ValueError):
                w.run(root, root / "missing_m09d", root / "new_candidate", "align")
            self.assertFalse((root / "new_candidate").exists())

    def test_pairs_require_gene_alias_agreement_and_unique_unordered_pairs(self):
        w = self.workflow
        mappings = [{"tree_label": "a", "gene_id": "geneA"}, {"tree_label": "b", "gene_id": "geneB"}]
        pair = {"alias_a": "a", "gene_a": "geneA", "alias_b": "b", "gene_b": "geneB"}
        self.assertEqual(w.validate_pairs([pair], mappings), [pair])
        with self.assertRaises(ValueError):
            w.validate_pairs([dict(pair, gene_b="geneX")], mappings)
        with self.assertRaises(ValueError):
            w.validate_pairs([pair, dict(alias_a="b", gene_a="geneB", alias_b="a", gene_b="geneA")], mappings)

    def test_pair_gate_does_not_equate_short_or_divergent_alignment_with_pass(self):
        w = self.workflow
        self.assertEqual(w.alignment_gate({"paired_codons": 600, "retained_fraction_shorter": 0.95, "protein_identity_fraction": 0.8}), "PASS_TECHNICAL_SCREEN_PENDING_VISUAL_REVIEW")
        self.assertIn("REVIEW", w.alignment_gate({"paired_codons": 20, "retained_fraction_shorter": 0.1, "protein_identity_fraction": 0.2}))


if __name__ == "__main__":
    unittest.main()
