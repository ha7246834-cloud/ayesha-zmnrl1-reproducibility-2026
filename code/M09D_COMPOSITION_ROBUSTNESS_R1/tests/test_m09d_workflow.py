from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_workflow as workflow  # noqa: E402


class CommandLockTests(unittest.TestCase):
    def test_r1_site_heterogeneous_command_is_pre_registered(self):
        command = workflow.build_iqtree_command(
            "/x/iqtree3", "R1", Path("dual.faa"), Path("r1/prefix")
        )
        self.assertEqual(command[0], "/x/iqtree3")
        self.assertEqual(command[command.index("--seqtype") + 1], "AA")
        self.assertEqual(command[command.index("-m") + 1], "LG+C20+F+G4")
        self.assertEqual(command[command.index("--seed") + 1], "20261108")
        self.assertNotIn("-p", command)
        self._assert_shared_lock(command)

    def test_r2_partitioned_masked_command_is_pre_registered(self):
        command = workflow.build_iqtree_command(
            "/x/iqtree3",
            "R2",
            Path("dual.masked.faa"),
            Path("r2/prefix"),
            partitions=Path("masked.nex"),
        )
        self.assertEqual(command[command.index("-m") + 1], "MFP+MERGE")
        self.assertEqual(command[command.index("--seed") + 1], "20261109")
        self.assertEqual(command[command.index("-p") + 1], "masked.nex")
        self._assert_shared_lock(command)

    def test_r3_full_length_masked_command_is_pre_registered(self):
        command = workflow.build_iqtree_command(
            "/x/iqtree3", "R3", Path("full.masked.faa"), Path("r3/prefix")
        )
        self.assertEqual(command[command.index("-m") + 1], "MFP")
        self.assertEqual(command[command.index("--seed") + 1], "20261110")
        self._assert_shared_lock(command)

    def test_r4_dayhoff6_command_is_pre_registered(self):
        command = workflow.build_iqtree_command(
            "/x/iqtree3", "R4", Path("dual.dayhoff6.faa"), Path("r4/prefix")
        )
        self.assertEqual(command[command.index("--seqtype") + 1], "MORPH{6}")
        self.assertEqual(command[command.index("-m") + 1], "GTR+FO+G4")
        self.assertEqual(command[command.index("--seed") + 1], "20261111")
        self._assert_shared_lock(command)

    def test_unknown_analysis_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "unknown analysis"):
            workflow.build_iqtree_command(
                "/x/iqtree3", "POST_HOC", Path("x"), Path("y")
            )

    def _assert_shared_lock(self, command):
        self.assertEqual(command[command.index("-T") + 1], "2")
        self.assertEqual(command[command.index("--runs") + 1], "5")
        self.assertEqual(command[command.index("-B") + 1], "1000")
        self.assertEqual(command[command.index("--alrt") + 1], "1000")
        self.assertIn("-bnni", command)
        self.assertIn("--boot-trees", command)


class IntegrityAndResumeTests(unittest.TestCase):
    def test_prepare_checkpoint_excludes_mutable_runtime_records(self):
        self.assertFalse(workflow.prepare_checkpoint_includes(Path("RUN_STATE.json")))
        self.assertFalse(workflow.prepare_checkpoint_includes(Path("PERFORMANCE_LOG.tsv")))
        self.assertFalse(
            workflow.prepare_checkpoint_includes(Path("CHECKPOINTS/R1_COMPLETE.json"))
        )
        self.assertTrue(
            workflow.prepare_checkpoint_includes(Path("M09D_COMPOSITION_DIAGNOSTICS.tsv"))
        )

    def test_iqtree_nan_probability_is_missing_not_zero(self):
        self.assertIsNone(workflow.parse_optional_probability("nan"))
        self.assertIsNone(workflow.parse_optional_probability("NA"))
        self.assertIsNone(workflow.parse_optional_probability(""))
        self.assertEqual(workflow.parse_optional_probability("0.0488069"), 0.0488069)

    def test_manifest_verification_detects_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.txt"
            payload.write_text("locked\n", encoding="utf-8")
            digest = workflow.sha256_file(payload)
            manifest = root / "SHA256SUMS"
            manifest.write_text(f"{digest}  payload.txt\n", encoding="utf-8")
            self.assertEqual(workflow.verify_manifest(root, manifest), 1)
            payload.write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum"):
                workflow.verify_manifest(root, manifest)

    def test_manifest_rejects_path_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "SHA256SUMS"
            manifest.write_text(f"{'0' * 64}  ../outside\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsafe manifest"):
                workflow.verify_manifest(root, manifest)

    def test_safe_staging_must_be_exclusive_hidden_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / "release-parent"
            parent.mkdir()
            good = parent / ".STAGING_20261108T000000Z_abcdef12"
            self.assertEqual(workflow.safe_staging_path(parent, good), good.resolve())
            with self.assertRaisesRegex(ValueError, "staging"):
                workflow.safe_staging_path(parent, parent / "RELEASE_bad")
            with self.assertRaisesRegex(ValueError, "staging"):
                workflow.safe_staging_path(parent, parent.parent / ".STAGING_escape")

    def test_provenance_binds_command_input_executable_and_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alignment = root / "a.faa"
            executable = root / "iqtree3"
            code = root / "workflow.py"
            alignment.write_text(">A\nACD\n", encoding="utf-8")
            executable.write_text("binary", encoding="utf-8")
            code.write_text("code", encoding="utf-8")
            got = workflow.build_provenance(
                [str(executable), "-s", str(alignment)],
                [alignment],
                [executable],
                [code],
            )
            self.assertEqual(got["schema"], "M09D_PROVENANCE_V1")
            self.assertEqual(got["inputs"][str(alignment)], workflow.sha256_file(alignment))
            changed = dict(got)
            changed["command"] = ["different"]
            self.assertFalse(workflow.provenance_matches(got, changed))

    def test_iqtree_metrics_parse_fixed_and_model_finder_reports(self):
        report = """
Best-fit model according to BIC: LG+F+R4
Log-likelihood of the tree: -1234.5678 (s.e. 0.0000)
Number of free parameters (#branches + #model parameters): 88
Akaike information criterion (AIC) score: 2645.1356
Corrected Akaike information criterion (AICc) score: 2700.1000
Bayesian information criterion (BIC) score: 2999.9000
"""
        got = workflow.parse_iqtree_metrics(report, declared_model="MFP")
        self.assertEqual(got["selected_model"], "LG+F+R4")
        self.assertEqual(got["log_likelihood"], -1234.5678)
        self.assertEqual(got["free_parameters"], 88)
        self.assertEqual(got["aic"], 2645.1356)
        self.assertEqual(got["aicc"], 2700.1)
        self.assertEqual(got["bic"], 2999.9)

    def test_iqtree_completion_requires_1000_bootstraps_and_five_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "analysis"
            for suffix in (".treefile", ".contree"):
                prefix.with_suffix(suffix).write_text("(A:0.1,B:0.1,C:0.1,D:0.1);\n")
            prefix.with_suffix(".iqtree").write_text(
                "Log-likelihood of the tree: -1\n", encoding="utf-8"
            )
            prefix.with_suffix(".log").write_text(
                "Analysis results written to:\n", encoding="utf-8"
            )
            prefix.with_suffix(".ufboot").write_text("x\n" * 1000, encoding="utf-8")
            prefix.with_suffix(".runtrees").write_text("x\n" * 5, encoding="utf-8")
            self.assertTrue(workflow.iqtree_artifacts_complete(prefix, {"A", "B", "C", "D"}))
            prefix.with_suffix(".ufboot").write_text("x\n" * 999, encoding="utf-8")
            self.assertFalse(workflow.iqtree_artifacts_complete(prefix, {"A", "B", "C", "D"}))

    def test_alignment_column_count_accepts_morphology_and_rejects_unequal_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "morph.faa"
            path.write_text(">A\n0123?\n>B\n5432?\n", encoding="utf-8")
            self.assertEqual(workflow.alignment_column_count(path), 5)
            path.write_text(">A\n0123?\n>B\n543?\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unequal"):
                workflow.alignment_column_count(path)


if __name__ == "__main__":
    unittest.main()
