"""Synthetic unit fixtures only; never biological evidence.

Mutations caught: wrong isoform acceptance; internal-stop removal; codon frame
shift; gap-only columns counted as data; zero/sentinel Ks emitted as a ratio;
and biological work admitted before the M09D completion gate.
"""
import importlib
import sys
from pathlib import Path

import unittest
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


class CoreTests(unittest.TestCase):
  @property
  def core(self):
    self.assertIsNotNone(importlib.util.find_spec("kaks_core"), "Ka/Ks core is not implemented")
    return importlib.import_module("kaks_core")


  def test_translation_checks_exact_isoform_and_terminal_stop(self):
    core = self.core
    assert core.validated_cds("ATGGCTTAA", "MA") == "ATGGCT"
    assert core.validated_cds("ATGGCT", "MA") == "ATGGCT"
    with self.assertRaises(ValueError):
        core.validated_cds("ATGGCTTAA", "MV")
    with self.assertRaises(ValueError):
        core.validated_cds("ATGTAAGCTTAA", "MA")
    with self.assertRaises(ValueError):
        core.validated_cds("ATGGCTTA", "MA")
    with self.assertRaises(ValueError):
        core.validated_cds("ATGNNNTAA", "MX")


  def test_backtranslation_is_frame_preserving_and_drops_only_gap_pairs(self):
    core = self.core
    alignment = {"seqA": "MA-G", "seqB": "M-VG"}
    cds = {"seqA": "ATGGCTGGT", "seqB": "ATGGTTGGT"}
    result, indices = core.independent_paired_codons(alignment, cds)
    assert result == {"seqA": "ATGGGT", "seqB": "ATGGGT"}
    assert indices == [1, 4]
    with self.assertRaises(ValueError):
        core.independent_paired_codons(alignment, {"seqA": "ATGGTTGGT", "seqB": cds["seqB"]})


  def test_alignment_metrics_do_not_count_gaps_as_identity(self):
    core = self.core
    metrics = core.alignment_metrics({"seqA": "MA--GG", "seqB": "MAVT-G"})
    assert metrics["paired_codons"] == 3
    assert metrics["protein_identity_fraction"] == 1.0
    assert metrics["retained_fraction_shorter"] == 0.75
    assert metrics["retained_fraction_longer"] == 0.6
    assert metrics["maximum_gap_run"] == 2


  def test_rate_screen_distinguishes_zero_missing_and_saturation_warning(self):
    core = self.core
    zero = core.rate_screen(0.0, 0.0, 0.0)
    assert zero["omega"] is None and "ZERO_KS" in zero["flags"]
    valid = core.rate_screen(0.1, 0.5, 0.02)
    assert valid["omega"] == 0.2
    assert valid["scientific_status"] == "REVIEW"
    assert "SATURATION_NOT_EXCLUDED" in valid["flags"]
    assert "HIGH_KS" in core.rate_screen(0.1, 2.5, 0.2)["flags"]
    assert core.rate_screen(0.1, 99.0, 0.2)["omega"] is None
    assert core.rate_screen(float("nan"), 0.5, 0.02)["omega"] is None
    assert core.rate_screen(0.1, 99.0, 0.2)["dS"] is None
    assert core.rate_screen(-0.1, 0.5, 0.02)["dN"] is None
    assert "KS_SE_INTERVAL_INCLUDES_ZERO" in core.rate_screen(0.1, 0.2, 0.2)["flags"]


  def test_completion_gate_requires_each_real_analysis(self):
    core = self.core
    with tempfile.TemporaryDirectory() as temp:
      self._completion_gate(core, Path(temp))

  def _completion_gate(self, core, tmp_path):
    import json
    (tmp_path / "CHECKPOINTS").mkdir()
    with self.assertRaises(ValueError):
        core.require_m09d_completion(tmp_path)
    for name in ("R1", "R2", "R3", "R4"):
        (tmp_path / "CHECKPOINTS" / f"{name}_COMPLETE.json").write_text(json.dumps({"analysis": name}))
    with self.assertRaises(ValueError):
        core.require_m09d_completion(tmp_path)
    # Checkpoint filenames alone never suffice; adjudication is mandatory.


  def test_empty_alignment_and_unequal_lengths_fail_closed(self):
    core = self.core
    with self.assertRaises(ValueError):
        core.alignment_metrics({"seqA": "---", "seqB": "---"})
    with self.assertRaises(ValueError):
        core.alignment_metrics({"seqA": "MA", "seqB": "M"})

  def test_completed_review_release_is_allowed_but_tampering_is_rejected(self):
    import json
    import hashlib
    core = self.core
    with tempfile.TemporaryDirectory() as directory:
      release = Path(directory)
      (release / 'CHECKPOINTS').mkdir()
      for analysis in ('R1', 'R2', 'R3', 'R4'):
        (release / 'CHECKPOINTS' / f'{analysis}_COMPLETE.json').write_text(json.dumps({'analysis': analysis}))
      (release / 'M09D_QC.json').write_text(json.dumps({'status': 'REVIEW', 'analyses': {x: {} for x in ('R1', 'R2', 'R3', 'R4')}}))
      matrix = release / 'M09D_ALL_23_FINAL_REVIEW.tsv'
      matrix.write_text('maize_family_alias\n' + ''.join(f'ZmNRG2_{i:03d}\n' for i in range(1,24)))
      lines = [f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(release)}' for p in release.rglob('*') if p.is_file()]
      (release / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n')
      self.assertEqual(core.require_m09d_completion(release)['status'], 'REVIEW')
      matrix.write_text('altered')
      with self.assertRaises(ValueError):
        core.require_m09d_completion(release)

if __name__ == "__main__":
    unittest.main()
