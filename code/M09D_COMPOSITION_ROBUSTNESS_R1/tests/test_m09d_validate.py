from __future__ import annotations

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_validate as validate  # noqa: E402


class IndependentValidatorTests(unittest.TestCase):
    def test_manifest_must_cover_every_non_manifest_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.txt"
            payload.write_text("locked\n", encoding="utf-8")
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            manifest = root / "SHA256SUMS.txt"
            manifest.write_text(f"{digest}  payload.txt\n", encoding="utf-8")
            self.assertEqual(validate.verify_complete_manifest(root, manifest), 1)
            (root / "unlisted.txt").write_text("not locked\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "coverage"):
                validate.verify_complete_manifest(root, manifest)

    def test_manifest_rejects_path_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "SHA256SUMS.txt"
            manifest.write_text(f"{'0' * 64}  ../escape.txt\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsafe"):
                validate.verify_complete_manifest(root, manifest)


if __name__ == "__main__":
    unittest.main()
