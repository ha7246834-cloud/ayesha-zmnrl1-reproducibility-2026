from __future__ import annotations

import sys
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_core as core  # noqa: E402


class TransformTests(unittest.TestCase):
    def test_occupancy_mask_uses_standard_residues_and_records_source_columns(self):
        records = {"A": "AC-D", "B": "A--D", "C": "ACXD", "D": "ACDD"}
        got = core.mask_alignment_by_occupancy(records, minimum_present=2)
        self.assertEqual(got["source_columns"], 4)
        self.assertEqual(got["retained_columns"], 3)
        self.assertEqual(got["retained_source_columns_1based"], [1, 2, 4])
        self.assertEqual(
            got["records"],
            {"A": "ACD", "B": "A-D", "C": "ACD", "D": "ACD"},
        )

    def test_partition_masking_regenerates_dynamic_endpoints(self):
        records = {
            "A": "AC-DFGH",
            "B": "A--D-GH",
            "C": "ACXDF-H",
            "D": "ACDDFGH",
        }
        got = core.mask_partitioned_alignment(
            records,
            [("PF04783", 1, 3), ("PF04782", 4, 7)],
            minimum_present=2,
        )
        self.assertEqual(got["retained_source_columns_1based"], [1, 2, 4, 5, 6, 7])
        self.assertEqual(
            got["partitions"],
            [
                {"name": "PF04783", "start": 1, "end": 2, "retained": 2},
                {"name": "PF04782", "start": 3, "end": 6, "retained": 4},
            ],
        )
        self.assertEqual({len(seq) for seq in got["records"].values()}, {6})

    def test_partition_masking_fails_if_a_partition_becomes_empty(self):
        records = {"A": "A--", "B": "A--", "C": "A--"}
        with self.assertRaisesRegex(ValueError, "empty partition"):
            core.mask_partitioned_alignment(
                records,
                [("P1", 1, 1), ("P2", 2, 3)],
                minimum_present=2,
            )

    def test_dayhoff6_recode_is_literal_and_preserves_missingness(self):
        source = "AGPSTDENQHKRILMVFWYC-?X"
        expected = "".join(
            ["0"] * 5
            + ["1"] * 4
            + ["2"] * 3
            + ["3"] * 4
            + ["4"] * 3
            + ["5"]
            + ["?"] * 3
        )
        got = core.dayhoff6_recode({"T": source})
        self.assertEqual(got, {"T": expected})
        self.assertEqual(len(got["T"]), len(source))

    def test_dayhoff6_recode_fails_on_unapproved_symbol(self):
        with self.assertRaisesRegex(ValueError, "unapproved"):
            core.dayhoff6_recode({"T": "ACD#"})


if __name__ == "__main__":
    unittest.main()

