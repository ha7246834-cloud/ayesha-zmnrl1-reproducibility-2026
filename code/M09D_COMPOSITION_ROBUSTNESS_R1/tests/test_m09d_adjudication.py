from __future__ import annotations

import sys
import unittest
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOL_DIR))
import m09d_core as core  # noqa: E402


class SplitTests(unittest.TestCase):
    def test_supported_splits_are_unrooted_and_require_two_support_metrics(self):
        newick = "((A:0.1,B:0.1)90/99:0.2,(C:0.1,(D:0.1,E:0.1)85/96:0.2)70:0.1);"
        got = core.parse_supported_splits(newick)
        self.assertEqual(got[("A", "B")], (90.0, 99.0))
        self.assertEqual(got[("D", "E")], (85.0, 96.0))
        self.assertNotIn(("C", "D", "E"), got)

    def test_reproducible_exact_alternative_requires_same_pair_twice(self):
        target = frozenset({"Zm", "Sb1"})
        trees = {
            "R1": {frozenset({"Zm", "Sb2"}), frozenset({"X", "Y"})},
            "R2": {frozenset({"Zm", "Sb2"})},
            "R3": {frozenset({"Sb1", "Sb3"})},
            "R4": set(),
        }
        got = core.reproducible_exact_alternatives(target, trees)
        self.assertEqual(got, {("Sb2", "Zm"): ["R1", "R2"]})


class LockedClassificationTests(unittest.TestCase):
    def _call(self, **changes):
        values = {
            "rbh": True,
            "original_exact": {"M09C_PRIMARY": True, "M09C_FULL": True},
            "new_exact": {"R1": True, "R2": True, "R3": True, "R4": False},
            "original_component": {"M09C_PRIMARY": True, "M09C_FULL": True},
            "new_component": {"R1": True, "R2": True, "R3": True, "R4": True},
            "reproducible_alternative_count": 0,
            "direct_caveat": False,
        }
        values.update(changes)
        return core.adjudicate_m09d_relationship(**values)

    def test_robust_candidate_requires_r1_and_an_occupancy_tree(self):
        self.assertEqual(self._call(), "ROBUST_COMPUTATIONAL_CANDIDATE")
        self.assertEqual(
            self._call(new_exact={"R1": False, "R2": True, "R3": True, "R4": True}),
            "SUPPORTED_WITH_MODEL_CAVEAT",
        )

    def test_direct_caveat_prevents_robust_class(self):
        self.assertEqual(
            self._call(direct_caveat=True), "SUPPORTED_WITH_MODEL_CAVEAT"
        )

    def test_reproducible_alternative_has_first_priority(self):
        self.assertEqual(
            self._call(reproducible_alternative_count=1),
            "CONTRADICTED",
        )

    def test_no_rbh_can_only_reach_component_level_class(self):
        self.assertEqual(
            self._call(rbh=False),
            "POSSIBLE_CO_ORTHOLOGY",
        )
        self.assertEqual(
            self._call(
                rbh=False,
                new_component={"R1": True, "R2": False, "R3": False, "R4": False},
            ),
            "UNRESOLVED",
        )

    def test_component_support_is_not_promoted_to_exact_candidate(self):
        self.assertEqual(
            self._call(
                new_exact={"R1": False, "R2": False, "R3": False, "R4": False}
            ),
            "POSSIBLE_CO_ORTHOLOGY",
        )

    def test_missing_new_analysis_prevents_robust_class(self):
        self.assertEqual(
            self._call(new_exact={"R1": True, "R2": True, "R3": True, "R4": None}),
            "SUPPORTED_WITH_MODEL_CAVEAT",
        )


if __name__ == "__main__":
    unittest.main()

