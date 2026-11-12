"""Scientific core for the M09D composition-robustness workflow.

This module deliberately avoids third-party dependencies so an independent
validator can recompute basic diagnostics without importing producer code.
"""

from __future__ import annotations

import math
import re
import statistics
from io import StringIO
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from Bio import Phylo


STANDARD_AA = tuple("ACDEFGHIKLMNPQRSTVWY")
STANDARD_AA_SET = frozenset(STANDARD_AA)
ALIGNED_MISSING = frozenset("-?XBZJOU*")
UNALIGNED_MISSING = frozenset("X*")
MAIZE_ID = re.compile(r"^ZmNRG2_[0-9]{3}$")
SORGHUM_ID = re.compile(r"^Sb_GeneID_[0-9]+$")


def parse_fasta(path: str | Path, *, aligned: bool = False) -> dict[str, str]:
    """Parse FASTA using the first header token as a unique identifier."""

    fasta_path = Path(path)
    records: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []
    permitted = STANDARD_AA_SET | (ALIGNED_MISSING if aligned else UNALIGNED_MISSING)

    def commit() -> None:
        nonlocal current_id, chunks
        if current_id is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"empty FASTA sequence for {current_id!r} in {fasta_path}")
        invalid = sorted(set(sequence) - permitted)
        if invalid:
            raise ValueError(
                f"invalid FASTA symbols for {current_id!r}: {''.join(invalid)}"
            )
        records[current_id] = sequence

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header at line {line_number}")
                identifier = header.split()[0]
                if identifier in records or identifier == current_id:
                    raise ValueError(f"duplicate FASTA identifier: {identifier}")
                current_id = identifier
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"FASTA sequence occurs before a header at line {line_number}"
                    )
                chunks.append(line)
    commit()
    if not records:
        raise ValueError(f"no FASTA records in {fasta_path}")
    return records


def validate_taxa(labels: Iterable[str]) -> dict[str, int]:
    """Require the approved M09D taxon design: 23 maize plus 17 Sorghum."""

    label_list = list(labels)
    if len(set(label_list)) != len(label_list):
        raise ValueError("duplicate taxon identifiers")
    maize = [label for label in label_list if MAIZE_ID.fullmatch(label)]
    sorghum = [label for label in label_list if SORGHUM_ID.fullmatch(label)]
    unknown = sorted(set(label_list) - set(maize) - set(sorghum))
    if len(maize) != 23 or len(sorghum) != 17 or len(label_list) != 40 or unknown:
        raise ValueError(
            "expected exactly 23 maize and 17 Sorghum taxa; "
            f"observed maize={len(maize)}, sorghum={len(sorghum)}, "
            f"total={len(label_list)}, unknown={unknown}"
        )
    return {"maize": 23, "sorghum": 17, "total": 40}


def _validate_aligned_records(records: Mapping[str, str]) -> int:
    if not records:
        raise ValueError("alignment contains no records")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"alignment sequences have unequal lengths: {sorted(lengths)}")
    columns = next(iter(lengths))
    if columns == 0:
        raise ValueError("alignment has zero columns")
    return columns


def alignment_metrics(records: Mapping[str, str]) -> dict[str, object]:
    """Return hand-auditable alignment occupancy and missingness metrics."""

    columns = _validate_aligned_records(records)
    per_taxon: dict[str, dict[str, float | int]] = {}
    standard_total = 0
    for taxon in sorted(records):
        sequence = records[taxon].upper()
        standard = sum(symbol in STANDARD_AA_SET for symbol in sequence)
        missing = columns - standard
        standard_total += standard
        per_taxon[taxon] = {
            "nongap_standard": standard,
            "gap_or_ambiguity": missing,
            "gap_or_ambiguity_fraction": missing / columns,
        }
    cells = len(records) * columns
    missing_total = cells - standard_total
    return {
        "taxa": len(records),
        "columns": columns,
        "standard_residues": standard_total,
        "gap_or_ambiguity": missing_total,
        "gap_or_ambiguity_fraction": missing_total / cells,
        "per_taxon": per_taxon,
    }


def _frequency_vector(sequence: str) -> tuple[list[float], int]:
    counts = [sequence.upper().count(amino_acid) for amino_acid in STANDARD_AA]
    observed = sum(counts)
    if observed == 0:
        return [0.0] * len(STANDARD_AA), 0
    return [count / observed for count in counts], observed


def _kl_divergence(left: Sequence[float], right: Sequence[float]) -> float:
    total = 0.0
    for p_value, q_value in zip(left, right):
        if p_value > 0.0:
            if q_value <= 0.0:
                raise ValueError("KL reference probability is zero")
            total += p_value * math.log2(p_value / q_value)
    return total


def composition_diagnostics(records: Mapping[str, str]) -> list[dict[str, object]]:
    """Measure each sequence against pooled amino-acid composition."""

    columns = _validate_aligned_records(records)
    pooled_counts = [0] * len(STANDARD_AA)
    vectors: dict[str, tuple[list[float], int]] = {}
    for taxon, sequence in records.items():
        vector, observed = _frequency_vector(sequence)
        vectors[taxon] = (vector, observed)
        for index, amino_acid in enumerate(STANDARD_AA):
            pooled_counts[index] += sequence.upper().count(amino_acid)
    pooled_total = sum(pooled_counts)
    if pooled_total == 0:
        raise ValueError("alignment has no standard amino-acid residues")
    pooled = [count / pooled_total for count in pooled_counts]

    rows: list[dict[str, object]] = []
    for taxon in sorted(records):
        vector, observed = vectors[taxon]
        if observed == 0:
            distance: float | None = None
        else:
            midpoint = [(left + right) / 2.0 for left, right in zip(vector, pooled)]
            divergence = 0.5 * _kl_divergence(vector, midpoint) + 0.5 * _kl_divergence(
                pooled, midpoint
            )
            distance = math.sqrt(max(divergence, 0.0))
        row: dict[str, object] = {
            "taxon": taxon,
            "observed_standard_residues": observed,
            "gap_or_ambiguity": columns - observed,
            "gap_or_ambiguity_fraction": (columns - observed) / columns,
            "jensen_shannon_distance": distance,
        }
        row.update({f"frequency_{aa}": vector[index] for index, aa in enumerate(STANDARD_AA)})
        rows.append(row)
    return rows


def benjamini_hochberg(values: Sequence[float | None]) -> list[float | None]:
    """Return Benjamini-Hochberg adjusted p-values in original order."""

    indexed: list[tuple[int, float]] = []
    for index, value in enumerate(values):
        if value is None:
            continue
        numeric = float(value)
        if not 0.0 <= numeric <= 1.0:
            raise ValueError(f"p-value outside [0, 1]: {numeric}")
        indexed.append((index, numeric))
    result: list[float | None] = [None] * len(values)
    if not indexed:
        return result
    ordered = sorted(indexed, key=lambda item: (item[1], item[0]))
    count = len(ordered)
    running = 1.0
    adjusted_by_index: dict[int, float] = {}
    for rank_from_end in range(count - 1, -1, -1):
        original_index, p_value = ordered[rank_from_end]
        rank = rank_from_end + 1
        running = min(running, p_value * count / rank)
        adjusted_by_index[original_index] = min(running, 1.0)
    for index, adjusted in adjusted_by_index.items():
        result[index] = adjusted
    return result


def _linear_quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot calculate a quantile of no values")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def terminal_branch_diagnostics(
    branch_lengths: Mapping[str, float],
) -> list[dict[str, object]]:
    """Flag terminal branches using the pre-registered robust threshold."""

    if not branch_lengths:
        raise ValueError("no terminal branch lengths supplied")
    values = [float(value) for value in branch_lengths.values()]
    if any(not math.isfinite(value) or value < 0.0 for value in values):
        raise ValueError("terminal branch lengths must be finite and non-negative")
    median = statistics.median(values)
    mad = statistics.median(abs(value - median) for value in values)
    if mad > 0.0:
        scaled_mad = 1.4826 * mad
        threshold = median + 3.0 * scaled_mad
        method = "MEDIAN_PLUS_3_SCALED_MAD"
    else:
        first_quartile = _linear_quantile(values, 0.25)
        third_quartile = _linear_quantile(values, 0.75)
        threshold = third_quartile + 1.5 * (third_quartile - first_quartile)
        scaled_mad = 0.0
        method = "TUKEY_UPPER_FENCE"
    return [
        {
            "taxon": taxon,
            "terminal_branch_length": float(branch_lengths[taxon]),
            "median": median,
            "mad": mad,
            "scaled_mad": scaled_mad,
            "long_branch_threshold": threshold,
            "threshold_method": method,
            "long_branch_flag": float(branch_lengths[taxon]) > threshold,
        }
        for taxon in sorted(branch_lengths)
    ]


def mask_alignment_by_occupancy(
    records: Mapping[str, str], *, minimum_present: int
) -> dict[str, object]:
    """Retain columns with standard residues in at least ``minimum_present`` taxa."""

    columns = _validate_aligned_records(records)
    if not isinstance(minimum_present, int) or not 1 <= minimum_present <= len(records):
        raise ValueError(
            f"minimum_present must be an integer in [1, {len(records)}]"
        )
    retained_zero_based = [
        column
        for column in range(columns)
        if sum(
            records[taxon][column].upper() in STANDARD_AA_SET for taxon in records
        )
        >= minimum_present
    ]
    if not retained_zero_based:
        raise ValueError("occupancy masking removed every alignment column")
    masked = {
        taxon: "".join(sequence[column] for column in retained_zero_based)
        for taxon, sequence in records.items()
    }
    return {
        "records": masked,
        "source_columns": columns,
        "retained_columns": len(retained_zero_based),
        "retained_source_columns_1based": [column + 1 for column in retained_zero_based],
        "minimum_present": minimum_present,
        "taxa": len(records),
    }


def mask_partitioned_alignment(
    records: Mapping[str, str],
    partitions: Sequence[tuple[str, int, int]],
    *,
    minimum_present: int,
) -> dict[str, object]:
    """Mask each 1-based inclusive partition and regenerate concatenated endpoints."""

    columns = _validate_aligned_records(records)
    if not partitions:
        raise ValueError("no partitions supplied")
    if not isinstance(minimum_present, int) or not 1 <= minimum_present <= len(records):
        raise ValueError(
            f"minimum_present must be an integer in [1, {len(records)}]"
        )
    names: set[str] = set()
    expected_start = 1
    retained_global: list[int] = []
    output_partitions: list[dict[str, object]] = []
    output_start = 1
    for name, start, end in partitions:
        if not name or name in names:
            raise ValueError(f"duplicate or empty partition name: {name!r}")
        names.add(name)
        if start != expected_start or end < start or end > columns:
            raise ValueError(
                f"invalid/non-contiguous partition {name}: {start}-{end}; "
                f"expected start {expected_start} within 1-{columns}"
            )
        retained_here = [
            column
            for column in range(start - 1, end)
            if sum(
                records[taxon][column].upper() in STANDARD_AA_SET
                for taxon in records
            )
            >= minimum_present
        ]
        if not retained_here:
            raise ValueError(f"empty partition after occupancy masking: {name}")
        retained_global.extend(retained_here)
        retained_count = len(retained_here)
        output_end = output_start + retained_count - 1
        output_partitions.append(
            {
                "name": name,
                "start": output_start,
                "end": output_end,
                "retained": retained_count,
            }
        )
        output_start = output_end + 1
        expected_start = end + 1
    if expected_start != columns + 1:
        raise ValueError(
            f"partitions do not cover the full alignment: end={expected_start - 1}, "
            f"columns={columns}"
        )
    masked = {
        taxon: "".join(sequence[column] for column in retained_global)
        for taxon, sequence in records.items()
    }
    return {
        "records": masked,
        "source_columns": columns,
        "retained_columns": len(retained_global),
        "retained_source_columns_1based": [column + 1 for column in retained_global],
        "minimum_present": minimum_present,
        "taxa": len(records),
        "partitions": output_partitions,
    }


DAYHOFF6_GROUPS = ("AGPST", "DENQ", "HKR", "ILMV", "FWY", "C")
DAYHOFF6_MAP = {
    amino_acid: str(group_index)
    for group_index, group in enumerate(DAYHOFF6_GROUPS)
    for amino_acid in group
}


def dayhoff6_recode(records: Mapping[str, str]) -> dict[str, str]:
    """Apply the locked six-state Dayhoff recoding; all missingness becomes ``?``."""

    if not records:
        raise ValueError("no sequences supplied for Dayhoff-6 recoding")
    recoded: dict[str, str] = {}
    approved_missing = ALIGNED_MISSING
    for taxon, sequence in records.items():
        output: list[str] = []
        for symbol in sequence.upper():
            if symbol in DAYHOFF6_MAP:
                output.append(DAYHOFF6_MAP[symbol])
            elif symbol in approved_missing:
                output.append("?")
            else:
                raise ValueError(
                    f"unapproved symbol for Dayhoff-6 recoding in {taxon}: {symbol!r}"
                )
        recoded[taxon] = "".join(output)
    return recoded


STRONG_SH_ALRT = 80.0
STRONG_UFBOOT = 95.0


def canonical_side(side: Iterable[str], universe: Iterable[str]) -> tuple[str, ...]:
    """Return a root-invariant canonical representation of an unrooted split."""

    universe_set = frozenset(universe)
    side_set = frozenset(side)
    if not side_set or not side_set < universe_set:
        raise ValueError("split side must be a nonempty proper subset of the universe")
    complement = universe_set - side_set
    left = tuple(sorted(side_set))
    right = tuple(sorted(complement))
    if len(left) != len(right):
        return left if len(left) < len(right) else right
    return min(left, right)


def _support_pair(clade) -> tuple[float, float] | None:
    label = clade.name
    if label is None or "/" not in str(label):
        return None
    pieces = str(label).split("/")
    if len(pieces) != 2:
        return None
    try:
        return float(pieces[0]), float(pieces[1])
    except ValueError:
        return None


def parse_split_supports(newick: str) -> dict[tuple[str, ...], tuple[float, float]]:
    """Parse explicit SH-aLRT/UFBoot split labels from an unrooted Newick."""

    tree = Phylo.read(StringIO(newick.strip()), "newick")
    labels = [leaf.name for leaf in tree.get_terminals()]
    if any(label is None for label in labels) or len(labels) != len(set(labels)):
        raise ValueError("tree leaf labels are missing or duplicated")
    universe = frozenset(labels)
    supports: dict[tuple[str, ...], tuple[float, float]] = {}
    for clade in tree.get_nonterminals(order="preorder"):
        descendants = frozenset(leaf.name for leaf in clade.get_terminals())
        complement = universe - descendants
        if descendants == universe or len(descendants) < 2 or len(complement) < 2:
            continue
        support = _support_pair(clade)
        if support is None:
            continue
        key = canonical_side(descendants, universe)
        previous = supports.get(key)
        if previous is None or min(support) > min(previous):
            supports[key] = support
    return supports


def parse_supported_splits(
    newick: str,
    *,
    sh_threshold: float = STRONG_SH_ALRT,
    ufboot_threshold: float = STRONG_UFBOOT,
) -> dict[tuple[str, ...], tuple[float, float]]:
    """Return strong unrooted splits under the locked 80/95 rule."""

    return {
        split: support
        for split, support in parse_split_supports(newick).items()
        if support[0] >= sh_threshold and support[1] >= ufboot_threshold
    }


def reproducible_exact_alternatives(
    target_pair: frozenset[str],
    strong_splits_by_tree: Mapping[str, set[frozenset[str]]],
) -> dict[tuple[str, ...], list[str]]:
    """Find the same strong alternative two-leaf pairing in at least two trees."""

    if len(target_pair) != 2:
        raise ValueError("target pair must contain exactly two taxa")
    occurrences: dict[frozenset[str], list[str]] = {}
    for tree_name, splits in strong_splits_by_tree.items():
        for split in splits:
            if len(split) != 2 or split == target_pair:
                continue
            if len(split & target_pair) != 1:
                continue
            occurrences.setdefault(frozenset(split), []).append(tree_name)
    return {
        tuple(sorted(alternative)): sorted(tree_names)
        for alternative, tree_names in sorted(
            occurrences.items(), key=lambda item: tuple(sorted(item[0]))
        )
        if len(set(tree_names)) >= 2
    }


def adjudicate_m09d_relationship(
    *,
    rbh: bool,
    original_exact: Mapping[str, bool],
    new_exact: Mapping[str, bool | None],
    original_component: Mapping[str, bool],
    new_component: Mapping[str, bool | None],
    reproducible_alternative_count: int,
    direct_caveat: bool,
) -> str:
    """Apply the pre-registered conservative M09D classification hierarchy."""

    required_original = {"M09C_PRIMARY", "M09C_FULL"}
    required_new = {"R1", "R2", "R3", "R4"}
    if set(original_exact) != required_original:
        raise ValueError("original_exact must contain M09C_PRIMARY and M09C_FULL")
    if set(original_component) != required_original:
        raise ValueError(
            "original_component must contain M09C_PRIMARY and M09C_FULL"
        )
    if set(new_exact) != required_new or set(new_component) != required_new:
        raise ValueError("new evidence must contain exactly R1, R2, R3 and R4")
    if reproducible_alternative_count < 0:
        raise ValueError("reproducible_alternative_count cannot be negative")

    if reproducible_alternative_count > 0:
        return "CONTRADICTED"

    original_exact_both = all(original_exact.values())
    exact_positive = {name for name, value in new_exact.items() if value is True}
    all_new_computed = all(value is not None for value in new_exact.values())
    has_r1 = "R1" in exact_positive
    has_occupancy = bool(exact_positive & {"R2", "R3"})

    if (
        rbh
        and original_exact_both
        and all_new_computed
        and len(exact_positive) >= 3
        and has_r1
        and has_occupancy
        and not direct_caveat
    ):
        return "ROBUST_COMPUTATIONAL_CANDIDATE"

    if (
        rbh
        and original_exact_both
        and len(exact_positive) >= 2
        and (has_r1 or has_occupancy)
    ):
        return "SUPPORTED_WITH_MODEL_CAVEAT"

    component_positive = sum(value is True for value in new_component.values())
    if all(original_component.values()) and component_positive >= 2:
        return "POSSIBLE_CO_ORTHOLOGY"
    return "UNRESOLVED"
