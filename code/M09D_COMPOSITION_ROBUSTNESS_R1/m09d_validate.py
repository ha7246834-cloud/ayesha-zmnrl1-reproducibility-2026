#!/usr/bin/env python3
"""Producer-independent validation for an M09D staging tree or immutable release."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Iterable

from Bio import Phylo
from PIL import Image


EXPECTED_MASTER_SHA256 = "0e9d08a9d3c66cfc37cf58eab28477bd7f90e6bdaa787cd1068396556b449f37"
EXPECTED_MAIZE = {f"ZmNRG2_{index:03d}" for index in range(1, 24)}
EXPECTED_NO_RBH = {
    "ZmNRG2_004",
    "ZmNRG2_005",
    "ZmNRG2_011",
    "ZmNRG2_013",
    "ZmNRG2_016",
    "ZmNRG2_020",
    "ZmNRG2_022",
}
ALLOWED_CLASSES = {
    "ROBUST_COMPUTATIONAL_CANDIDATE",
    "SUPPORTED_WITH_MODEL_CAVEAT",
    "POSSIBLE_CO_ORTHOLOGY",
    "UNRESOLVED",
    "CONTRADICTED",
}
REQUIRED_ROOT_FILES = {
    "M09D_QC.json",
    "M09D_SOURCE_LOCK.tsv",
    "M09D_COMPOSITION_DIAGNOSTICS.tsv",
    "M09D_ALIGNMENT_SENSITIVITY.tsv",
    "M09D_MODEL_COMPARISON.tsv",
    "M09D_ALL_23_FINAL_REVIEW.tsv",
    "M09D_SEVEN_NO_RBH_CASES.tsv",
    "M09D_ORTHOLOGY_EVIDENCE_MATRIX.tsv",
    "M09D_SCIENTIFIC_LIMITATIONS.md",
    "M09D_MANUSCRIPT_METHODS_RESULTS_DRAFT.md",
    "M09D_MANUSCRIPT_LEGENDS.md",
    "M09D_SUPPLEMENTARY_INDEX.tsv",
    "M09D_REPRODUCIBILITY.md",
    "M09D_EXECUTION_SUMMARY.md",
    "M09D_SOFTWARE_AND_PARAMETERS.tsv",
    "M09D_RESOURCE_AND_PERFORMANCE_SUMMARY.tsv",
    "RUN_STATE.json",
    "RESUME_COMMAND.txt",
    "SOURCE_INTEGRITY_AUDIT.json",
}


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    pieces: list[str] = []

    def finish() -> None:
        nonlocal pieces
        if identifier is not None:
            sequence = "".join(pieces)
            if not sequence:
                raise ValueError(f"empty FASTA sequence: {identifier}")
            records[identifier] = sequence

    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            token = line[1:].split()[0] if line[1:].split() else ""
            if not token or token in records or token == identifier:
                raise ValueError(f"missing or duplicate FASTA identifier at line {line_number}")
            identifier, pieces = token, []
        else:
            if identifier is None:
                raise ValueError("FASTA sequence precedes first identifier")
            pieces.append(line)
    finish()
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def tree_leaves(path: str | Path) -> set[str]:
    tree = Phylo.read(str(path), "newick")
    labels = [leaf.name for leaf in tree.get_terminals()]
    if any(label is None for label in labels) or len(labels) != len(set(labels)):
        raise ValueError(f"missing or duplicate tree leaf labels: {path}")
    return set(labels)


def verify_complete_manifest(root: str | Path, manifest: str | Path) -> int:
    release = Path(root).resolve()
    manifest_path = Path(manifest).resolve()
    try:
        manifest_path.relative_to(release)
    except ValueError as exc:
        raise ValueError("manifest is outside release") from exc
    listed: list[str] = []
    for line_number, raw in enumerate(manifest_path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        pieces = raw.split(maxsplit=1)
        if len(pieces) != 2 or len(pieces[0]) != 64:
            raise ValueError(f"malformed manifest row {line_number}")
        relative = pieces[1].lstrip(" *")
        candidate = (release / relative).resolve()
        try:
            candidate.relative_to(release)
        except ValueError as exc:
            raise ValueError(f"unsafe manifest path at row {line_number}: {relative}") from exc
        if not candidate.is_file():
            raise ValueError(f"manifest file is absent: {relative}")
        if sha256_file(candidate) != pieces[0].lower():
            raise ValueError(f"manifest checksum mismatch: {relative}")
        listed.append(relative)
    actual = {
        path.relative_to(release).as_posix()
        for path in release.rglob("*")
        if path.is_file() and path.resolve() != manifest_path
    }
    if not listed or len(listed) != len(set(listed)):
        raise ValueError("manifest is empty or has duplicate rows")
    if set(listed) != actual:
        missing = sorted(actual - set(listed))
        extra = sorted(set(listed) - actual)
        raise ValueError(f"manifest coverage mismatch: unlisted={missing}, nonexistent={extra}")
    return len(listed)


def _verify_upstream_manifest(root: Path, manifest: Path, expected_count: int) -> int:
    listed: list[str] = []
    for line_number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        pieces = raw.split(maxsplit=1)
        if len(pieces) != 2 or len(pieces[0]) != 64:
            raise ValueError(f"malformed upstream manifest row {line_number}: {manifest}")
        relative = pieces[1].lstrip(" *")
        candidate = (root / relative).resolve()
        try:
            candidate.relative_to(root.resolve())
        except ValueError as exc:
            raise ValueError(f"unsafe upstream manifest path: {relative}") from exc
        if not candidate.is_file() or sha256_file(candidate) != pieces[0].lower():
            raise ValueError(f"upstream manifest verification failed: {relative}")
        listed.append(relative)
    if len(listed) != expected_count or len(set(listed)) != expected_count:
        raise ValueError(
            f"upstream manifest count differs: {manifest}: {len(listed)} != {expected_count}"
        )
    return len(listed)


def _require_files(release: Path, relative_paths: Iterable[str]) -> None:
    absent = [relative for relative in relative_paths if not (release / relative).is_file()]
    if absent:
        raise ValueError(f"mandatory files are absent: {absent}")


def _verify_source_lock(release: Path) -> dict[str, object]:
    current = read_tsv(release / "M09D_SOURCE_LOCK.tsv")
    initial = read_tsv(release / "M09D_INITIAL_SOURCE_LOCK.tsv")
    if current != initial:
        raise ValueError("M09D source lock differs from its initial frozen copy")
    if not current:
        raise ValueError("M09D source lock is empty")
    role_paths: dict[str, Path] = {}
    for row in current:
        path = Path(row["path"])
        target = path if path.is_absolute() else release / path
        if not target.is_file() or sha256_file(target) != row["sha256"]:
            raise ValueError(f"source-lock verification failed: {row['role']}")
        if target.stat().st_size != int(row["bytes"]):
            raise ValueError(f"source-lock byte count failed: {row['role']}")
        role_paths[row["role"]] = target
    required_roles = {"M09C_MANIFEST", "M09B_MANIFEST", "CANONICAL_MASTER"}
    if not required_roles <= set(role_paths):
        raise ValueError("mandatory upstream roles are absent from source lock")
    m09c_count = _verify_upstream_manifest(
        role_paths["M09C_MANIFEST"].parent, role_paths["M09C_MANIFEST"], 120
    )
    m09b_count = _verify_upstream_manifest(
        role_paths["M09B_MANIFEST"].parent, role_paths["M09B_MANIFEST"], 61
    )
    return {
        "rows": len(current),
        "status": "PASS",
        "m09c_manifest_entries": m09c_count,
        "m09b_manifest_entries": m09b_count,
    }


def _verify_alignments(release: Path) -> dict[str, object]:
    input_dir = release / "PHYLOGENY/INPUTS"
    alignment_dir = release / "PHYLOGENY/ALIGNMENTS"
    dual = parse_fasta(input_dir / "M09D_M09C_DUAL_DOMAIN_40_UNCHANGED.faa")
    full = parse_fasta(input_dir / "M09D_M09C_FULL_LENGTH_40_UNCHANGED.faa")
    if set(dual) != set(full) or len(dual) != 40:
        raise ValueError("locked input alignments do not share exactly 40 taxa")
    maize = {taxon for taxon in dual if taxon.startswith("ZmNRG2_")}
    sorghum = set(dual) - maize
    if maize != EXPECTED_MAIZE or len(sorghum) != 17:
        raise ValueError("locked alignment taxon census is not 23 maize plus 17 Sorghum")
    if {len(sequence) for sequence in dual.values()} != {506}:
        raise ValueError("locked dual-domain alignment is not 506 columns")
    if {len(sequence) for sequence in full.values()} != {1611}:
        raise ValueError("locked full-length alignment is not 1611 columns")
    widths = {}
    expected_widths = {
        "R1": ("M09D_R1_DUAL_DOMAIN_UNCHANGED.faa", 506),
        "R2": ("M09D_R2_DUAL_DOMAIN_OCCUPANCY20.faa", 363),
        "R3": ("M09D_R3_FULL_LENGTH_OCCUPANCY20.faa", 735),
        "R4": ("M09D_R4_DUAL_DOMAIN_DAYHOFF6.faa", 506),
    }
    for analysis, (name, expected_width) in expected_widths.items():
        records = parse_fasta(alignment_dir / name)
        if set(records) != set(dual) or {len(sequence) for sequence in records.values()} != {expected_width}:
            raise ValueError(f"{analysis} transformed alignment dimensions/taxa are invalid")
        widths[analysis] = expected_width
        if analysis == "R4" and any(set(sequence) - set("012345?-") for sequence in records.values()):
            raise ValueError("R4 recoded alignment contains an invalid state")
    return {"taxa": 40, "maize": 23, "sorghum": 17, "widths": widths, "labels": sorted(dual)}


def _verify_trees(release: Path, expected_labels: set[str]) -> dict[str, object]:
    configuration = {
        "R1": "R1_SITE_HETEROGENEOUS/M09D_R1",
        "R2": "R2_DUAL_OCCUPANCY20/M09D_R2",
        "R3": "R3_FULL_OCCUPANCY20/M09D_R3",
        "R4": "R4_DAYHOFF6/M09D_R4",
    }
    output: dict[str, object] = {}
    for analysis, relative_prefix in configuration.items():
        prefix = release / "PHYLOGENY/TREES" / relative_prefix
        required = [prefix.with_suffix(suffix) for suffix in (".treefile", ".contree", ".iqtree", ".log", ".ufboot", ".runtrees")]
        if any(not path.is_file() or path.stat().st_size == 0 for path in required):
            raise ValueError(f"{analysis} IQ-TREE artifact set is incomplete")
        if tree_leaves(required[0]) != expected_labels or tree_leaves(required[1]) != expected_labels:
            raise ValueError(f"{analysis} tree leaves differ from the locked 40 taxa")
        boot_count = sum(bool(line.strip()) for line in required[4].read_text(encoding="utf-8").splitlines())
        run_count = sum(bool(line.strip()) for line in required[5].read_text(encoding="utf-8").splitlines())
        if boot_count != 1000 or run_count != 5:
            raise ValueError(f"{analysis} replicate/search counts are invalid: {boot_count}/{run_count}")
        metrics = json.loads((prefix.parent / "METRICS.json").read_text(encoding="utf-8"))
        if not str(metrics.get("status", "")).startswith("PASS_COMPLETE"):
            raise ValueError(f"{analysis} metrics status is not complete")
        if metrics.get("treefile_sha256") != sha256_file(required[0]):
            raise ValueError(f"{analysis} tree hash differs from metrics")
        output[analysis] = {
            "treefile_sha256": sha256_file(required[0]),
            "ufboot_trees": boot_count,
            "independent_runs": run_count,
        }
    return output


def _verify_tables(release: Path) -> dict[str, object]:
    all23 = read_tsv(release / "M09D_ALL_23_FINAL_REVIEW.tsv")
    seven = read_tsv(release / "M09D_SEVEN_NO_RBH_CASES.tsv")
    evidence = read_tsv(release / "M09D_ORTHOLOGY_EVIDENCE_MATRIX.tsv")
    if len(all23) != 23 or {row["maize_family_alias"] for row in all23} != EXPECTED_MAIZE:
        raise ValueError("final review does not account for exactly 23 maize genes")
    if len(seven) != 7 or {row["maize_family_alias"] for row in seven} != EXPECTED_NO_RBH:
        raise ValueError("NO_RBH table does not contain the predetermined seven maize genes")
    if any(row["reciprocal_best_hit_status"] != "NO" for row in seven):
        raise ValueError("NO_RBH table contains a reciprocal-best-hit row")
    if len(evidence) != 24 or len({(row["maize_family_alias"], row["sorghum_candidate_id"]) for row in evidence}) != 24:
        raise ValueError("orthology matrix does not contain 24 unique pair hypotheses")
    source_evidence = release / "PHYLOGENY/SOURCE_EVIDENCE"
    inherited_candidates = read_tsv(source_evidence / "M09C_CANDIDATE_PAIR_ADJUDICATION.tsv")
    inherited_no_rbh = read_tsv(source_evidence / "M09C_NO_RBH_SEVEN_CASES.tsv")
    family_anchors = read_tsv(source_evidence / "M09B_NRG2_FAMILY_FAMILY_ANCHORS.tsv")
    if len(inherited_candidates) != 16 or len(inherited_no_rbh) != 7:
        raise ValueError("inherited M09C candidate/NO_RBH counts are not 16/7")
    anchor_pairs = {
        (row.get("maize_family_alias", ""), row.get("sorghum_gene_id", ""))
        for row in family_anchors
    }
    if len(family_anchors) != 26 or len(anchor_pairs) != 26:
        raise ValueError("inherited M09B family anchors are not 26 unique pairs")
    classes = [row["final_computational_classification"] for row in all23]
    if set(classes) - ALLOWED_CLASSES:
        raise ValueError("final review contains an unknown classification")
    composition = read_tsv(release / "M09D_COMPOSITION_DIAGNOSTICS.tsv")
    if len(composition) != 80:
        raise ValueError(f"composition diagnostics have {len(composition)} rows, expected 80")
    sensitivity = read_tsv(release / "M09D_ALIGNMENT_SENSITIVITY.tsv")
    retained = {row["analysis"]: int(row["retained_columns"]) for row in sensitivity}
    if retained != {"R1": 506, "R2": 363, "R3": 735, "R4": 506}:
        raise ValueError(f"alignment-sensitivity counts differ: {retained}")
    return {
        "all23": len(all23),
        "seven_no_rbh": len(seven),
        "pair_hypotheses": len(evidence),
        "inherited_candidate_pairs": len(inherited_candidates),
        "inherited_no_rbh_cases": len(inherited_no_rbh),
        "inherited_family_anchor_pairs": len(anchor_pairs),
        "classification_counts": {name: classes.count(name) for name in sorted(set(classes))},
    }


def _verify_figures(release: Path) -> dict[str, object]:
    figure_dir = release / "PHYLOGENY/FIGURES"
    qc = json.loads((figure_dir / "M09D_FIGURE_QC.json").read_text(encoding="utf-8"))
    attestation = json.loads((figure_dir / "M09D_VISUAL_QA_ATTESTATION.json").read_text(encoding="utf-8"))
    if attestation.get("status") != "PASS_VISUAL_QA":
        raise ValueError("visual-QA attestation is absent or not PASS")
    prefixes = [f"M09D_{analysis}_UNROOTED_ROBUSTNESS_TREE" for analysis in ("R1", "R2", "R3", "R4")]
    prefixes.append("M09D_ALL_23_ROBUSTNESS_MATRIX")
    tree_source = read_tsv(figure_dir / "M09D_FIGURE_SOURCE.tsv")
    matrix_source = read_tsv(figure_dir / "M09D_ROBUSTNESS_MATRIX_SOURCE.tsv")
    if len(tree_source) != 160 or len(matrix_source) != 23:
        raise ValueError("figure source-table row counts are not 160 and 23")
    by_analysis: dict[str, list[dict[str, str]]] = {}
    for row in tree_source:
        by_analysis.setdefault(row["analysis"], []).append(row)
    sequence_map = {
        row["tree_label"]: row
        for row in read_tsv(release / "PHYLOGENY/SOURCE_EVIDENCE/M09C_SEQUENCE_MAP.tsv")
    }
    tree_locations = {
        "R1": release / "PHYLOGENY/TREES/R1_SITE_HETEROGENEOUS/M09D_R1.treefile",
        "R2": release / "PHYLOGENY/TREES/R2_DUAL_OCCUPANCY20/M09D_R2.treefile",
        "R3": release / "PHYLOGENY/TREES/R3_FULL_OCCUPANCY20/M09D_R3.treefile",
        "R4": release / "PHYLOGENY/TREES/R4_DAYHOFF6/M09D_R4.treefile",
    }
    for analysis, tree_path in tree_locations.items():
        rows = by_analysis.get(analysis, [])
        labels = {row["tree_label"] for row in rows}
        if len(rows) != 40 or labels != tree_leaves(tree_path):
            raise ValueError(f"{analysis} figure source labels do not match its tree")
        expected_hash = sha256_file(tree_path)
        if {row["tree_sha256"] for row in rows} != {expected_hash}:
            raise ValueError(f"{analysis} figure source tree hash is inconsistent")
        for row in rows:
            source = sequence_map.get(row["tree_label"])
            if source is None or any(
                row[field] != source[field]
                for field in ("species", "family_alias", "gene_id", "protein_id", "source_fasta_id", "length_aa", "sequence_sha256")
            ):
                raise ValueError(f"{analysis} figure source accession mapping is inconsistent")
        svg_text = (figure_dir / f"M09D_{analysis}_UNROOTED_ROBUSTNESS_TREE.svg").read_text(encoding="utf-8")
        if any(label not in svg_text for label in labels):
            raise ValueError(f"{analysis} SVG is missing a source label")
    if {row["maize_family_alias"] for row in matrix_source} != EXPECTED_MAIZE:
        raise ValueError("robustness-matrix source does not account for 23 maize aliases")
    matrix_svg = (figure_dir / "M09D_ALL_23_ROBUSTNESS_MATRIX.svg").read_text(encoding="utf-8")
    if any(alias not in matrix_svg for alias in EXPECTED_MAIZE):
        raise ValueError("robustness-matrix SVG is missing a maize label")

    for prefix in prefixes:
        for suffix in (".svg", ".pdf", ".png"):
            path = figure_dir / f"{prefix}{suffix}"
            if not path.is_file() or path.stat().st_size <= 1000:
                raise ValueError(f"figure artifact is absent or too small: {path.name}")
        if not (figure_dir / f"{prefix}.pdf").read_bytes().startswith(b"%PDF"):
            raise ValueError(f"figure PDF signature is invalid: {prefix}")
        with Image.open(figure_dir / f"{prefix}.png") as image:
            if image.width < 3000 or image.height < 2000:
                raise ValueError(f"figure PNG is below the resolution floor: {prefix}")
            dpi = image.info.get("dpi", (0.0, 0.0))
            if len(dpi) < 2 or min(float(dpi[0]), float(dpi[1])) < 590.0:
                raise ValueError(f"figure PNG does not report approximately 600 dpi: {prefix}")
    return {"figure_sets": len(prefixes), "artifacts": len(prefixes) * 3, "programmatic_qc": qc.get("status")}


def validate_release(project_root: Path, release: Path, require_manifest: bool) -> dict[str, object]:
    release = release.resolve()
    project_root = project_root.resolve()
    if not release.is_dir():
        raise FileNotFoundError(f"release/staging directory is absent: {release}")
    _require_files(release, REQUIRED_ROOT_FILES)
    master = project_root / "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/MASTER_CURRENT.docx"
    master_hash = sha256_file(master)
    if master_hash != EXPECTED_MASTER_SHA256:
        raise ValueError(f"canonical master checksum mismatch: {master_hash}")
    source_lock = _verify_source_lock(release)
    alignments = _verify_alignments(release)
    trees = _verify_trees(release, set(alignments["labels"]))
    tables = _verify_tables(release)
    figures = _verify_figures(release)
    tests = json.loads((release / "PHYLOGENY/TEST_REPORTS/M09D_TEST_RESULTS.json").read_text(encoding="utf-8"))
    if tests.get("status") != "PASS" or tests["m09d"]["tests"] < 1 or tests["m09c_inherited"]["tests"] < 1:
        raise ValueError("fresh regression-test report is not PASS")
    qc = json.loads((release / "M09D_QC.json").read_text(encoding="utf-8"))
    if qc.get("status") not in {"PASS", "REVIEW"}:
        raise ValueError(f"M09D QC status is invalid for release: {qc.get('status')}")
    manifest_entries: int | str = "NOT_REQUESTED"
    if require_manifest:
        manifest_entries = verify_complete_manifest(release, release / "SHA256SUMS.txt")
    return {
        "schema": "M09D_INDEPENDENT_VALIDATION_V1",
        "status": "PASS",
        "release_reported_status": qc["status"],
        "release": str(release),
        "canonical_master_sha256": master_hash,
        "source_lock": source_lock,
        "alignments": {key: value for key, value in alignments.items() if key != "labels"},
        "trees": trees,
        "tables": tables,
        "figures": figures,
        "tests": {"m09d": tests["m09d"]["tests"], "m09c_inherited": tests["m09c_inherited"]["tests"]},
        "manifest_entries": manifest_entries,
    }


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.{os.getpid()}.tmp"
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    parser.add_argument("release", type=Path)
    parser.add_argument("--require-manifest", action="store_true")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    try:
        result = validate_release(arguments.project_root, arguments.release, arguments.require_manifest)
        exit_code = 0
    except Exception as exc:  # fail-closed CLI boundary
        result = {
            "schema": "M09D_INDEPENDENT_VALIDATION_V1",
            "status": "FAIL",
            "release": str(arguments.release.resolve()),
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        exit_code = 1
    if arguments.output:
        atomic_write_json(arguments.output, result)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
