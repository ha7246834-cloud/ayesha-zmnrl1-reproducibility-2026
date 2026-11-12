"""M09D composition-robustness orchestration."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from Bio import Phylo

import m09d_core as core


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_text(path: str | Path, text: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", dir=destination.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_write_json(path: str | Path, value: object) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def verify_manifest(root: str | Path, manifest: str | Path) -> int:
    """Verify a GNU-style SHA256 manifest while rejecting path escapes."""

    release = Path(root).resolve()
    manifest_path = Path(manifest)
    rows = []
    for line_number, raw in enumerate(
        manifest_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
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
            raise ValueError(f"manifest file missing: {relative}")
        observed = sha256_file(candidate)
        if observed != pieces[0].lower():
            raise ValueError(
                f"checksum mismatch for {relative}: expected {pieces[0]}, observed {observed}"
            )
        rows.append(relative)
    if not rows or len(rows) != len(set(rows)):
        raise ValueError("manifest is empty or contains duplicate paths")
    return len(rows)


def safe_staging_path(parent: str | Path, candidate: str | Path) -> Path:
    release_parent = Path(parent).resolve()
    staging = Path(candidate).resolve()
    if staging.parent != release_parent or not staging.name.startswith(".STAGING_"):
        raise ValueError("staging path must be an exclusive hidden child of release parent")
    return staging


def build_provenance(
    command: Sequence[str],
    inputs: Iterable[str | Path],
    executables: Iterable[str | Path],
    code_files: Iterable[str | Path],
) -> dict[str, object]:
    return {
        "schema": "M09D_PROVENANCE_V1",
        "command": list(command),
        "inputs": {str(Path(path)): sha256_file(path) for path in inputs},
        "executables": {str(Path(path)): sha256_file(path) for path in executables},
        "implementation": {str(Path(path)): sha256_file(path) for path in code_files},
    }


def provenance_matches(left: object, right: object) -> bool:
    return left == right


ANALYSIS_LOCK = {
    "R1": {"seqtype": "AA", "model": "LG+C20+F+G4", "seed": 20261108},
    "R2": {"seqtype": "AA", "model": "MFP+MERGE", "seed": 20261109},
    "R3": {"seqtype": "AA", "model": "MFP", "seed": 20261110},
    "R4": {"seqtype": "MORPH{6}", "model": "GTR+FO+G4", "seed": 20261111},
}


EXPECTED_MASTER_SHA256 = "0e9d08a9d3c66cfc37cf58eab28477bd7f90e6bdaa787cd1068396556b449f37"
M09C_RELATIVE = Path(
    "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/QA/"
    "M09C_INDEPENDENT_ORTHOLOGY_PHYLOGENY/RELEASE_20261107T133337Z_f5edc166"
)
M09B_RELATIVE = Path(
    "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/QA/"
    "M09B_INTER_SPECIES_DUAL_SYNTENY/RELEASE_20261106T075710Z_c0a4e014"
)
RELEASE_PARENT_RELATIVE = Path(
    "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/QA/M09D_COMPOSITION_ROBUSTNESS"
)
IQTREE3 = Path("/home/hammadse/miniconda3/envs/wrky_q1/bin/iqtree3")
PYTHON = Path("/home/hammadse/miniconda3/envs/wrky_q1/bin/python")
TIME_EXECUTABLE = Path("/usr/bin/time")
TOOL_DIR = Path(__file__).resolve().parent


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_optional_probability(value: str) -> float | None:
    """Normalize IQ-TREE missing/non-finite probability tokens to ``None``."""

    token = value.strip()
    if not token or token.upper() == "NA":
        return None
    probability = float(token)
    if not math.isfinite(probability):
        return None
    if not 0.0 <= probability <= 1.0:
        raise ValueError(f"probability outside [0, 1]: {probability}")
    return probability


def write_tsv(
    path: str | Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]
) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", dir=destination.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_fasta(path: str | Path, records: Mapping[str, str]) -> None:
    lines: list[str] = []
    for identifier, sequence in records.items():
        lines.append(f">{identifier}")
        lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    atomic_write_text(path, "\n".join(lines) + "\n")


def alignment_column_count(path: str | Path) -> int:
    """Return aligned FASTA width without imposing protein/morphology alphabet."""

    lengths: list[int] = []
    current = False
    pieces: list[str] = []
    identifiers: set[str] = set()

    def finish() -> None:
        nonlocal pieces
        if current:
            if not pieces:
                raise ValueError("empty aligned FASTA record")
            lengths.append(len("".join(pieces)))

    for line_number, raw in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            fields = line[1:].split()
            if not fields or fields[0] in identifiers:
                raise ValueError(
                    f"missing or duplicate aligned FASTA identifier at line {line_number}"
                )
            identifiers.add(fields[0])
            current = True
            pieces = []
        else:
            if not current:
                raise ValueError("aligned FASTA sequence before first header")
            pieces.append(line)
    finish()
    if not lengths:
        raise ValueError("empty aligned FASTA")
    if len(set(lengths)) != 1:
        raise ValueError(f"unequal aligned FASTA row lengths: {sorted(set(lengths))}")
    return lengths[0]


def _tree_leaf_set(path: Path) -> set[str]:
    tree = Phylo.read(str(path), "newick")
    labels = [leaf.name for leaf in tree.get_terminals()]
    if any(label is None for label in labels) or len(labels) != len(set(labels)):
        raise ValueError(f"missing or duplicate tree labels: {path}")
    return set(labels)


def _tree_terminal_lengths(path: Path) -> dict[str, float]:
    tree = Phylo.read(str(path), "newick")
    output: dict[str, float] = {}
    for leaf in tree.get_terminals():
        if leaf.name is None or leaf.name in output or leaf.branch_length is None:
            raise ValueError(f"invalid terminal branch in {path}")
        output[leaf.name] = float(leaf.branch_length)
    return output


def source_integrity_audit(project_root: Path) -> dict[str, object]:
    root = project_root.resolve()
    m09c = root / M09C_RELATIVE
    m09b = root / M09B_RELATIVE
    master = root / "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/MASTER_CURRENT.docx"
    m09c_count = verify_manifest(m09c, m09c / "SHA256SUMS")
    m09b_count = verify_manifest(m09b, m09b / "SHA256SUMS")
    if m09c_count != 120:
        raise ValueError(f"M09C manifest count is {m09c_count}, expected 120")
    if m09b_count != 61:
        raise ValueError(f"M09B manifest count is {m09b_count}, expected 61")
    master_hash = sha256_file(master)
    if master_hash != EXPECTED_MASTER_SHA256:
        raise ValueError(f"canonical master checksum mismatch: {master_hash}")

    dual = core.parse_fasta(
        m09c / "PHYLOGENY/ALIGNMENTS/M09C_DUAL_DOMAIN_40_CONCATENATED.faa",
        aligned=True,
    )
    full = core.parse_fasta(
        m09c / "PHYLOGENY/ALIGNMENTS/M09C_FULL_LENGTH_40_MAFFT.faa",
        aligned=True,
    )
    taxon_counts = core.validate_taxa(dual)
    if set(dual) != set(full):
        raise ValueError("M09C dual-domain and full-length taxon sets differ")
    if {len(value) for value in dual.values()} != {506}:
        raise ValueError("M09C dual-domain alignment is not exactly 506 columns")
    if {len(value) for value in full.values()} != {1611}:
        raise ValueError("M09C full-length alignment is not exactly 1611 columns")

    primary_tree = m09c / "PHYLOGENY/TREES/PRIMARY/M09C_PRIMARY.treefile"
    full_tree = m09c / "PHYLOGENY/TREES/FULL_LENGTH/M09C_FULL_LENGTH.treefile"
    if _tree_leaf_set(primary_tree) != set(dual) or _tree_leaf_set(full_tree) != set(dual):
        raise ValueError("M09C tree leaves do not match the locked 40 taxa")

    sequence_map = read_tsv(m09c / "PHYLOGENY/INPUTS/M09C_SEQUENCE_MAP.tsv")
    all_23 = read_tsv(m09c / "M09C_ALL_23_ORTHOLOGY_REVIEW.tsv")
    candidates = read_tsv(m09c / "M09C_CANDIDATE_PAIR_ADJUDICATION.tsv")
    no_rbh = read_tsv(m09c / "M09C_NO_RBH_SEVEN_CASES.tsv")
    family_anchors = read_tsv(
        m09c / "PHYLOGENY/SOURCE_EVIDENCE/M09B_NRG2_FAMILY_FAMILY_ANCHORS.tsv"
    )
    if len(sequence_map) != 40 or {row["tree_label"] for row in sequence_map} != set(dual):
        raise ValueError("M09C sequence map does not reconcile to 40 tree labels")
    if len(all_23) != 23 or len(candidates) != 16 or len(no_rbh) != 7:
        raise ValueError(
            f"M09C table counts differ: all23={len(all_23)}, "
            f"candidates={len(candidates)}, no_rbh={len(no_rbh)}"
        )
    unique_anchor_pairs = {
        (row.get("maize_family_alias", ""), row.get("sorghum_gene_id", ""))
        for row in family_anchors
    }
    if len(family_anchors) != 26 and len(unique_anchor_pairs) != 26:
        raise ValueError(f"M09B family anchor count is not 26: {len(family_anchors)}")
    qc = json.loads((m09c / "M09C_QC.json").read_text(encoding="utf-8"))
    if qc.get("status") != "REVIEW":
        raise ValueError(f"unexpected authoritative M09C status: {qc.get('status')}")
    return {
        "status": "PASS",
        "verified_utc": utc_now(),
        "m09c_release": str(m09c),
        "m09c_manifest_entries": m09c_count,
        "m09c_manifest_sha256": sha256_file(m09c / "SHA256SUMS"),
        "m09b_release": str(m09b),
        "m09b_manifest_entries": m09b_count,
        "m09b_manifest_sha256": sha256_file(m09b / "SHA256SUMS"),
        "canonical_master_sha256": master_hash,
        "taxa": taxon_counts,
        "dual_domain_columns": 506,
        "full_length_columns": 1611,
        "candidate_pairs": 16,
        "no_rbh_cases": 7,
        "family_anchor_rows": len(family_anchors),
    }


def _new_staging(project_root: Path) -> tuple[Path, Path]:
    parent = project_root / RELEASE_PARENT_RELATIVE
    parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = uuid.uuid4().hex[:8]
    token = f"{stamp}_{suffix}"
    staging = safe_staging_path(parent, parent / f".STAGING_{token}")
    final = parent.resolve() / f"RELEASE_{token}"
    if staging.exists() or final.exists():
        raise FileExistsError("new M09D staging or release path already exists")
    staging.mkdir()
    atomic_write_json(
        staging / "RUN_STATE.json",
        {
            "schema": "M09D_RUN_STATE_V1",
            "created_utc": utc_now(),
            "project_root": str(project_root.resolve()),
            "staging": str(staging),
            "final_release": str(final),
            "status": "STAGING_CREATED",
        },
    )
    atomic_write_text(
        staging / "RESUME_COMMAND.txt",
        f"bash {TOOL_DIR / 'MASTER_RUN_M09D.sh'} {project_root.resolve()} --resume {staging}\n",
    )
    return staging, final


def _copy_locked_inputs(project_root: Path, staging: Path) -> dict[str, Path]:
    m09c = project_root / M09C_RELATIVE
    inputs = staging / "PHYLOGENY/INPUTS"
    baseline = staging / "PHYLOGENY/BASELINES"
    evidence = staging / "PHYLOGENY/SOURCE_EVIDENCE"
    inputs.mkdir(parents=True, exist_ok=True)
    baseline.mkdir(parents=True, exist_ok=True)
    evidence.mkdir(parents=True, exist_ok=True)
    mapping = {
        "dual": (
            m09c / "PHYLOGENY/ALIGNMENTS/M09C_DUAL_DOMAIN_40_CONCATENATED.faa",
            inputs / "M09D_M09C_DUAL_DOMAIN_40_UNCHANGED.faa",
        ),
        "full": (
            m09c / "PHYLOGENY/ALIGNMENTS/M09C_FULL_LENGTH_40_MAFFT.faa",
            inputs / "M09D_M09C_FULL_LENGTH_40_UNCHANGED.faa",
        ),
        "partitions": (
            m09c / "PHYLOGENY/ALIGNMENTS/M09C_DUAL_DOMAIN_PARTITIONS.nex",
            inputs / "M09D_M09C_DUAL_DOMAIN_PARTITIONS_UNCHANGED.nex",
        ),
        "sequence_map": (
            m09c / "PHYLOGENY/INPUTS/M09C_SEQUENCE_MAP.tsv",
            evidence / "M09C_SEQUENCE_MAP.tsv",
        ),
        "domain_extraction": (
            m09c / "PHYLOGENY/INPUTS/M09C_DOMAIN_EXTRACTION.tsv",
            evidence / "M09C_DOMAIN_EXTRACTION.tsv",
        ),
        "primary_tree": (
            m09c / "PHYLOGENY/TREES/PRIMARY/M09C_PRIMARY.treefile",
            baseline / "M09C_PRIMARY.treefile",
        ),
        "full_tree": (
            m09c / "PHYLOGENY/TREES/FULL_LENGTH/M09C_FULL_LENGTH.treefile",
            baseline / "M09C_FULL_LENGTH.treefile",
        ),
        "primary_composition": (
            m09c / "PHYLOGENY/TREES/PRIMARY/M09C_PRIMARY_COMPOSITION_TEST.tsv",
            baseline / "M09C_PRIMARY_COMPOSITION_TEST.tsv",
        ),
        "full_composition": (
            m09c / "PHYLOGENY/TREES/FULL_LENGTH/M09C_FULL_LENGTH_COMPOSITION_TEST.tsv",
            baseline / "M09C_FULL_LENGTH_COMPOSITION_TEST.tsv",
        ),
        "all23": (
            m09c / "M09C_ALL_23_ORTHOLOGY_REVIEW.tsv",
            evidence / "M09C_ALL_23_ORTHOLOGY_REVIEW.tsv",
        ),
        "candidates": (
            m09c / "M09C_CANDIDATE_PAIR_ADJUDICATION.tsv",
            evidence / "M09C_CANDIDATE_PAIR_ADJUDICATION.tsv",
        ),
        "no_rbh": (
            m09c / "M09C_NO_RBH_SEVEN_CASES.tsv",
            evidence / "M09C_NO_RBH_SEVEN_CASES.tsv",
        ),
        "hypotheses": (
            m09c / "M09C_HYPOTHESIS_RECONSTRUCTION_AND_DUPLICATE_AUDIT.tsv",
            evidence / "M09C_HYPOTHESIS_RECONSTRUCTION_AND_DUPLICATE_AUDIT.tsv",
        ),
        "family_anchors": (
            m09c / "PHYLOGENY/SOURCE_EVIDENCE/M09B_NRG2_FAMILY_FAMILY_ANCHORS.tsv",
            evidence / "M09B_NRG2_FAMILY_FAMILY_ANCHORS.tsv",
        ),
    }
    for source, destination in mapping.values():
        shutil.copy2(source, destination)
        if sha256_file(source) != sha256_file(destination):
            raise ValueError(f"copy checksum mismatch: {source}")
    return {name: destination for name, (_source, destination) in mapping.items()}


def build_iqtree_command(
    executable: str | Path,
    analysis_id: str,
    alignment: Path,
    prefix: Path,
    *,
    partitions: Path | None = None,
) -> list[str]:
    """Build one of the four pre-registered IQ-TREE command profiles."""

    if analysis_id not in ANALYSIS_LOCK:
        raise ValueError(f"unknown analysis: {analysis_id}")
    if analysis_id == "R2" and partitions is None:
        raise ValueError("R2 requires the regenerated partition file")
    if analysis_id != "R2" and partitions is not None:
        raise ValueError(f"{analysis_id} does not accept a partition file")
    lock = ANALYSIS_LOCK[analysis_id]
    command = [str(executable), "-s", str(alignment)]
    if partitions is not None:
        command.extend(["-p", str(partitions)])
    command.extend(
        [
            "--seqtype",
            str(lock["seqtype"]),
            "-m",
            str(lock["model"]),
            "-B",
            "1000",
            "--alrt",
            "1000",
            "-bnni",
            "--boot-trees",
            "--runs",
            "5",
            "-T",
            "2",
            "--seed",
            str(lock["seed"]),
            "--prefix",
            str(prefix),
        ]
    )
    return command


def _write_source_lock(project_root: Path, staging: Path, copied: Mapping[str, Path]) -> None:
    m09c = project_root / M09C_RELATIVE
    m09b = project_root / M09B_RELATIVE
    master = project_root / "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/MASTER_CURRENT.docx"
    upstream = {
        "M09C_MANIFEST": m09c / "SHA256SUMS",
        "M09C_QC": m09c / "M09C_QC.json",
        "M09B_MANIFEST": m09b / "SHA256SUMS",
        "CANONICAL_MASTER": master,
        "M09D_CORE_CODE": TOOL_DIR / "m09d_core.py",
        "M09D_WORKFLOW_CODE": TOOL_DIR / "m09d_workflow.py",
        "M09D_DESIGN": project_root
        / "docs/superpowers/specs/2026-11-07-m09d-composition-robustness-design.md",
        "M09D_PLAN": project_root
        / "docs/superpowers/plans/2026-11-07-m09d-composition-robustness.md",
    }
    rows: list[dict[str, object]] = []
    for role, path in sorted(upstream.items()):
        rows.append(
            {
                "role": role,
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "lock_type": "UPSTREAM_READ_ONLY",
            }
        )
    for role, path in sorted(copied.items()):
        rows.append(
            {
                "role": f"COPIED_{role.upper()}",
                "path": str(path.relative_to(staging)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "lock_type": "STAGING_SOURCE_LOCKED_COPY",
            }
        )
    fields = ["role", "path", "sha256", "bytes", "lock_type"]
    write_tsv(staging / "M09D_SOURCE_LOCK.tsv", rows, fields)
    shutil.copy2(staging / "M09D_SOURCE_LOCK.tsv", staging / "M09D_INITIAL_SOURCE_LOCK.tsv")


def _prepare_transforms(staging: Path, copied: Mapping[str, Path]) -> dict[str, object]:
    alignments = staging / "PHYLOGENY/ALIGNMENTS"
    alignments.mkdir(parents=True, exist_ok=True)
    dual = core.parse_fasta(copied["dual"], aligned=True)
    full = core.parse_fasta(copied["full"], aligned=True)
    core.validate_taxa(dual)
    if set(dual) != set(full):
        raise ValueError("transform input taxon sets differ")

    partition_text = copied["partitions"].read_text(encoding="utf-8")
    if "PF04783 = 1-70" not in partition_text or "PF04782 = 71-506" not in partition_text:
        raise ValueError("locked M09C partition endpoints are not 1-70 and 71-506")

    masked_dual = core.mask_partitioned_alignment(
        dual,
        [("PF04783", 1, 70), ("PF04782", 71, 506)],
        minimum_present=20,
    )
    masked_full = core.mask_alignment_by_occupancy(full, minimum_present=20)
    dayhoff = core.dayhoff6_recode(dual)

    paths = {
        "R1": alignments / "M09D_R1_DUAL_DOMAIN_UNCHANGED.faa",
        "R2": alignments / "M09D_R2_DUAL_DOMAIN_OCCUPANCY20.faa",
        "R3": alignments / "M09D_R3_FULL_LENGTH_OCCUPANCY20.faa",
        "R4": alignments / "M09D_R4_DUAL_DOMAIN_DAYHOFF6.faa",
        "R2_PARTITIONS": alignments / "M09D_R2_MASKED_PARTITIONS.nex",
    }
    write_fasta(paths["R1"], dual)
    write_fasta(paths["R2"], masked_dual["records"])
    write_fasta(paths["R3"], masked_full["records"])
    write_fasta(paths["R4"], dayhoff)
    partition_rows = masked_dual["partitions"]
    atomic_write_text(
        paths["R2_PARTITIONS"],
        "#nexus\nbegin sets;\n"
        + "".join(
            f"  charset {row['name']} = {row['start']}-{row['end']};\n"
            for row in partition_rows
        )
        + "end;\n",
    )

    column_rows: list[dict[str, object]] = []
    for retained_index, source_column in enumerate(
        masked_dual["retained_source_columns_1based"], start=1
    ):
        column_rows.append(
            {
                "analysis": "R2",
                "retained_column_1based": retained_index,
                "source_column_1based": source_column,
                "source_partition": "PF04783" if source_column <= 70 else "PF04782",
            }
        )
    for retained_index, source_column in enumerate(
        masked_full["retained_source_columns_1based"], start=1
    ):
        column_rows.append(
            {
                "analysis": "R3",
                "retained_column_1based": retained_index,
                "source_column_1based": source_column,
                "source_partition": "FULL_LENGTH",
            }
        )
    write_tsv(
        staging / "PHYLOGENY/ALIGNMENTS/M09D_RETAINED_COLUMN_MAP.tsv",
        column_rows,
        [
            "analysis",
            "retained_column_1based",
            "source_column_1based",
            "source_partition",
        ],
    )

    sensitivity_rows = [
        {
            "analysis": "R1",
            "matrix": "UNCHANGED_DUAL_DOMAIN",
            "taxa": 40,
            "source_columns": 506,
            "retained_columns": 506,
            "removed_columns": 0,
            "retention_fraction": "1.000000",
            "minimum_standard_residues": "NOT_APPLICABLE",
            "sequences_removed": 0,
            "comparability": "SAME_MATRIX_AS_M09C_PRIMARY_DIFFERENT_MODEL",
        },
        {
            "analysis": "R2",
            "matrix": "DUAL_DOMAIN_OCCUPANCY20_PARTITIONWISE",
            "taxa": 40,
            "source_columns": 506,
            "retained_columns": masked_dual["retained_columns"],
            "removed_columns": 506 - int(masked_dual["retained_columns"]),
            "retention_fraction": f"{int(masked_dual['retained_columns']) / 506:.6f}",
            "minimum_standard_residues": 20,
            "sequences_removed": 0,
            "comparability": "NOT_COMPARABLE_DIFFERENT_DATA",
        },
        {
            "analysis": "R3",
            "matrix": "FULL_LENGTH_OCCUPANCY20",
            "taxa": 40,
            "source_columns": 1611,
            "retained_columns": masked_full["retained_columns"],
            "removed_columns": 1611 - int(masked_full["retained_columns"]),
            "retention_fraction": f"{int(masked_full['retained_columns']) / 1611:.6f}",
            "minimum_standard_residues": 20,
            "sequences_removed": 0,
            "comparability": "NOT_COMPARABLE_DIFFERENT_DATA",
        },
        {
            "analysis": "R4",
            "matrix": "UNCHANGED_DUAL_DOMAIN_DAYHOFF6_RECODED",
            "taxa": 40,
            "source_columns": 506,
            "retained_columns": 506,
            "removed_columns": 0,
            "retention_fraction": "1.000000",
            "minimum_standard_residues": "NOT_APPLICABLE",
            "sequences_removed": 0,
            "comparability": "NOT_COMPARABLE_DIFFERENT_STATE_SPACE",
        },
    ]
    write_tsv(
        staging / "M09D_ALIGNMENT_SENSITIVITY.tsv",
        sensitivity_rows,
        list(sensitivity_rows[0]),
    )
    return {
        "paths": {key: str(value) for key, value in paths.items()},
        "r2_retained_columns": masked_dual["retained_columns"],
        "r2_partitions": partition_rows,
        "r3_retained_columns": masked_full["retained_columns"],
        "r2_column_map_rows": len(masked_dual["retained_source_columns_1based"]),
        "r3_column_map_rows": len(masked_full["retained_source_columns_1based"]),
    }


def _prepare_composition_diagnostics(
    staging: Path, copied: Mapping[str, Path]
) -> dict[str, object]:
    dual = core.parse_fasta(copied["dual"], aligned=True)
    full = core.parse_fasta(copied["full"], aligned=True)
    primary_comp = {
        row["taxon"]: row for row in read_tsv(copied["primary_composition"])
    }
    full_comp = {row["taxon"]: row for row in read_tsv(copied["full_composition"])}
    domain_rows = read_tsv(copied["domain_extraction"])
    concerns: dict[str, list[str]] = {taxon: [] for taxon in dual}
    for row in domain_rows:
        if row["mandatory_completeness_concern"] == "YES":
            concerns[row["tree_label"]].append(
                f"{row['domain']}_coverage={row['profile_coverage']}_segments={row['segment_count']}"
            )
    if sum(len(values) for values in concerns.values()) != 5:
        raise ValueError("domain completeness concern count is not exactly five")
    primary_branch = {
        row["taxon"]: row
        for row in core.terminal_branch_diagnostics(
            _tree_terminal_lengths(copied["primary_tree"])
        )
    }
    full_branch = {
        row["taxon"]: row
        for row in core.terminal_branch_diagnostics(
            _tree_terminal_lengths(copied["full_tree"])
        )
    }
    sequence_map = {
        row["tree_label"]: row for row in read_tsv(copied["sequence_map"])
    }

    rows: list[dict[str, object]] = []
    for alignment_name, records, composition, branches in (
        ("DUAL_DOMAIN", dual, primary_comp, primary_branch),
        ("FULL_LENGTH", full, full_comp, full_branch),
    ):
        for row in core.composition_diagnostics(records):
            taxon = str(row["taxon"])
            original = composition[taxon]
            branch = branches[taxon]
            output = {
                "alignment": alignment_name,
                "taxon": taxon,
                "species": sequence_map[taxon]["species"],
                "alignment_columns": len(records[taxon]),
                "observed_standard_residues": row["observed_standard_residues"],
                "gap_or_ambiguity_fraction": f"{float(row['gap_or_ambiguity_fraction']):.8f}",
                "greater_than_50_percent_gap_or_ambiguity": (
                    "YES" if float(row["gap_or_ambiguity_fraction"]) > 0.5 else "NO"
                ),
                "jensen_shannon_distance": (
                    "NOT_COMPUTED"
                    if row["jensen_shannon_distance"] is None
                    else f"{float(row['jensen_shannon_distance']):.10f}"
                ),
                "original_iqtree_composition_status": original["status"].upper(),
                "original_iqtree_composition_p_value": original["p_value"],
                "terminal_branch_length": f"{float(branch['terminal_branch_length']):.10g}",
                "long_branch_threshold": f"{float(branch['long_branch_threshold']):.10g}",
                "long_branch_threshold_method": branch["threshold_method"],
                "long_branch_flag": "YES" if branch["long_branch_flag"] else "NO",
                "domain_completeness_concern": (
                    ";".join(concerns[taxon]) if concerns[taxon] else "NONE"
                ),
            }
            for amino_acid in core.STANDARD_AA:
                output[f"frequency_{amino_acid}"] = f"{float(row[f'frequency_{amino_acid}']):.10f}"
            rows.append(output)
    fields = list(rows[0])
    write_tsv(staging / "M09D_COMPOSITION_DIAGNOSTICS.tsv", rows, fields)
    primary_fail = [
        row["taxon"]
        for row in rows
        if row["alignment"] == "DUAL_DOMAIN"
        and row["original_iqtree_composition_status"] == "FAILED"
    ]
    full_fail = [
        row["taxon"]
        for row in rows
        if row["alignment"] == "FULL_LENGTH"
        and row["original_iqtree_composition_status"] == "FAILED"
    ]
    high_missing = [
        f"{row['alignment']}:{row['taxon']}"
        for row in rows
        if row["greater_than_50_percent_gap_or_ambiguity"] == "YES"
    ]
    return {
        "rows": len(rows),
        "primary_composition_failures": primary_fail,
        "full_length_composition_failures": full_fail,
        "greater_than_50_percent_gap_or_ambiguity_rows": len(high_missing),
        "greater_than_50_percent_gap_or_ambiguity_labels": high_missing,
        "domain_completeness_concerns": {
            taxon: values for taxon, values in concerns.items() if values
        },
    }


PERFORMANCE_FIELDS = [
    "step",
    "command",
    "exit_code",
    "wall_seconds",
    "max_rss_kb",
    "completed_utc",
]


def _parse_time_v(path: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    if not path.is_file():
        return result
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("\tMaximum resident set size"):
            result["max_rss_kb"] = int(line.rsplit(":", 1)[1].strip())
    return result


def _run_timed(
    command: Sequence[str],
    *,
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
    time_path: Path,
) -> dict[str, object]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.run(
            [str(TIME_EXECUTABLE), "-v", "-o", str(time_path), *command],
            cwd=cwd,
            stdout=stdout,
            stderr=stderr,
            check=False,
        )
    timing = _parse_time_v(time_path)
    return {
        "step": "",
        "command": " ".join(str(item) for item in command),
        "exit_code": process.returncode,
        "wall_seconds": f"{time.monotonic() - started:.3f}",
        "max_rss_kb": timing.get("max_rss_kb", ""),
        "completed_utc": utc_now(),
    }


def _append_performance(staging: Path, row: Mapping[str, object]) -> None:
    path = staging / "PERFORMANCE_LOG.tsv"
    rows = read_tsv(path) if path.is_file() else []
    rows.append(dict(row))
    write_tsv(path, rows, PERFORMANCE_FIELDS)


def _run_symmetry_tests(
    staging: Path, copied: Mapping[str, Path]
) -> dict[str, object]:
    output = staging / "PHYLOGENY/SYMMETRY"
    output.mkdir(parents=True, exist_ok=True)
    jobs = [
        (
            "DUAL_DOMAIN",
            copied["dual"],
            copied["partitions"],
            20261108,
        ),
        ("FULL_LENGTH", copied["full"], None, 20261110),
    ]
    combined: list[dict[str, object]] = []
    for name, alignment, partitions, seed in jobs:
        prefix = output / f"M09D_{name}_SYMTEST"
        command = [
            str(IQTREE3),
            "-s",
            str(alignment),
        ]
        if partitions is not None:
            command.extend(["-p", str(partitions)])
        command.extend(
            [
                "--seqtype",
                "AA",
                "--symtest-only",
                "-T",
                "2",
                "--seed",
                str(seed),
                "--prefix",
                str(prefix),
            ]
        )
        atomic_write_json(output / f"{name}.COMMAND.json", {"command": command})
        csv_path = prefix.with_suffix(".symtest.csv")
        if not csv_path.is_file():
            performance = _run_timed(
                command,
                cwd=staging,
                stdout_path=output / f"{name}.stdout.txt",
                stderr_path=output / f"{name}.stderr.txt",
                time_path=output / f"{name}.time.txt",
            )
            performance["step"] = f"symtest_{name.lower()}"
            _append_performance(staging, performance)
            if performance["exit_code"] != 0 or not csv_path.is_file():
                raise RuntimeError(f"IQ-TREE symmetry test failed for {name}")
        raw_rows: list[dict[str, str]] = []
        with csv_path.open(encoding="utf-8", newline="") as handle:
            filtered = (line for line in handle if not line.startswith("#"))
            raw_rows = list(csv.DictReader(filtered))
        for family, p_field in (
            ("SYMMETRY", "SymPval"),
            ("MARGINAL", "MarPval"),
            ("INTERNAL", "IntPval"),
        ):
            adjusted = core.benjamini_hochberg(
                [parse_optional_probability(row[p_field]) for row in raw_rows]
            )
            for row, q_value in zip(raw_rows, adjusted):
                combined.append(
                    {
                        "alignment": name,
                        "partition": row["Name"],
                        "test_family": family,
                        "raw_p_value": row[p_field],
                        "bh_q_value_within_alignment_family": (
                            "NA" if q_value is None else f"{q_value:.10g}"
                        ),
                        "significant_pair_count": row[
                            {"SYMMETRY": "SymSig", "MARGINAL": "MarSig", "INTERNAL": "IntSig"}[family]
                        ],
                        "non_significant_pair_count": row[
                            {"SYMMETRY": "SymNon", "MARGINAL": "MarNon", "INTERNAL": "IntNon"}[family]
                        ],
                        "q_lt_0.05": (
                            "NOT_COMPUTED" if q_value is None else "YES" if q_value < 0.05 else "NO"
                        ),
                        "taxon_pair_p_values": "NOT_EXPORTED_BY_IQTREE_3.1.2",
                    }
                )
    write_tsv(
        staging / "PHYLOGENY/SYMMETRY/M09D_SYMMETRY_TESTS.tsv",
        combined,
        list(combined[0]),
    )
    return {
        "rows": len(combined),
        "significant_partition_family_tests_q_lt_0.05": sum(
            row["q_lt_0.05"] == "YES" for row in combined
        ),
        "pair_specific_p_values": "NOT_AVAILABLE_FROM_IQTREE_OUTPUT",
    }


def _run_r4_capability_gate(staging: Path) -> dict[str, object]:
    output = staging / "PHYLOGENY/CAPABILITY/R4_DAYHOFF6"
    output.mkdir(parents=True, exist_ok=True)
    synthetic = output / "MORPH6_SYNTHETIC.faa"
    records = {
        f"T{index + 1}": ("012345"[index:] + "012345"[:index]) * 10
        for index in range(6)
    }
    write_fasta(synthetic, records)
    prefix = output / "MORPH6_GTR_FO_G4_CAPABILITY"
    command = [
        str(IQTREE3),
        "-s",
        str(synthetic),
        "--seqtype",
        "MORPH{6}",
        "-m",
        "GTR+FO+G4",
        "-n",
        "0",
        "-T",
        "1",
        "--seed",
        "20261111",
        "--prefix",
        str(prefix),
    ]
    atomic_write_json(output / "COMMAND.json", {"command": command})
    performance = _run_timed(
        command,
        cwd=staging,
        stdout_path=output / "stdout.txt",
        stderr_path=output / "stderr.txt",
        time_path=output / "time.txt",
    )
    performance["step"] = "r4_dayhoff6_capability_gate"
    _append_performance(staging, performance)
    passed = performance["exit_code"] == 0 and prefix.with_suffix(".iqtree").is_file()
    result = {
        "status": "PASS_SUPPORTED" if passed else "NOT_COMPUTED_UNSUPPORTED",
        "supported": passed,
        "command": command,
        "exit_code": performance["exit_code"],
        "synthetic_only": True,
        "model": "GTR+FO+G4",
        "seqtype": "MORPH{6}",
    }
    atomic_write_json(output / "CAPABILITY_RESULT.json", result)
    return result


def prepare_checkpoint_includes(relative: Path) -> bool:
    """Return whether a prepare artifact is immutable across later run stages."""

    if relative.as_posix() in {"RUN_STATE.json", "PERFORMANCE_LOG.tsv"}:
        return False
    if relative.parts and relative.parts[0] == "CHECKPOINTS":
        return False
    return True


def prepare_workflow(project_root: Path, staging: Path | None = None) -> Path:
    root = project_root.resolve()
    if staging is None:
        staging, _final = _new_staging(root)
    else:
        staging = safe_staging_path(root / RELEASE_PARENT_RELATIVE, staging)
        if not staging.is_dir():
            raise FileNotFoundError(f"resume staging does not exist: {staging}")
    checkpoint = staging / "CHECKPOINTS/PREPARE_COMPLETE.json"
    if checkpoint.is_file():
        record = json.loads(checkpoint.read_text(encoding="utf-8"))
        for relative, expected in record["outputs"].items():
            path = staging / relative
            if not path.is_file() or sha256_file(path) != expected:
                raise ValueError(f"prepare checkpoint output mismatch: {relative}")
        return staging

    audit = source_integrity_audit(root)
    atomic_write_json(staging / "SOURCE_INTEGRITY_AUDIT.json", audit)
    copied = _copy_locked_inputs(root, staging)
    _write_source_lock(root, staging, copied)
    transforms = _prepare_transforms(staging, copied)
    composition = _prepare_composition_diagnostics(staging, copied)
    symmetry = _run_symmetry_tests(staging, copied)
    r4_capability = _run_r4_capability_gate(staging)
    hypotheses = {
        "registered_before_biological_tree_execution": True,
        "registered_utc": utc_now(),
        "questions": [
            "whether composition warnings co-occur with missingness, incomplete domains or long branches",
            "whether exact-pair and component splits persist under LG+C20+F+G4",
            "whether they persist after taxon-blind 20-of-40 occupancy masking",
            "whether they persist under literal Dayhoff-6 recoding",
        ],
        "strong_support_rule": "SH-aLRT >=80 and UFBoot >=95",
        "taxon_removal": "NONE",
        "rooting": "UNROOTED",
        "maximum_threads": 2,
        "classification_rules": "docs/superpowers/specs/2026-11-07-m09d-composition-robustness-design.md",
        "pair_specific_symmetry_limitation": "IQ-TREE_3.1.2_EXPORTS_PARTITION_MAXIMUM_SUMMARIES_NOT_INDIVIDUAL_PAIR_P_VALUES",
    }
    atomic_write_json(staging / "M09D_PRE_REGISTERED_HYPOTHESES.json", hypotheses)
    summary = {
        "status": "PASS_PREPARED_FOR_ROBUSTNESS_TREES",
        "completed_utc": utc_now(),
        "source_audit": audit,
        "transforms": transforms,
        "composition": composition,
        "symmetry": symmetry,
        "r4_capability": r4_capability,
    }
    atomic_write_json(staging / "M09D_PREPARE_SUMMARY.json", summary)
    outputs = {}
    for path in sorted(staging.rglob("*")):
        relative = path.relative_to(staging)
        if path.is_file() and prepare_checkpoint_includes(relative):
            outputs[str(relative)] = sha256_file(path)
    atomic_write_json(
        checkpoint,
        {
            "schema": "M09D_PREPARE_CHECKPOINT_V1",
            "created_utc": utc_now(),
            "outputs": outputs,
        },
    )
    state = json.loads((staging / "RUN_STATE.json").read_text(encoding="utf-8"))
    state["status"] = "PREPARE_COMPLETE"
    state["updated_utc"] = utc_now()
    atomic_write_json(staging / "RUN_STATE.json", state)
    print(f"M09D_STAGING={staging}", flush=True)
    print(f"M09D_PREPARE=PASS", flush=True)
    print(f"M09D_R4_CAPABILITY={r4_capability['status']}", flush=True)
    return staging


def parse_iqtree_metrics(report: str, *, declared_model: str) -> dict[str, object]:
    patterns = {
        "log_likelihood": r"Log-likelihood of the tree:\s*([-+0-9.eE]+)",
        "free_parameters": r"Number of free parameters \([^\n]+\):\s*([0-9]+)",
        "aic": r"Akaike information criterion \(AIC\) score:\s*([-+0-9.eE]+)",
        "aicc": r"Corrected Akaike information criterion \(AICc\) score:\s*([-+0-9.eE]+)",
        "bic": r"Bayesian information criterion \(BIC\) score:\s*([-+0-9.eE]+)",
    }
    selected = re.search(r"Best-fit model according to BIC:\s*(.+)", report)
    values: dict[str, object] = {
        "declared_model": declared_model,
        "selected_model": selected.group(1).strip() if selected else declared_model,
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, report)
        if match is None:
            values[key] = "NOT_REPORTED"
        elif key == "free_parameters":
            values[key] = int(match.group(1))
        else:
            values[key] = float(match.group(1))
    return values


def iqtree_artifacts_complete(prefix: Path, expected_taxa: set[str]) -> bool:
    required = [
        prefix.with_suffix(".treefile"),
        prefix.with_suffix(".contree"),
        prefix.with_suffix(".iqtree"),
        prefix.with_suffix(".log"),
        prefix.with_suffix(".ufboot"),
        prefix.with_suffix(".runtrees"),
    ]
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        return False
    try:
        if _tree_leaf_set(required[0]) != expected_taxa:
            return False
        if _tree_leaf_set(required[1]) != expected_taxa:
            return False
    except (ValueError, OSError):
        return False
    if sum(bool(line.strip()) for line in required[4].read_text(encoding="utf-8").splitlines()) != 1000:
        return False
    if sum(bool(line.strip()) for line in required[5].read_text(encoding="utf-8").splitlines()) != 5:
        return False
    if "Analysis results written to:" not in required[3].read_text(
        encoding="utf-8", errors="replace"
    ):
        return False
    return True


def _analysis_configuration(staging: Path, analysis_id: str) -> dict[str, object]:
    base = staging / "PHYLOGENY/ALIGNMENTS"
    names = {
        "R1": "R1_SITE_HETEROGENEOUS",
        "R2": "R2_DUAL_OCCUPANCY20",
        "R3": "R3_FULL_OCCUPANCY20",
        "R4": "R4_DAYHOFF6",
    }
    alignments = {
        "R1": base / "M09D_R1_DUAL_DOMAIN_UNCHANGED.faa",
        "R2": base / "M09D_R2_DUAL_DOMAIN_OCCUPANCY20.faa",
        "R3": base / "M09D_R3_FULL_LENGTH_OCCUPANCY20.faa",
        "R4": base / "M09D_R4_DUAL_DOMAIN_DAYHOFF6.faa",
    }
    if analysis_id not in names:
        raise ValueError(f"unknown analysis: {analysis_id}")
    directory = staging / "PHYLOGENY/TREES" / names[analysis_id]
    prefix = directory / f"M09D_{analysis_id}"
    partitions = (
        base / "M09D_R2_MASKED_PARTITIONS.nex" if analysis_id == "R2" else None
    )
    return {
        "name": names[analysis_id],
        "directory": directory,
        "alignment": alignments[analysis_id],
        "partitions": partitions,
        "prefix": prefix,
    }


def _verify_prepare_checkpoint(staging: Path) -> None:
    checkpoint = staging / "CHECKPOINTS/PREPARE_COMPLETE.json"
    if not checkpoint.is_file():
        raise ValueError("M09D prepare checkpoint is absent")
    record = json.loads(checkpoint.read_text(encoding="utf-8"))
    for relative, expected in record["outputs"].items():
        path = staging / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise ValueError(f"prepare checkpoint output mismatch: {relative}")


def _verify_live_code_lock(staging: Path) -> None:
    rows = read_tsv(staging / "M09D_INITIAL_SOURCE_LOCK.tsv")
    expected = {
        row["role"]: row["sha256"]
        for row in rows
        if row["role"] in {"M09D_CORE_CODE", "M09D_WORKFLOW_CODE"}
    }
    observed = {
        "M09D_CORE_CODE": sha256_file(TOOL_DIR / "m09d_core.py"),
        "M09D_WORKFLOW_CODE": sha256_file(TOOL_DIR / "m09d_workflow.py"),
    }
    if expected != observed:
        raise ValueError(
            "live M09D implementation differs from the staging initial source lock"
        )


def run_analysis(staging: Path, analysis_id: str) -> dict[str, object]:
    staging = staging.resolve()
    _verify_prepare_checkpoint(staging)
    _verify_live_code_lock(staging)
    configuration = _analysis_configuration(staging, analysis_id)
    directory = Path(configuration["directory"])
    directory.mkdir(parents=True, exist_ok=True)
    alignment = Path(configuration["alignment"])
    partitions_value = configuration["partitions"]
    partitions = Path(partitions_value) if partitions_value is not None else None
    prefix = Path(configuration["prefix"])

    if analysis_id == "R4":
        capability = json.loads(
            (
                staging
                / "PHYLOGENY/CAPABILITY/R4_DAYHOFF6/CAPABILITY_RESULT.json"
            ).read_text(encoding="utf-8")
        )
        if not capability.get("supported"):
            result = {
                "analysis": analysis_id,
                "status": "NOT_COMPUTED_UNSUPPORTED",
                "capability": capability,
            }
            atomic_write_json(directory / "METRICS.json", result)
            return result

    command = build_iqtree_command(
        IQTREE3,
        analysis_id,
        alignment,
        prefix,
        partitions=partitions,
    )
    inputs = [alignment] + ([partitions] if partitions is not None else [])
    provenance = build_provenance(
        command,
        inputs,
        [IQTREE3, TIME_EXECUTABLE],
        [TOOL_DIR / "m09d_core.py", TOOL_DIR / "m09d_workflow.py"],
    )
    provenance_path = directory / "PROVENANCE.json"
    if provenance_path.is_file():
        existing = json.loads(provenance_path.read_text(encoding="utf-8"))
        if not provenance_matches(existing, provenance):
            raise ValueError(f"{analysis_id} resume provenance does not match")
    else:
        output_candidates = list(directory.glob(f"{prefix.name}.*"))
        if output_candidates:
            raise ValueError(f"{analysis_id} has outputs without provenance")
        atomic_write_json(provenance_path, provenance)
        atomic_write_json(
            directory / "COMMAND.json",
            {
                "analysis": analysis_id,
                "command": command,
                "threads": 2,
                "runs": 5,
                "ufboot_replicates": 1000,
                "shalrt_replicates": 1000,
                "rooting": "UNROOTED",
            },
        )
    expected_taxa = set(
        core.parse_fasta(
            staging / "PHYLOGENY/INPUTS/M09D_M09C_DUAL_DOMAIN_40_UNCHANGED.faa",
            aligned=True,
        )
    )
    if iqtree_artifacts_complete(prefix, expected_taxa):
        report = prefix.with_suffix(".iqtree").read_text(
            encoding="utf-8", errors="replace"
        )
        metrics = parse_iqtree_metrics(
            report, declared_model=str(ANALYSIS_LOCK[analysis_id]["model"])
        )
        metrics.update(
            {
                "analysis": analysis_id,
                "status": "PASS_COMPLETE_REUSED",
                "taxa": 40,
                "ufboot_trees": 1000,
                "independent_runs": 5,
                "treefile_sha256": sha256_file(prefix.with_suffix(".treefile")),
            }
        )
        atomic_write_json(directory / "METRICS.json", metrics)
        return metrics

    attempt = len(list(directory.glob("ATTEMPT_*.time.txt"))) + 1
    attempt_tag = f"ATTEMPT_{attempt:03d}"
    performance = _run_timed(
        command,
        cwd=directory,
        stdout_path=directory / f"{attempt_tag}.stdout.txt",
        stderr_path=directory / f"{attempt_tag}.stderr.txt",
        time_path=directory / f"{attempt_tag}.time.txt",
    )
    performance["step"] = f"iqtree_{analysis_id.lower()}_{attempt_tag.lower()}"
    _append_performance(staging, performance)
    if performance["exit_code"] != 0:
        raise RuntimeError(
            f"IQ-TREE {analysis_id} failed with exit code {performance['exit_code']}"
        )
    if not iqtree_artifacts_complete(prefix, expected_taxa):
        raise RuntimeError(f"IQ-TREE {analysis_id} returned success but outputs are incomplete")
    report = prefix.with_suffix(".iqtree").read_text(encoding="utf-8", errors="replace")
    metrics = parse_iqtree_metrics(
        report, declared_model=str(ANALYSIS_LOCK[analysis_id]["model"])
    )
    metrics.update(
        {
            "analysis": analysis_id,
            "status": "PASS_COMPLETE",
            "taxa": 40,
            "alignment_columns": alignment_column_count(alignment),
            "ufboot_trees": 1000,
            "independent_runs": 5,
            "treefile_sha256": sha256_file(prefix.with_suffix(".treefile")),
            "completed_utc": utc_now(),
        }
    )
    atomic_write_json(directory / "METRICS.json", metrics)
    atomic_write_json(
        staging / f"CHECKPOINTS/{analysis_id}_COMPLETE.json",
        {
            "schema": "M09D_IQTREE_CHECKPOINT_V1",
            "analysis": analysis_id,
            "provenance_sha256": sha256_file(provenance_path),
            "outputs": {
                path.name: sha256_file(path)
                for path in sorted(directory.glob(f"{prefix.name}.*"))
                if path.is_file()
            },
        },
    )
    print(
        f"M09D_{analysis_id}=PASS model={metrics['selected_model']} "
        f"tree_sha256={metrics['treefile_sha256']}",
        flush=True,
    )
    return metrics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    parser.add_argument("command", choices=["prepare", "run"])
    parser.add_argument("--staging", type=Path)
    parser.add_argument("--analysis", choices=sorted(ANALYSIS_LOCK))
    arguments = parser.parse_args(argv)
    if arguments.command == "prepare":
        prepare_workflow(arguments.project_root, arguments.staging)
        return 0
    if arguments.command == "run":
        if arguments.staging is None or arguments.analysis is None:
            parser.error("run requires --staging and --analysis")
        run_analysis(arguments.staging, arguments.analysis)
        return 0
    raise AssertionError("unreachable command")


if __name__ == "__main__":
    raise SystemExit(main())
