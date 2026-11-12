"""M09D scientific re-adjudication and immutable release finalization."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import m09d_core as core
import m09d_figures as figures
import m09d_workflow as workflow


NEW_ANALYSES = ("R1", "R2", "R3", "R4")
ALL_TREES = ("M09C_PRIMARY", "M09C_FULL", *NEW_ANALYSES)


def _fmt_number(value: float) -> str:
    return f"{value:g}"


def support_for(
    group: Iterable[str],
    universe: Iterable[str],
    supports: Mapping[tuple[str, ...], tuple[float, float]],
) -> str:
    key = core.canonical_side(group, universe)
    support = supports.get(key)
    if support is None:
        return "NOT_PRESENT"
    return f"{_fmt_number(support[0])}/{_fmt_number(support[1])}"


def rollup_relationships(classes: Iterable[str]) -> str:
    values = list(classes)
    if not values:
        raise ValueError("cannot roll up no relationship classifications")
    priority = (
        "CONTRADICTED",
        "UNRESOLVED",
        "POSSIBLE_CO_ORTHOLOGY",
        "SUPPORTED_WITH_MODEL_CAVEAT",
        "ROBUST_COMPUTATIONAL_CANDIDATE",
    )
    for value in priority:
        if value in values:
            return value
    raise ValueError(f"unknown relationship classification(s): {values}")


def _tree_paths(staging: Path) -> dict[str, Path]:
    return {
        "M09C_PRIMARY": staging / "PHYLOGENY/BASELINES/M09C_PRIMARY.treefile",
        "M09C_FULL": staging / "PHYLOGENY/BASELINES/M09C_FULL_LENGTH.treefile",
        "R1": staging
        / "PHYLOGENY/TREES/R1_SITE_HETEROGENEOUS/M09D_R1.treefile",
        "R2": staging / "PHYLOGENY/TREES/R2_DUAL_OCCUPANCY20/M09D_R2.treefile",
        "R3": staging / "PHYLOGENY/TREES/R3_FULL_OCCUPANCY20/M09D_R3.treefile",
        "R4": staging / "PHYLOGENY/TREES/R4_DAYHOFF6/M09D_R4.treefile",
    }


def _load_tree_supports(staging: Path) -> tuple[
    set[str],
    dict[str, dict[tuple[str, ...], tuple[float, float]]],
    dict[str, set[frozenset[str]]],
]:
    paths = _tree_paths(staging)
    universe: set[str] | None = None
    all_supports: dict[str, dict[tuple[str, ...], tuple[float, float]]] = {}
    strong: dict[str, set[frozenset[str]]] = {}
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"mandatory tree is absent: {name}: {path}")
        leaves = workflow._tree_leaf_set(path)
        if universe is None:
            universe = leaves
            core.validate_taxa(universe)
        elif leaves != universe:
            raise ValueError(f"tree leaf-set mismatch: {name}")
        parsed = core.parse_split_supports(path.read_text(encoding="utf-8"))
        all_supports[name] = parsed
        strong[name] = {
            frozenset(split)
            for split, support in parsed.items()
            if support[0] >= core.STRONG_SH_ALRT and support[1] >= core.STRONG_UFBOOT
        }
    assert universe is not None
    return universe, all_supports, strong


def _yes_no_strong(label: str) -> bool:
    return label != "NOT_PRESENT"


def _label_is_strong(label: str) -> bool:
    if label == "NOT_PRESENT":
        return False
    sh_alrt, ufboot = (float(value) for value in label.split("/"))
    return sh_alrt >= core.STRONG_SH_ALRT and ufboot >= core.STRONG_UFBOOT


def _relationship_row(
    *,
    maize_alias: str,
    maize_gene: str,
    sorghum: str,
    rbh: bool,
    component_text: str,
    universe: set[str],
    supports: Mapping[str, Mapping[tuple[str, ...], tuple[float, float]]],
    strong: Mapping[str, set[frozenset[str]]],
    primary_failures: set[str],
    full_failures: set[str],
    domain_evidence: str,
    collinearity_evidence: str,
    inherited_class: str,
) -> dict[str, object]:
    pair = frozenset({maize_alias, sorghum})
    component = frozenset(value for value in component_text.split(";") if value)
    if len(pair) != 2 or not pair <= universe:
        raise ValueError(f"invalid relationship pair: {maize_alias}, {sorghum}")
    if not pair <= component or not component <= universe:
        raise ValueError(f"invalid collinearity component for {maize_alias}, {sorghum}")

    exact_labels = {
        name: support_for(pair, universe, supports[name]) for name in ALL_TREES
    }
    component_labels = {
        name: support_for(component, universe, supports[name]) for name in ALL_TREES
    }
    new_exact = {name: _label_is_strong(exact_labels[name]) for name in NEW_ANALYSES}
    new_component = {
        name: _label_is_strong(component_labels[name]) for name in NEW_ANALYSES
    }
    original_exact = {
        name: _label_is_strong(exact_labels[name])
        for name in ("M09C_PRIMARY", "M09C_FULL")
    }
    original_component = {
        name: _label_is_strong(component_labels[name])
        for name in ("M09C_PRIMARY", "M09C_FULL")
    }
    alternatives = core.reproducible_exact_alternatives(
        pair, {name: strong[name] for name in NEW_ANALYSES}
    )
    caveats: list[str] = []
    if pair & primary_failures:
        caveats.append("PAIR_MEMBER_FAILED_ORIGINAL_PRIMARY_COMPOSITION_TEST")
    if original_exact["M09C_PRIMARY"] != new_exact["R2"]:
        caveats.append("EXACT_SUPPORT_CHANGED_AFTER_DUAL_OCCUPANCY_MASK")
    if original_exact["M09C_FULL"] != new_exact["R3"]:
        caveats.append("EXACT_SUPPORT_CHANGED_AFTER_FULL_LENGTH_OCCUPANCY_MASK")
    if pair & full_failures and original_exact["M09C_FULL"] != new_exact["R3"]:
        caveats.append("FULL_LENGTH_COMPOSITION_FAILURE_WITH_NONCONCORDANT_MASKED_RESULT")
    classification = core.adjudicate_m09d_relationship(
        rbh=rbh,
        original_exact=original_exact,
        new_exact=new_exact,
        original_component=original_component,
        new_component=new_component,
        reproducible_alternative_count=len(alternatives),
        direct_caveat=bool(caveats),
    )
    limitations = [
        "COMPUTATIONAL_RELATIONSHIP_NOT_EXPERIMENTAL_OR_FUNCTIONAL_PROOF",
        "UNROOTED_TWO_SPECIES_FAMILY_TREE",
        "PAIR_SPECIFIC_SYMMETRY_P_VALUES_NOT_EXPORTED_BY_IQTREE",
    ]
    if not rbh:
        limitations.append("NO_RECIPROCAL_RANK1_PREVENTS_SUPPORTED_CANDIDATE_CLASS")
    if caveats:
        limitations.extend(caveats)
    return {
        "maize_family_alias": maize_alias,
        "maize_gene_id": maize_gene,
        "sorghum_candidate_id": sorghum,
        "reciprocal_best_hit_status": "YES" if rbh else "NO",
        "domain_evidence": domain_evidence,
        "collinearity_evidence": collinearity_evidence,
        "collinearity_component": component_text,
        "m09c_inherited_classification": inherited_class,
        **{f"{name.lower()}_exact_support_shalrt_ufboot": exact_labels[name] for name in ALL_TREES},
        **{
            f"{name.lower()}_component_support_shalrt_ufboot": component_labels[name]
            for name in ALL_TREES
        },
        "r1_exact_strong": "YES" if new_exact["R1"] else "NO",
        "r2_exact_strong": "YES" if new_exact["R2"] else "NO",
        "r3_exact_strong": "YES" if new_exact["R3"] else "NO",
        "r4_exact_strong": "YES" if new_exact["R4"] else "NO",
        "new_component_strong_tree_count": sum(new_component.values()),
        "reproducible_strong_exact_alternatives": (
            ";".join(
                f"{','.join(alternative)}@{','.join(tree_names)}"
                for alternative, tree_names in alternatives.items()
            )
            if alternatives
            else "NONE"
        ),
        "composition_and_mask_caveats": ";".join(caveats) if caveats else "NONE",
        "final_computational_classification": classification,
        "scientific_limitations": ";".join(dict.fromkeys(limitations)),
    }


def adjudicate_relationships(staging: Path) -> dict[str, object]:
    universe, supports, strong = _load_tree_supports(staging)
    evidence = staging / "PHYLOGENY/SOURCE_EVIDENCE"
    diagnostic_rows = workflow.read_tsv(staging / "M09D_COMPOSITION_DIAGNOSTICS.tsv")
    primary_failures = {
        row["taxon"]
        for row in diagnostic_rows
        if row["alignment"] == "DUAL_DOMAIN"
        and row["original_iqtree_composition_status"] == "FAILED"
    }
    full_failures = {
        row["taxon"]
        for row in diagnostic_rows
        if row["alignment"] == "FULL_LENGTH"
        and row["original_iqtree_composition_status"] == "FAILED"
    }
    candidates = workflow.read_tsv(
        evidence / "M09C_CANDIDATE_PAIR_ADJUDICATION.tsv"
    )
    no_rbh = workflow.read_tsv(evidence / "M09C_NO_RBH_SEVEN_CASES.tsv")
    relationship_rows: list[dict[str, object]] = []
    for row in candidates:
        relationship_rows.append(
            _relationship_row(
                maize_alias=row["maize_family_alias"],
                maize_gene=row["maize_gene_id"],
                sorghum=row["sorghum_candidate_id"],
                rbh=True,
                component_text=row["collinearity_component"],
                universe=universe,
                supports=supports,
                strong=strong,
                primary_failures=primary_failures,
                full_failures=full_failures,
                domain_evidence=row["domain_evidence"],
                collinearity_evidence=row["collinearity_evidence"],
                inherited_class=row["final_computational_classification"],
            )
        )
    for row in no_rbh:
        for sorghum in row["sorghum_collinear_family_partners"].split(";"):
            relationship_rows.append(
                _relationship_row(
                    maize_alias=row["maize_family_alias"],
                    maize_gene=row["maize_gene_id"],
                    sorghum=sorghum,
                    rbh=False,
                    component_text=row["collinearity_component"],
                    universe=universe,
                    supports=supports,
                    strong=strong,
                    primary_failures=primary_failures,
                    full_failures=full_failures,
                    domain_evidence=row["domain_evidence"],
                    collinearity_evidence=row["collinearity_evidence"],
                    inherited_class=row["final_computational_classification"],
                )
            )
    if len(relationship_rows) != 24:
        raise ValueError(f"expected 24 pair hypotheses, observed {len(relationship_rows)}")
    relationship_rows.sort(
        key=lambda row: (str(row["maize_family_alias"]), str(row["sorghum_candidate_id"]))
    )
    workflow.write_tsv(
        staging / "M09D_ORTHOLOGY_EVIDENCE_MATRIX.tsv",
        relationship_rows,
        list(relationship_rows[0]),
    )

    by_maize: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in relationship_rows:
        by_maize[str(row["maize_family_alias"])].append(row)
    all23_source = workflow.read_tsv(evidence / "M09C_ALL_23_ORTHOLOGY_REVIEW.tsv")
    all23: list[dict[str, object]] = []
    for source in all23_source:
        alias = source["maize_family_alias"]
        relationships = by_maize[alias]
        classes = [str(row["final_computational_classification"]) for row in relationships]
        final_class = rollup_relationships(classes)
        all23.append(
            {
                "maize_family_alias": alias,
                "maize_gene_id": source["maize_gene_id"],
                "sorghum_candidate_ids": ";".join(
                    str(row["sorghum_candidate_id"]) for row in relationships
                ),
                "relationship_count": len(relationships),
                "reciprocal_best_hit_status": source["reciprocal_best_hit_status"],
                "m09c_final_classification": source["final_computational_classification"],
                "m09d_pair_classifications": ";".join(
                    f"{row['sorghum_candidate_id']}:{row['final_computational_classification']}"
                    for row in relationships
                ),
                "composition_and_mask_caveats": ";".join(
                    sorted(
                        {
                            str(row["composition_and_mask_caveats"])
                            for row in relationships
                            if row["composition_and_mask_caveats"] != "NONE"
                        }
                    )
                )
                or "NONE",
                "final_computational_classification": final_class,
                "scientific_limitations": ";".join(
                    sorted(
                        {
                            str(row["scientific_limitations"])
                            for row in relationships
                        }
                    )
                ),
            }
        )
    all23.sort(key=lambda row: str(row["maize_family_alias"]))
    if len(all23) != 23 or set(by_maize) != {f"ZmNRG2_{index:03d}" for index in range(1, 24)}:
        raise ValueError("M09D final review does not account for all 23 maize genes")
    workflow.write_tsv(
        staging / "M09D_ALL_23_FINAL_REVIEW.tsv", all23, list(all23[0])
    )
    seven = [row for row in all23 if row["reciprocal_best_hit_status"] == "NO"]
    if len(seven) != 7:
        raise ValueError(f"M09D NO_RBH rollup has {len(seven)} rows, expected 7")
    workflow.write_tsv(
        staging / "M09D_SEVEN_NO_RBH_CASES.tsv", seven, list(seven[0])
    )
    return {
        "relationship_rows": relationship_rows,
        "all23_rows": all23,
        "seven_rows": seven,
        "relationship_classification_counts": dict(
            sorted(Counter(row["final_computational_classification"] for row in relationship_rows).items())
        ),
        "all23_classification_counts": dict(
            sorted(Counter(row["final_computational_classification"] for row in all23).items())
        ),
        "focal_zmnrg2_001": next(
            row for row in all23 if row["maize_family_alias"] == "ZmNRG2_001"
        ),
        "primary_composition_failures": sorted(primary_failures),
        "full_length_composition_failures": sorted(full_failures),
    }


def _warning_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return sorted(
        {
            line.strip()
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip().startswith("WARNING:")
        }
    )


def build_model_comparison(project_root: Path, staging: Path) -> list[dict[str, object]]:
    m09c = project_root / workflow.M09C_RELATIVE
    baseline_dir = staging / "PHYLOGENY/BASELINES"
    baseline_sources = {
        "M09C_PRIMARY": (
            m09c / "PHYLOGENY/TREES/PRIMARY/M09C_PRIMARY.iqtree",
            m09c / "PHYLOGENY/TREES/PRIMARY/M09C_PRIMARY.log",
            "MFP+MERGE",
            "UNCHANGED_DUAL_DOMAIN_506_AA",
            "DUAL_506_AA",
        ),
        "M09C_FULL": (
            m09c / "PHYLOGENY/TREES/FULL_LENGTH/M09C_FULL_LENGTH.iqtree",
            m09c / "PHYLOGENY/TREES/FULL_LENGTH/M09C_FULL_LENGTH.log",
            "MFP",
            "UNCHANGED_FULL_LENGTH_1611_AA",
            "FULL_1611_AA",
        ),
    }
    rows: list[dict[str, object]] = []
    for name, (report_path, log_path, declared, matrix, group) in baseline_sources.items():
        shutil.copy2(report_path, baseline_dir / report_path.name)
        shutil.copy2(log_path, baseline_dir / log_path.name)
        metrics = workflow.parse_iqtree_metrics(
            report_path.read_text(encoding="utf-8", errors="replace"),
            declared_model=declared,
        )
        rows.append(
            {
                "analysis": name,
                "matrix": matrix,
                "taxa": 40,
                "columns": 506 if name == "M09C_PRIMARY" else 1611,
                "declared_model": metrics["declared_model"],
                "selected_model": metrics["selected_model"],
                "log_likelihood": metrics["log_likelihood"],
                "free_parameters": metrics["free_parameters"],
                "aic": metrics["aic"],
                "aicc": metrics["aicc"],
                "bic": metrics["bic"],
                "information_criterion_comparison_group": group,
                "information_criterion_interpretation": (
                    "COMPARABLE_ONLY_WITH_OTHER_ROWS_IN_SAME_GROUP"
                    if group == "DUAL_506_AA"
                    else "NO_NEW_SAME_MATRIX_COMPARATOR"
                ),
                "warnings": "; ".join(_warning_lines(log_path)) or "NONE",
            }
        )
    config = {
        "R1": ("R1_SITE_HETEROGENEOUS", "UNCHANGED_DUAL_DOMAIN_506_AA", "DUAL_506_AA"),
        "R2": ("R2_DUAL_OCCUPANCY20", "DUAL_DOMAIN_OCCUPANCY20_363_AA", "DIFFERENT_DATA_R2"),
        "R3": ("R3_FULL_OCCUPANCY20", "FULL_LENGTH_OCCUPANCY20_735_AA", "DIFFERENT_DATA_R3"),
        "R4": ("R4_DAYHOFF6", "UNCHANGED_DUAL_DOMAIN_506_DAYHOFF6", "DIFFERENT_STATE_SPACE_R4"),
    }
    for analysis, (directory, matrix, group) in config.items():
        tree_dir = staging / "PHYLOGENY/TREES" / directory
        metrics = json.loads((tree_dir / "METRICS.json").read_text(encoding="utf-8"))
        if not str(metrics.get("status", "")).startswith("PASS_COMPLETE"):
            raise ValueError(f"mandatory {analysis} metrics are not complete")
        prefix = tree_dir / f"M09D_{analysis}"
        rows.append(
            {
                "analysis": analysis,
                "matrix": matrix,
                "taxa": metrics["taxa"],
                "columns": metrics["alignment_columns"],
                "declared_model": metrics["declared_model"],
                "selected_model": metrics["selected_model"],
                "log_likelihood": metrics["log_likelihood"],
                "free_parameters": metrics["free_parameters"],
                "aic": metrics["aic"],
                "aicc": metrics["aicc"],
                "bic": metrics["bic"],
                "information_criterion_comparison_group": group,
                "information_criterion_interpretation": (
                    "COMPARABLE_ONLY_WITH_OTHER_ROWS_IN_SAME_GROUP"
                    if group == "DUAL_506_AA"
                    else "NOT_COMPARABLE_DIFFERENT_DATA_OR_STATE_SPACE"
                ),
                "warnings": "; ".join(_warning_lines(prefix.with_suffix(".log"))) or "NONE",
            }
        )
    workflow.write_tsv(staging / "M09D_MODEL_COMPARISON.tsv", rows, list(rows[0]))
    return rows


def generate_figures(
    staging: Path, adjudication: Mapping[str, object], model_rows: list[dict[str, object]]
) -> dict[str, object]:
    sequence_map = workflow.read_tsv(
        staging / "PHYLOGENY/SOURCE_EVIDENCE/M09C_SEQUENCE_MAP.tsv"
    )
    class_by_maize = {
        str(row["maize_family_alias"]): str(row["final_computational_classification"])
        for row in adjudication["all23_rows"]
    }
    source_rows = []
    for row in sequence_map:
        label = row["tree_label"]
        source_rows.append(
            {
                "tree_label": label,
                "species": row["species"],
                "display_label": label,
                "family_alias": row["family_alias"],
                "gene_id": row["gene_id"],
                "protein_id": row["protein_id"],
                "source_fasta_id": row["source_fasta_id"],
                "length_aa": row["length_aa"],
                "sequence_sha256": row["sequence_sha256"],
                "classification": class_by_maize.get(label, "CONTEXT_ONLY"),
            }
        )
    model_by_analysis = {str(row["analysis"]): row for row in model_rows}
    titles = {
        "R1": "M09D R1 — unchanged dual-domain alignment, LG+C20+F+G4",
        "R2": f"M09D R2 — dual-domain occupancy mask, {model_by_analysis['R2']['selected_model']}",
        "R3": f"M09D R3 — full-length occupancy mask, {model_by_analysis['R3']['selected_model']}",
        "R4": "M09D R4 — Dayhoff-6 dual-domain recoding, GTR+FO+G4",
    }
    figure_dir = staging / "PHYLOGENY/FIGURES"
    validations: list[dict[str, object]] = []
    trace_rows: list[dict[str, object]] = []
    tree_paths = _tree_paths(staging)
    for analysis in NEW_ANALYSES:
        prefix = figure_dir / f"M09D_{analysis}_UNROOTED_ROBUSTNESS_TREE"
        artifacts = figures.render_phylogram(
            tree_paths[analysis], source_rows, prefix, titles[analysis]
        )
        validation = figures.verify_figure_artifacts(
            tree_paths[analysis], source_rows, artifacts
        )
        validation["artifacts"] = [
            str(Path(path).relative_to(staging)) for path in validation["artifacts"]
        ]
        validation["analysis"] = analysis
        validation["tree_sha256"] = workflow.sha256_file(tree_paths[analysis])
        validations.append(validation)
        for row in source_rows:
            trace_rows.append(
                {
                    "analysis": analysis,
                    "tree_sha256": validation["tree_sha256"],
                    **row,
                }
            )

    relationship_by_maize: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in adjudication["relationship_rows"]:
        relationship_by_maize[str(row["maize_family_alias"])].append(row)

    def matrix_state(alias: str, analysis: str) -> str:
        rows = relationship_by_maize[alias]
        exact = [
            _label_is_strong(str(row[f"{analysis}_exact_support_shalrt_ufboot"]))
            for row in rows
        ]
        component = [
            _label_is_strong(str(row[f"{analysis}_component_support_shalrt_ufboot"]))
            for row in rows
        ]
        if all(exact):
            return "EXACT_STRONG"
        if any(exact):
            return "MIXED"
        if any(component):
            return "COMPONENT_STRONG"
        return "WEAK_OR_ABSENT"

    matrix_rows: list[dict[str, object]] = []
    for row in adjudication["all23_rows"]:
        alias = str(row["maize_family_alias"])
        matrix_rows.append(
            {
                "maize_family_alias": alias,
                "maize_gene_id": row["maize_gene_id"],
                "sorghum_candidate_ids": row["sorghum_candidate_ids"],
                "reciprocal_best_hit_status": row["reciprocal_best_hit_status"],
                "m09c_final_classification": row["m09c_final_classification"],
                "m09c_primary_state": matrix_state(alias, "m09c_primary"),
                "m09c_full_state": matrix_state(alias, "m09c_full"),
                "r1_state": matrix_state(alias, "r1"),
                "r2_state": matrix_state(alias, "r2"),
                "r3_state": matrix_state(alias, "r3"),
                "r4_state": matrix_state(alias, "r4"),
                "final_computational_classification": row["final_computational_classification"],
                "composition_and_mask_caveats": row["composition_and_mask_caveats"],
            }
        )
    workflow.write_tsv(
        figure_dir / "M09D_ROBUSTNESS_MATRIX_SOURCE.tsv",
        matrix_rows,
        list(matrix_rows[0]),
    )
    matrix_artifacts = figures.render_robustness_matrix(
        matrix_rows,
        figure_dir / "M09D_ALL_23_ROBUSTNESS_MATRIX",
        "M09D — 23-gene cross-model phylogenetic robustness matrix",
    )
    matrix_validation = figures.verify_matrix_artifacts(matrix_rows, matrix_artifacts)
    matrix_validation["artifacts"] = [
        str(Path(path).relative_to(staging))
        for path in matrix_validation["artifacts"]
    ]
    workflow.write_tsv(
        figure_dir / "M09D_FIGURE_SOURCE.tsv", trace_rows, list(trace_rows[0])
    )
    result = {
        "status": "PASS_PROGRAMMATIC_QC_PENDING_VISUAL_ATTESTATION",
        "phylogeny_figure_sets": len(validations),
        "robustness_matrix_figure_sets": 1,
        "figure_sets": len(validations) + 1,
        "validation": validations,
        "matrix_validation": matrix_validation,
        "trace_rows": len(trace_rows),
        "matrix_rows": len(matrix_rows),
    }
    workflow.atomic_write_json(figure_dir / "M09D_FIGURE_QC.json", result)
    return result


def _run_test_suite(command: list[str], cwd: Path, output: Path) -> dict[str, object]:
    process = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    text = process.stdout + process.stderr
    workflow.atomic_write_text(output, text)
    match = re.search(r"Ran\s+(\d+)\s+tests?", text)
    count = int(match.group(1)) if match else 0
    if process.returncode != 0 or count == 0 or not text.rstrip().endswith("OK"):
        raise RuntimeError(f"test suite failed: {' '.join(command)}; see {output}")
    return {"command": command, "exit_code": process.returncode, "tests": count, "status": "PASS"}


def run_regression_tests(project_root: Path, staging: Path) -> dict[str, object]:
    report_dir = staging / "PHYLOGENY/TEST_REPORTS"
    report_dir.mkdir(parents=True, exist_ok=True)
    m09d = _run_test_suite(
        [str(workflow.PYTHON), "-m", "unittest", "discover", "-v", "-s", str(workflow.TOOL_DIR / "tests")],
        project_root,
        report_dir / "M09D_UNITTEST.txt",
    )
    m09c = _run_test_suite(
        [
            str(workflow.PYTHON),
            "-m",
            "unittest",
            "discover",
            "-v",
            "-s",
            str(
                project_root
                / "11_MANUSCRIPT/UBUNTU_MASTER_TOOLS/M09C_INDEPENDENT_ORTHOLOGY_PHYLOGENY_R1/tests"
            ),
        ],
        project_root,
        report_dir / "M09C_INHERITED_UNITTEST.txt",
    )
    result = {"status": "PASS", "m09d": m09d, "m09c_inherited": m09c}
    workflow.atomic_write_json(report_dir / "M09D_TEST_RESULTS.json", result)
    return result


def _copy_reproducibility_code(project_root: Path, staging: Path) -> None:
    destination = staging / "PHYLOGENY/CODE"
    destination.mkdir(parents=True, exist_ok=True)
    for name in (
        "m09d_core.py",
        "m09d_workflow.py",
        "m09d_figures.py",
        "m09d_finalize.py",
        "m09d_validate.py",
        "MASTER_RUN_M09D.sh",
        "README.md",
    ):
        source = workflow.TOOL_DIR / name
        if not source.is_file():
            raise FileNotFoundError(f"reproducibility code is absent: {source}")
        shutil.copy2(source, destination / name)
    tests_destination = destination / "tests"
    tests_destination.mkdir(exist_ok=True)
    for source in sorted((workflow.TOOL_DIR / "tests").glob("test_*.py")):
        shutil.copy2(source, tests_destination / source.name)
    documentation = staging / "DOCUMENTATION"
    documentation.mkdir(exist_ok=True)
    shutil.copy2(
        project_root
        / "docs/superpowers/specs/2026-11-07-m09d-composition-robustness-design.md",
        documentation / "M09D_APPROVED_DESIGN.md",
    )
    shutil.copy2(
        project_root
        / "docs/superpowers/plans/2026-11-07-m09d-composition-robustness.md",
        documentation / "M09D_IMPLEMENTATION_PLAN.md",
    )


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _software_version(command: list[str]) -> str:
    process = subprocess.run(command, text=True, capture_output=True, check=False)
    output = (process.stdout + process.stderr).strip()
    if process.returncode != 0 or not output:
        raise RuntimeError(f"software version command failed: {' '.join(command)}")
    return " ".join(output.split())


def write_software_and_parameters(staging: Path, model_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    import Bio
    import matplotlib

    rows: list[dict[str, object]] = [
        {
            "record_type": "SOFTWARE",
            "analysis": "ALL",
            "software_or_parameter": "IQ-TREE",
            "value": _software_version([str(workflow.IQTREE3), "--version"]),
            "path_or_scope": str(workflow.IQTREE3),
        },
        {
            "record_type": "SOFTWARE",
            "analysis": "ALL",
            "software_or_parameter": "Python",
            "value": sys.version.replace("\n", " "),
            "path_or_scope": sys.executable,
        },
        {
            "record_type": "SOFTWARE",
            "analysis": "FIGURES_AND_VALIDATION",
            "software_or_parameter": "Biopython",
            "value": Bio.__version__,
            "path_or_scope": "Python_package",
        },
        {
            "record_type": "SOFTWARE",
            "analysis": "FIGURES",
            "software_or_parameter": "Matplotlib",
            "value": matplotlib.__version__,
            "path_or_scope": "Python_package",
        },
    ]
    for row in model_rows:
        analysis = str(row["analysis"])
        rows.extend(
            [
                {
                    "record_type": "PARAMETER",
                    "analysis": analysis,
                    "software_or_parameter": "declared_model",
                    "value": row["declared_model"],
                    "path_or_scope": row["matrix"],
                },
                {
                    "record_type": "RESULT",
                    "analysis": analysis,
                    "software_or_parameter": "selected_model",
                    "value": row["selected_model"],
                    "path_or_scope": row["matrix"],
                },
            ]
        )
    for parameter, value in (
        ("maximum_threads", 2),
        ("independent_ml_searches", 5),
        ("ultrafast_bootstrap_replicates", 1000),
        ("sh_alrt_replicates", 1000),
        ("bootstrap_nni", "ENABLED"),
        ("rooting", "UNROOTED"),
        ("strong_support_threshold", "SH-aLRT>=80_AND_UFBoot>=95"),
    ):
        rows.append(
            {
                "record_type": "PARAMETER",
                "analysis": "R1_R2_R3_R4",
                "software_or_parameter": parameter,
                "value": value,
                "path_or_scope": "PRE_REGISTERED",
            }
        )
    workflow.write_tsv(
        staging / "M09D_SOFTWARE_AND_PARAMETERS.tsv", rows, list(rows[0])
    )
    return rows


def write_resource_and_performance_summary(staging: Path) -> list[dict[str, object]]:
    memory: dict[str, int] = {}
    for raw in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if ":" not in raw:
            continue
        key, value = raw.split(":", 1)
        token = value.strip().split()[0]
        if token.isdigit():
            memory[key] = int(token)
    disk = shutil.disk_usage(staging)
    rows: list[dict[str, object]] = [
        {"record_type": "RESOURCE_SNAPSHOT", "step": "logical_cpus_available", "wall_seconds": "NA", "max_rss_kb": "NA", "value": os.cpu_count() or "UNKNOWN", "unit": "CPUS"},
        {"record_type": "RESOURCE_SNAPSHOT", "step": "MemTotal", "wall_seconds": "NA", "max_rss_kb": "NA", "value": memory.get("MemTotal", "UNKNOWN"), "unit": "KB"},
        {"record_type": "RESOURCE_SNAPSHOT", "step": "MemAvailable", "wall_seconds": "NA", "max_rss_kb": "NA", "value": memory.get("MemAvailable", "UNKNOWN"), "unit": "KB"},
        {"record_type": "RESOURCE_SNAPSHOT", "step": "SwapTotal", "wall_seconds": "NA", "max_rss_kb": "NA", "value": memory.get("SwapTotal", "UNKNOWN"), "unit": "KB"},
        {"record_type": "RESOURCE_SNAPSHOT", "step": "release_filesystem_free", "wall_seconds": "NA", "max_rss_kb": "NA", "value": disk.free, "unit": "BYTES"},
    ]
    for row in workflow.read_tsv(staging / "PERFORMANCE_LOG.tsv"):
        rows.append(
            {
                "record_type": "EXECUTION_PERFORMANCE",
                "step": row["step"],
                "wall_seconds": row["wall_seconds"],
                "max_rss_kb": row["max_rss_kb"],
                "value": row["exit_code"],
                "unit": "EXIT_CODE",
            }
        )
    workflow.write_tsv(
        staging / "M09D_RESOURCE_AND_PERFORMANCE_SUMMARY.tsv",
        rows,
        list(rows[0]),
    )
    return rows


def _classification_counts(rows: Iterable[Mapping[str, object]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                str(row["final_computational_classification"]) for row in rows
            ).items()
        )
    )


def _markdown_table(headers: list[str], rows: Iterable[Iterable[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(str(value).replace("|", "\\|").replace("\n", " ") for value in row)
            + " |"
        )
    return "\n".join(lines)


def _scientific_status(
    adjudication: Mapping[str, object], symmetry_rows: list[dict[str, str]], model_rows: list[dict[str, object]]
) -> tuple[str, list[str]]:
    counts = _classification_counts(adjudication["all23_rows"])
    reasons: list[str] = []
    if counts.get("UNRESOLVED", 0):
        reasons.append(f"{counts['UNRESOLVED']} maize outcomes remain UNRESOLVED")
    if counts.get("CONTRADICTED", 0):
        reasons.append(f"{counts['CONTRADICTED']} maize outcomes are CONTRADICTED")
    caveat_rows = [
        row
        for row in adjudication["all23_rows"]
        if str(row["composition_and_mask_caveats"]) != "NONE"
    ]
    if caveat_rows:
        reasons.append(f"{len(caveat_rows)} maize outcomes retain direct composition/masking caveats")
    significant = [row for row in symmetry_rows if row["q_lt_0.05"] == "YES"]
    if significant:
        reasons.append(f"{len(significant)} alignment-by-test-family symmetry summaries have BH q<0.05")
    if any("overfitting" in str(row["warnings"]).lower() for row in model_rows):
        reasons.append("the R1 profile-mixture fit reported near-zero weights and an overfitting warning")
    reasons.append("IQ-TREE did not export individual taxon-pair symmetry p-values for pair-level caveat assignment")
    return ("REVIEW" if reasons else "PASS"), reasons


def write_scientific_and_manuscript_reports(
    project_root: Path,
    staging: Path,
    adjudication: Mapping[str, object],
    model_rows: list[dict[str, object]],
    figure_qc: Mapping[str, object],
    test_results: Mapping[str, object],
) -> dict[str, object]:
    composition = workflow.read_tsv(staging / "M09D_COMPOSITION_DIAGNOSTICS.tsv")
    sensitivity = workflow.read_tsv(staging / "M09D_ALIGNMENT_SENSITIVITY.tsv")
    symmetry = workflow.read_tsv(
        staging / "PHYLOGENY/SYMMETRY/M09D_SYMMETRY_TESTS.tsv"
    )
    source_audit = json.loads((staging / "SOURCE_INTEGRITY_AUDIT.json").read_text(encoding="utf-8"))
    all23 = list(adjudication["all23_rows"])
    seven = list(adjudication["seven_rows"])
    relationships = list(adjudication["relationship_rows"])
    class_counts = _classification_counts(all23)
    pair_class_counts = _classification_counts(relationships)
    status, review_reasons = _scientific_status(adjudication, symmetry, model_rows)
    primary_failures = list(adjudication["primary_composition_failures"])
    full_failures = list(adjudication["full_length_composition_failures"])
    gap_rows = [row for row in composition if row["greater_than_50_percent_gap_or_ambiguity"] == "YES"]
    gap_taxa = sorted({row["taxon"] for row in gap_rows})
    long_branches = [row for row in composition if row["long_branch_flag"] == "YES"]
    domain_concern_taxa = sorted(
        {row["taxon"] for row in composition if row["domain_completeness_concern"] != "NONE"}
    )
    significant_symmetry = [row for row in symmetry if row["q_lt_0.05"] == "YES"]
    focal = next(row for row in relationships if row["maize_family_alias"] == "ZmNRG2_001")
    class_text = ", ".join(
        f"{count} as {name}" for name, count in class_counts.items()
    )
    pair_class_text = ", ".join(
        f"{count} as {name}" for name, count in pair_class_counts.items()
    )

    limitations = f"""# M09D Scientific Limitations

M09D status: **{status}**. This stage provides composition-aware computational robustness evidence only. It does not establish experimental orthology, biochemical function, functional conservation, or stress-response equivalence.

## Model and composition limitations

- The unchanged dual-domain analysis failed the original IQ-TREE composition chi-square test for {len(primary_failures)} taxa: {', '.join(primary_failures)}.
- The unchanged full-length analysis failed that test for {len(full_failures)} taxa: {', '.join(full_failures)}.
- IQ-TREE symmetry summaries yielded {len(significant_symmetry)} BH-significant alignment/partition-by-test-family results. Individual taxon-pair p-values were not exported by IQ-TREE 3.1.2, so pair-specific symmetry failure cannot be assigned from these outputs.
- R1 uses a bounded C20 profile mixture. It addresses across-site heterogeneity but is not a full nonstationary across-lineage model; its near-zero mixture weights triggered an overfitting warning.
- Dayhoff-6 recoding reduces compositional detail and is a sensitivity analysis, not a preferred biological model.
- AIC/BIC comparisons are valid only for fits to the identical character matrix. Masked and recoded analyses are explicitly non-comparable by information criteria.

## Alignment and sampling limitations

- {len(gap_rows)} alignment rows, representing {len(gap_taxa)} unique taxa, exceed 50% gaps or ambiguity. The row/taxon distinction is retained because Sb_GeneID_8076242 is flagged in both matrices.
- Five source-locked domain-completeness concerns are retained across {len(domain_concern_taxa)} taxa: {', '.join(domain_concern_taxa)}. No sequence was removed or padded.
- The analysis contains exactly 23 maize and 17 Sorghum representatives and no defensible outgroup; every tree is therefore unrooted. Direction, ancestral state, and duplication timing are not inferred.
- Component-level support can reflect duplicated or paralog-associated collinearity and is not equivalent to an exact one-to-one relationship.

## Classification limitations

- Reciprocal rank-1 homology, domain architecture, collinearity and phylogenetic support are jointly interpreted. Bootstrap support alone is never treated as proof of orthology.
- Seven NO_RBH cases are ineligible for either supported-candidate class by the pre-registered rules.
- The release remains AUTHOR_REVIEW. No content was promoted into `MASTER_CURRENT.docx`.

## Reasons for REVIEW

""" + "\n".join(f"- {reason}." for reason in review_reasons) + "\n"
    workflow.atomic_write_text(staging / "M09D_SCIENTIFIC_LIMITATIONS.md", limitations)

    selected_models = {str(row["analysis"]): str(row["selected_model"]) for row in model_rows}
    methods_results = f"""# M09D Author-Review Manuscript Evidence Draft

**Status:** AUTHOR_REVIEW_ONLY; scientific release status **{status}**. This text is evidence preparation and was not inserted into the canonical manuscript.

## Materials and Methods

The composition-aware robustness analysis reused the source-locked M09C set of 40 representatives (23 *Zea mays* and 17 *Sorghum bicolor*) and the original unrooted 506-residue partitioned PF04783/PF04782 and 1,611-residue full-length alignments. No sequence was substituted, removed, padded, or re-aligned. Per-taxon amino-acid composition, gap/ambiguity fraction, Jensen-Shannon distance from pooled composition, source domain coverage and terminal branch length were reconciled to the original IQ-TREE composition tests. IQ-TREE symmetry, marginal-symmetry and internal-symmetry summaries were evaluated with Benjamini-Hochberg correction separately by alignment and test family; IQ-TREE did not export individual pairwise p-values.

Four pre-registered, unrooted sensitivity analyses retained all 40 taxa. R1 fitted LG+C20+F+G4 to the unchanged dual-domain alignment. R2 applied a taxon-blind, partition-wise occupancy threshold of at least 20 standard amino acids among 40 taxa and retained {next(row['retained_columns'] for row in sensitivity if row['analysis'] == 'R2')} columns; ModelFinder/partition merging selected {selected_models['R2']}. R3 applied the same threshold to the full-length alignment, retained {next(row['retained_columns'] for row in sensitivity if row['analysis'] == 'R3')} columns and selected {selected_models['R3']}. R4 recoded the unchanged dual-domain alignment into the pre-specified Dayhoff-6 groups and fitted GTR+FO+G4. Each new analysis used five independent maximum-likelihood searches, 1,000 ultrafast bootstrap replicates, 1,000 SH-aLRT replicates, bootstrap NNI optimization, fixed seeds and at most two threads. Strong bipartition support was defined a priori as SH-aLRT ≥80 and UFBoot ≥95.

Relationships were re-adjudicated with root-invariant bipartitions using reciprocal rank-1 evidence, ordered dual-domain evidence, M09B inter-species collinearity components, exact-pair support, component support, alternative strongly supported pairs, composition warnings and occupancy sensitivity. Component support was not promoted to exact one-to-one orthology.

## Results

All four robustness analyses completed with 40 taxa, five searches and the requested 1,000/1,000 support replicates. R2 retained 363/506 columns (PF04783: 59; PF04782: 304), whereas R3 retained 735/1,611 columns. The original dual-domain composition test failed for {len(primary_failures)} taxa and the full-length test failed for {len(full_failures)} taxa. Across the two alignments, {len(gap_rows)} flagged alignment rows represented {len(gap_taxa)} unique taxa with greater than 50% gaps or ambiguity. No taxon was deleted.

At the 23-maize-gene roll-up level, M09D classified {class_text}. At the 24 separately retained pair-hypothesis level, the corresponding counts were {pair_class_text}. The focal ZmNRG2_001–Sb_GeneID_8069100 relationship was classified **{focal['final_computational_classification']}**; its exact-pair support was M09C-primary {focal['m09c_primary_exact_support_shalrt_ufboot']}, M09C-full {focal['m09c_full_exact_support_shalrt_ufboot']}, R1 {focal['r1_exact_support_shalrt_ufboot']}, R2 {focal['r2_exact_support_shalrt_ufboot']}, R3 {focal['r3_exact_support_shalrt_ufboot']} and R4 {focal['r4_exact_support_shalrt_ufboot']} (SH-aLRT/UFBoot).

All seven original NO_RBH maize cases remained explicitly separated from reciprocal-best-hit candidates. Their final outcomes are listed in `M09D_SEVEN_NO_RBH_CASES.tsv`; none was interpreted as experimentally proven orthology.

## Discussion

Concordance across unchanged, occupancy-masked, profile-mixture and reduced-state analyses strengthens some computational relationships, whereas changes in exact-pair support identify model- or missing-data-sensitive cases. The composition and symmetry diagnostics demonstrate that the homogeneous stationary models are imperfect descriptions of these data, but those warnings alone do not prove that a particular topology is wrong. Consequently, relationships carrying composition, masking or alternative-topology caveats remain qualified. Possible co-orthology denotes a duplication-aware component-level hypothesis and not an exact one-to-one assignment.

## Scientific limitations

See `M09D_SCIENTIFIC_LIMITATIONS.md`. Experimental validation, broader taxon sampling with a defensible outgroup and composition-aware nonstationary inference remain outside this stage.
"""
    workflow.atomic_write_text(
        staging / "M09D_MANUSCRIPT_METHODS_RESULTS_DRAFT.md", methods_results
    )

    legends = """# M09D Author-Review Figure Legends

**Figure M09D-1 (R1).** Unrooted maximum-likelihood phylogram of the 23 maize and 17 Sorghum source-locked representatives inferred from the unchanged 506-column PF04783/PF04782 dual-domain alignment under LG+C20+F+G4. Internal labels report SH-aLRT/UFBoot percentages. Branch lengths are substitutions per site. The rectangular orientation does not imply a root.

**Figure M09D-2 (R2).** Unrooted maximum-likelihood occupancy-masked dual-domain sensitivity tree. All taxa were retained; sites were retained partition-wise when at least 20 of 40 sequences contained a standard amino acid. The exact selected partition scheme is reported in the IQ-TREE output. Symbols and support labels follow Figure M09D-1.

**Figure M09D-3 (R3).** Unrooted maximum-likelihood full-length occupancy-masked sensitivity tree. The taxon-blind mask retained 735 of 1,611 columns and all 40 taxa. Symbols and support labels follow Figure M09D-1.

**Figure M09D-4 (R4).** Unrooted Dayhoff-6 recoding sensitivity tree inferred from the unchanged 506-column dual-domain alignment under GTR+FO+G4. Reduced-state branch lengths are model-specific and should not be compared numerically with amino-acid-tree branch lengths.

**Figure M09D-5.** Cross-model robustness matrix for all 23 maize family members. E denotes a strongly supported exact pair, C denotes strong support for the full M09B collinearity component, and a centered dot denotes weak or absent support under the pre-registered SH-aLRT ≥80 and UFBoot ≥95 rule. Multi-partner rows are summarized conservatively; `M09D_ORTHOLOGY_EVIDENCE_MATRIX.tsv` preserves all 24 individual pair hypotheses.
"""
    workflow.atomic_write_text(staging / "M09D_MANUSCRIPT_LEGENDS.md", legends)

    summary = f"""# M09D Execution Summary

## Outcome

**{status}** — all four pre-registered robustness analyses and independent re-adjudication computations completed, but material composition/model caveats remain. This is not a forced PASS.

## Source integrity

- M09C manifest: {source_audit['m09c_manifest_entries']}/120 verified; SHA256 `{source_audit['m09c_manifest_sha256']}`.
- M09B manifest: {source_audit['m09b_manifest_entries']}/61 verified; SHA256 `{source_audit['m09b_manifest_sha256']}`.
- Source-locked taxa: 23 maize + 17 Sorghum = 40.
- Canonical master SHA256: `{source_audit['canonical_master_sha256']}` (unchanged).

## Analyses executed

- R1 unchanged dual-domain LG+C20+F+G4.
- R2 partition-wise 20/40 occupancy mask plus ModelFinder/partition merging.
- R3 full-length 20/40 occupancy mask plus ModelFinder.
- R4 literal Dayhoff-6 recoding under GTR+FO+G4.
- All used 2 threads, five ML searches, 1,000 UFBoot, 1,000 SH-aLRT, `-bnni`, fixed seeds and unrooted interpretation.
- Root-invariant adjudication of 16 RBH candidate pairs and eight pair hypotheses across the seven NO_RBH maize cases.

## Classification totals

23-gene roll-up: {json.dumps(class_counts, sort_keys=True)}.

24 pair hypotheses: {json.dumps(pair_class_counts, sort_keys=True)}.

## Not executed

No new MAFFT, DIAMOND, HMMER, MCScanX or sequence extraction; no taxon deletion; no rooted inference; no unsupported outgroup; no functional-conservation claim; no manuscript promotion; no Git push.

## Review reasons

""" + "\n".join(f"- {reason}." for reason in review_reasons) + "\n"
    workflow.atomic_write_text(staging / "M09D_EXECUTION_SUMMARY.md", summary)

    index_rows = [
        {"artifact": "M09D_QC.json", "description": "Machine-readable terminal scientific status and counts", "review_state": "AUTHOR_REVIEW"},
        {"artifact": "M09D_COMPOSITION_DIAGNOSTICS.tsv", "description": "Per-alignment, per-taxon composition, missingness, domain and branch diagnostics", "review_state": "VALIDATED"},
        {"artifact": "M09D_ALIGNMENT_SENSITIVITY.tsv", "description": "Pre-registered transform and site-retention summary", "review_state": "VALIDATED"},
        {"artifact": "M09D_MODEL_COMPARISON.tsv", "description": "Model, likelihood and information-criterion provenance with comparability guardrails", "review_state": "VALIDATED"},
        {"artifact": "M09D_RESOURCE_AND_PERFORMANCE_SUMMARY.tsv", "description": "WSL resource snapshot and actual timed-step wall/RSS records", "review_state": "VALIDATED"},
        {"artifact": "M09D_ORTHOLOGY_EVIDENCE_MATRIX.tsv", "description": "All 24 pair hypotheses with six-tree exact/component support", "review_state": "AUTHOR_REVIEW"},
        {"artifact": "M09D_ALL_23_FINAL_REVIEW.tsv", "description": "Conservative 23-maize-gene roll-up", "review_state": "AUTHOR_REVIEW"},
        {"artifact": "M09D_SEVEN_NO_RBH_CASES.tsv", "description": "Seven predetermined NO_RBH outcomes", "review_state": "AUTHOR_REVIEW"},
        {"artifact": "PHYLOGENY/TREES/", "description": "Actual R1-R4 IQ-TREE alignments, treefiles, support trees, reports and logs", "review_state": "VALIDATED"},
        {"artifact": "PHYLOGENY/FIGURES/", "description": "Five traceable figure sets in SVG, vector PDF and 600-dpi PNG plus source TSVs", "review_state": "AUTHOR_REVIEW"},
        {"artifact": "M09D_MANUSCRIPT_METHODS_RESULTS_DRAFT.md", "description": "Methods, Results and Discussion draft", "review_state": "AUTHOR_REVIEW_ONLY"},
        {"artifact": "M09D_MANUSCRIPT_LEGENDS.md", "description": "Figure legends", "review_state": "AUTHOR_REVIEW_ONLY"},
        {"artifact": "M09D_SCIENTIFIC_LIMITATIONS.md", "description": "Explicit model, alignment and inference limitations", "review_state": "AUTHOR_REVIEW_ONLY"},
        {"artifact": "PHYLOGENY/TEST_REPORTS/", "description": "Fresh M09D and inherited M09C regression evidence", "review_state": "VALIDATED"},
        {"artifact": "PHYLOGENY/CODE/", "description": "Versioned workflow, finalizer, validator, figure code and tests", "review_state": "REPRODUCIBILITY"},
    ]
    workflow.write_tsv(
        staging / "M09D_SUPPLEMENTARY_INDEX.tsv", index_rows, list(index_rows[0])
    )

    state = json.loads((staging / "RUN_STATE.json").read_text(encoding="utf-8"))
    final_release = state["final_release"]
    reproducibility = f"""# M09D Reproducibility Commands

The commands below use only the immutable M09C/M09B sources and the versioned M09D code. All IQ-TREE commands are emitted in per-analysis `COMMAND.json` files. The maximum thread count is two.

```bash
PROJECT=/home/hammadse/AYESHA_ZmNRL1_WORKFLOW
TOOL=$PROJECT/11_MANUSCRIPT/UBUNTU_MASTER_TOOLS/M09D_COMPOSITION_ROBUSTNESS_R1
bash $TOOL/MASTER_RUN_M09D.sh $PROJECT
```

Resume an interrupted, unpublished staging directory without replacing successful IQ-TREE outputs:

```bash
bash {workflow.TOOL_DIR / 'MASTER_RUN_M09D.sh'} {project_root} --resume /absolute/path/to/.STAGING_<token>
```

Verify the published release independently:

```bash
{workflow.PYTHON} {workflow.TOOL_DIR / 'm09d_validate.py'} {project_root} {final_release} --require-manifest
```

The release was constructed atomically only after source verification, tree/output reconciliation, fresh regression tests, programmatic figure checks, visual attestation and producer-independent validation.
"""
    workflow.atomic_write_text(staging / "M09D_REPRODUCIBILITY.md", reproducibility)

    qc = {
        "schema": "M09D_QC_V1",
        "status": status,
        "scientific_status_label": "COMPOSITION_AWARE_ROBUSTNESS_COMPUTED_AUTHOR_REVIEW_REQUIRED" if status == "REVIEW" else "COMPOSITION_AWARE_ROBUSTNESS_VALIDATED",
        "completed_utc": _utc_now(),
        "source_integrity": source_audit,
        "input_taxa": {"total": 40, "maize": 23, "sorghum": 17},
        "analyses": {
            row["analysis"]: {
                "matrix": row["matrix"],
                "taxa": row["taxa"],
                "columns": row["columns"],
                "declared_model": row["declared_model"],
                "selected_model": row["selected_model"],
                "warnings": row["warnings"],
            }
            for row in model_rows
        },
        "composition": {
            "primary_failures": primary_failures,
            "primary_failure_count": len(primary_failures),
            "full_length_failures": full_failures,
            "full_length_failure_count": len(full_failures),
            "greater_than_50_percent_gap_or_ambiguity_rows": len(gap_rows),
            "greater_than_50_percent_gap_or_ambiguity_unique_taxa": len(gap_taxa),
            "long_branch_rows": len(long_branches),
            "domain_completeness_concern_taxa": domain_concern_taxa,
            "symmetry_significant_summary_rows": len(significant_symmetry),
            "pair_specific_symmetry_p_values": "NOT_EXPORTED_BY_IQTREE_3.1.2",
        },
        "alignment_sensitivity": {row["analysis"]: row for row in sensitivity},
        "adjudication": {
            "pair_hypotheses": len(relationships),
            "pair_classification_counts": pair_class_counts,
            "all_23_classification_counts": class_counts,
            "focal_zmnrg2_001": focal,
            "seven_no_rbh": [
                {
                    "maize_family_alias": row["maize_family_alias"],
                    "sorghum_candidate_ids": row["sorghum_candidate_ids"],
                    "final_computational_classification": row["final_computational_classification"],
                }
                for row in seven
            ],
        },
        "figures": dict(figure_qc),
        "tests": dict(test_results),
        "review_reasons": review_reasons,
        "independent_validation": {"status": "PENDING_VISUAL_QA_AND_PREPUBLICATION_VALIDATION"},
        "canonical_master_sha256": workflow.sha256_file(
            project_root / "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/MASTER_CURRENT.docx"
        ),
        "canonical_master_modified": False,
        "manuscript_promotion": "NOT_PERFORMED_AUTHOR_REVIEW_ONLY",
    }
    workflow.atomic_write_json(staging / "M09D_QC.json", qc)
    build_summary = {
        "status": "PASS_BUILD_COMPLETE_PENDING_VISUAL_QA",
        "scientific_status": status,
        "classification_counts": class_counts,
        "pair_classification_counts": pair_class_counts,
        "review_reasons": review_reasons,
        "figure_sets": figure_qc["figure_sets"],
        "m09d_tests": test_results["m09d"]["tests"],
        "m09c_tests": test_results["m09c_inherited"]["tests"],
    }
    workflow.atomic_write_json(staging / "M09D_FINAL_BUILD_SUMMARY.json", build_summary)
    return {"status": status, "qc": qc, "summary": build_summary}


def _assert_all_analyses_complete(staging: Path) -> None:
    expected_taxa = set(
        core.parse_fasta(
            staging / "PHYLOGENY/INPUTS/M09D_M09C_DUAL_DOMAIN_40_UNCHANGED.faa",
            aligned=True,
        )
    )
    core.validate_taxa(expected_taxa)
    for analysis in NEW_ANALYSES:
        configuration = workflow._analysis_configuration(staging, analysis)
        prefix = Path(configuration["prefix"])
        if not workflow.iqtree_artifacts_complete(prefix, expected_taxa):
            raise ValueError(f"mandatory IQ-TREE analysis is incomplete: {analysis}")
        metrics_path = prefix.parent / "METRICS.json"
        if not metrics_path.is_file():
            raise ValueError(f"mandatory metrics are absent: {analysis}")
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        if not str(metrics.get("status", "")).startswith("PASS_COMPLETE"):
            raise ValueError(f"mandatory metrics are not complete: {analysis}")


def build_release(project_root: Path, staging: Path) -> dict[str, object]:
    root = project_root.resolve()
    staging = workflow.safe_staging_path(root / workflow.RELEASE_PARENT_RELATIVE, staging)
    workflow._verify_prepare_checkpoint(staging)
    workflow._verify_live_code_lock(staging)
    fresh_audit = workflow.source_integrity_audit(root)
    original_audit = json.loads((staging / "SOURCE_INTEGRITY_AUDIT.json").read_text(encoding="utf-8"))
    for key in (
        "m09c_manifest_entries",
        "m09c_manifest_sha256",
        "m09b_manifest_entries",
        "m09b_manifest_sha256",
        "canonical_master_sha256",
        "candidate_pairs",
        "no_rbh_cases",
        "family_anchor_rows",
    ):
        if fresh_audit[key] != original_audit[key]:
            raise ValueError(f"source integrity changed since prepare: {key}")
    _assert_all_analyses_complete(staging)
    adjudication = adjudicate_relationships(staging)
    model_rows = build_model_comparison(root, staging)
    figure_qc = generate_figures(staging, adjudication, model_rows)
    _copy_reproducibility_code(root, staging)
    test_results = run_regression_tests(root, staging)
    software_rows = write_software_and_parameters(staging, model_rows)
    resource_rows = write_resource_and_performance_summary(staging)
    reports = write_scientific_and_manuscript_reports(
        root, staging, adjudication, model_rows, figure_qc, test_results
    )
    state = json.loads((staging / "RUN_STATE.json").read_text(encoding="utf-8"))
    state.update(
        {
            "status": "BUILD_COMPLETE_PENDING_VISUAL_QA",
            "scientific_status": reports["status"],
            "updated_utc": _utc_now(),
        }
    )
    workflow.atomic_write_json(staging / "RUN_STATE.json", state)
    result = {
        "status": "PASS_BUILD_COMPLETE_PENDING_VISUAL_QA",
        "scientific_status": reports["status"],
        "staging": str(staging),
        "figure_directory": str(staging / "PHYLOGENY/FIGURES"),
        "figure_sets": figure_qc["figure_sets"],
        "m09d_tests": test_results["m09d"]["tests"],
        "m09c_inherited_tests": test_results["m09c_inherited"]["tests"],
        "software_parameter_rows": len(software_rows),
        "resource_performance_rows": len(resource_rows),
    }
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return result


def _run_independent_validator(
    project_root: Path,
    release: Path,
    output: Path,
    *,
    require_manifest: bool,
) -> dict[str, object]:
    command = [
        str(workflow.PYTHON),
        str(workflow.TOOL_DIR / "m09d_validate.py"),
        str(project_root),
        str(release),
        "--output",
        str(output),
    ]
    if require_manifest:
        command.append("--require-manifest")
    process = subprocess.run(command, text=True, capture_output=True, check=False)
    if process.returncode != 0 or not output.is_file():
        raise RuntimeError(
            "producer-independent M09D validation failed: "
            + (process.stdout + process.stderr)[-4000:]
        )
    result = json.loads(output.read_text(encoding="utf-8"))
    if result.get("status") != "PASS":
        raise RuntimeError(f"producer-independent M09D validation did not PASS: {result}")
    return result


def _write_complete_manifest(staging: Path) -> int:
    manifest = staging / "SHA256SUMS.txt"
    if manifest.exists():
        raise FileExistsError("refusing to overwrite an existing M09D manifest")
    files = [
        path
        for path in sorted(staging.rglob("*"))
        if path.is_file() and path != manifest
    ]
    if not files:
        raise ValueError("cannot publish an empty M09D release")
    lines = [
        f"{workflow.sha256_file(path)}  {path.relative_to(staging).as_posix()}"
        for path in files
    ]
    workflow.atomic_write_text(manifest, "\n".join(lines) + "\n")
    verified = workflow.verify_manifest(staging, manifest)
    if verified != len(files):
        raise ValueError(f"manifest verification count mismatch: {verified} != {len(files)}")
    return verified


def _append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.is_file() else ""
    if marker in existing:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="") as handle:
        if existing and not existing.endswith("\n"):
            handle.write("\n")
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def _update_external_ledgers(
    project_root: Path,
    final_release: Path,
    scientific_status: str,
    manifest_entries: int,
    manifest_sha256: str,
) -> None:
    marker = final_release.name
    timestamp = _utc_now()
    ledger_entry = f"""

## {timestamp} — M09D composition-aware robustness release

- Release: `{final_release}`
- Scientific status: **{scientific_status}**
- Manifest: {manifest_entries} entries; SHA256 `{manifest_sha256}`
- Actual analyses: R1 LG+C20+F+G4, R2 partition-wise dual-domain occupancy mask, R3 full-length occupancy mask, R4 Dayhoff-6; all 40 taxa, five searches, 1,000 UFBoot, 1,000 SH-aLRT, at most two threads.
- Manuscript state: AUTHOR_REVIEW package only; canonical master not modified.
- Marker: `{marker}`
"""
    progress_entry = f"""

## M09D — {timestamp}

- `{marker}` published atomically with scientific status **{scientific_status}**.
- Source-locked M09B/M09C evidence preserved; all 23 maize and 17 Sorghum representatives retained.
- Composition-aware orthology re-adjudication and author-review manuscript evidence are available in `{final_release}`.
- `MASTER_CURRENT.docx` remains unchanged.
"""
    _append_once(project_root / "logs_QC/M09D_EXECUTION_LEDGER.md", marker, ledger_entry)
    _append_once(project_root / "logs_QC/MASTER_PROGRESS.md", marker, progress_entry)


def publish_release(project_root: Path, staging: Path) -> dict[str, object]:
    root = project_root.resolve()
    staging = workflow.safe_staging_path(root / workflow.RELEASE_PARENT_RELATIVE, staging)
    workflow._verify_prepare_checkpoint(staging)
    workflow._verify_live_code_lock(staging)
    _assert_all_analyses_complete(staging)
    state = json.loads((staging / "RUN_STATE.json").read_text(encoding="utf-8"))
    if state.get("status") != "BUILD_COMPLETE_PENDING_VISUAL_QA":
        raise ValueError(f"M09D build is not ready for visual QA: {state.get('status')}")
    attestation_path = staging / "PHYLOGENY/FIGURES/M09D_VISUAL_QA_ATTESTATION.json"
    if not attestation_path.is_file():
        raise ValueError("visual-QA attestation is absent; fail closed")
    attestation = json.loads(attestation_path.read_text(encoding="utf-8"))
    if attestation.get("status") != "PASS_VISUAL_QA" or int(attestation.get("figure_sets_reviewed", 0)) != 5:
        raise ValueError("visual-QA attestation is incomplete or not PASS")
    figure_qc_path = staging / "PHYLOGENY/FIGURES/M09D_FIGURE_QC.json"
    figure_qc = json.loads(figure_qc_path.read_text(encoding="utf-8"))
    figure_qc["status"] = "PASS_PROGRAMMATIC_AND_VISUAL_QC"
    figure_qc["visual_attestation"] = "PHYLOGENY/FIGURES/M09D_VISUAL_QA_ATTESTATION.json"
    workflow.atomic_write_json(figure_qc_path, figure_qc)
    state.update({"status": "READY_FOR_ATOMIC_PUBLISH", "updated_utc": _utc_now()})
    workflow.atomic_write_json(staging / "RUN_STATE.json", state)

    temporary_validation = Path("/tmp") / f"M09D_PREVALIDATION_{staging.name}.json"
    _run_independent_validator(root, staging, temporary_validation, require_manifest=False)
    qc_path = staging / "M09D_QC.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["independent_validation"] = {
        "status": "PASS_PRE_MANIFEST",
        "validator": "PHYLOGENY/CODE/m09d_validate.py",
        "validated_utc": _utc_now(),
    }
    qc["figures"] = figure_qc
    qc["canonical_master_sha256"] = workflow.sha256_file(
        root / "11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/MASTER_CURRENT.docx"
    )
    if qc["canonical_master_sha256"] != workflow.EXPECTED_MASTER_SHA256:
        raise ValueError("canonical master changed before M09D publication")
    workflow.atomic_write_json(qc_path, qc)
    in_release_validation = staging / "PHYLOGENY/TEST_REPORTS/INDEPENDENT_VALIDATION_PRE_MANIFEST.json"
    _run_independent_validator(root, staging, in_release_validation, require_manifest=False)

    manifest_entries = _write_complete_manifest(staging)
    manifest_sha256 = workflow.sha256_file(staging / "SHA256SUMS.txt")
    final_release = Path(state["final_release"]).resolve()
    if final_release.exists():
        raise FileExistsError(f"final M09D release already exists: {final_release}")
    if final_release.parent != staging.parent or not final_release.name.startswith("RELEASE_"):
        raise ValueError("final M09D release path is not the registered sibling target")
    os.replace(staging, final_release)

    final_validation_path = Path("/tmp") / f"M09D_FINAL_VALIDATION_{final_release.name}.json"
    final_validation = _run_independent_validator(
        root, final_release, final_validation_path, require_manifest=True
    )
    _update_external_ledgers(
        root,
        final_release,
        str(qc["status"]),
        manifest_entries,
        manifest_sha256,
    )
    result = {
        "status": "PASS_ATOMIC_PUBLICATION",
        "scientific_status": qc["status"],
        "release": str(final_release),
        "manifest_entries": manifest_entries,
        "manifest_sha256": manifest_sha256,
        "independent_validation": final_validation,
        "independent_validation_report": str(final_validation_path),
        "canonical_master_sha256": qc["canonical_master_sha256"],
    }
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    parser.add_argument("command", choices=["build", "publish"])
    parser.add_argument("--staging", required=True, type=Path)
    arguments = parser.parse_args(argv)
    if arguments.command == "build":
        build_release(arguments.project_root, arguments.staging)
    else:
        publish_release(arguments.project_root, arguments.staging)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
