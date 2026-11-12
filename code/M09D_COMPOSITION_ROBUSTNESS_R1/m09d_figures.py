"""Traceable publication phylograms for M09D unrooted robustness trees."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo


MAIZE_COLOR = "#1B9E77"
SORGHUM_COLOR = "#D95F02"
BRANCH_COLOR = "#303030"
CLASS_COLORS = {
    "ROBUST_COMPUTATIONAL_CANDIDATE": "#2E7D32",
    "SUPPORTED_WITH_MODEL_CAVEAT": "#558B2F",
    "POSSIBLE_CO_ORTHOLOGY": "#1565C0",
    "UNRESOLVED": "#EF6C00",
    "CONTRADICTED": "#C62828",
    "CONTEXT_ONLY": "#BDBDBD",
}

MATRIX_STATE_COLORS = {
    "EXACT_STRONG": "#2166AC",
    "STRONG": "#2166AC",
    "COMPONENT_STRONG": "#67A9CF",
    "MIXED": "#FDAE61",
    "WEAK_OR_ABSENT": "#E0E0E0",
    "NOT_PRESENT": "#E0E0E0",
}


def _source_map(tree, rows: Iterable[Mapping[str, object]]) -> dict[str, dict[str, object]]:
    materialized = [dict(row) for row in rows]
    labels = [str(row.get("tree_label", "")) for row in materialized]
    leaves = [leaf.name for leaf in tree.get_terminals()]
    if (
        any(not label for label in labels)
        or len(labels) != len(set(labels))
        or len(leaves) != len(set(leaves))
        or set(labels) != set(leaves)
    ):
        raise ValueError("figure source labels do not match unique tree labels")
    return {str(row["tree_label"]): row for row in materialized}


def _coordinates(tree):
    depths = tree.depths()
    if not depths or max(depths.values()) <= 0:
        depths = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    y_values = {
        leaf: float(len(terminals) - index) for index, leaf in enumerate(terminals)
    }

    def place(clade):
        if clade in y_values:
            return y_values[clade]
        child_y = [place(child) for child in clade.clades]
        y_values[clade] = (min(child_y) + max(child_y)) / 2.0
        return y_values[clade]

    place(tree.root)
    return depths, y_values


def _nice_scale(span: float) -> float:
    if span <= 0:
        return 0.1
    raw = span / 5.0
    exponent = math.floor(math.log10(raw))
    fraction = raw / 10**exponent
    base = 1 if fraction < 1.5 else 2 if fraction < 3.5 else 5
    return base * 10**exponent


def render_phylogram(
    tree_path: str | Path,
    source_rows: Iterable[Mapping[str, object]],
    output_prefix: str | Path,
    title: str,
) -> list[Path]:
    """Render one branch-length phylogram as editable SVG, PDF and 600-DPI PNG."""

    tree = Phylo.read(str(tree_path), "newick")
    source = _source_map(tree, source_rows)
    depths, y_values = _coordinates(tree)
    terminals = tree.get_terminals()
    max_x = max(depths[leaf] for leaf in terminals)
    span = max(max_x, 1e-6)
    marker_x = max_x + span * 0.025
    label_x = max_x + span * 0.055
    height = max(8.0, len(terminals) * 0.36 + 2.8)
    figure, axis = plt.subplots(figsize=(15.5, height))

    for parent in tree.find_clades(order="preorder"):
        if not parent.clades:
            continue
        children_y = [y_values[child] for child in parent.clades]
        axis.plot(
            [depths[parent], depths[parent]],
            [min(children_y), max(children_y)],
            color=BRANCH_COLOR,
            linewidth=0.85,
        )
        for child in parent.clades:
            axis.plot(
                [depths[parent], depths[child]],
                [y_values[child], y_values[child]],
                color=BRANCH_COLOR,
                linewidth=0.85,
            )
            if child.clades and child.name and "/" in str(child.name):
                axis.text(
                    (depths[parent] + depths[child]) / 2,
                    y_values[child] + 0.16,
                    str(child.name),
                    fontsize=5.8,
                    color="#444444",
                    ha="center",
                    va="bottom",
                )

    for leaf in terminals:
        row = source[str(leaf.name)]
        maize = str(row.get("species", "")) == "Zea mays"
        species_color = MAIZE_COLOR if maize else SORGHUM_COLOR
        axis.scatter(
            [depths[leaf]],
            [y_values[leaf]],
            marker="o" if maize else "s",
            s=22,
            color=species_color,
            zorder=3,
        )
        classification = str(row.get("classification", "CONTEXT_ONLY"))
        axis.scatter(
            [marker_x],
            [y_values[leaf]],
            marker="|",
            s=28,
            linewidths=2.3,
            color=CLASS_COLORS.get(classification, CLASS_COLORS["CONTEXT_ONLY"]),
            zorder=3,
        )
        label = str(row.get("display_label") or leaf.name)
        axis.text(
            label_x,
            y_values[leaf],
            label,
            fontsize=7.2,
            color=species_color,
            ha="left",
            va="center",
        )

    scale = _nice_scale(span)
    scale_y = 0.10
    axis.plot([0, scale], [scale_y, scale_y], color="black", linewidth=1.25)
    axis.plot([0, 0], [scale_y - 0.07, scale_y + 0.07], color="black", linewidth=1)
    axis.plot(
        [scale, scale], [scale_y - 0.07, scale_y + 0.07], color="black", linewidth=1
    )
    axis.text(
        scale / 2,
        scale_y - 0.18,
        f"{scale:g} substitutions/site",
        fontsize=7,
        ha="center",
        va="top",
    )

    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", color=MAIZE_COLOR,
            label="Zea mays (B73 RefGen_v5)", markersize=6,
        ),
        Line2D(
            [0], [0], marker="s", linestyle="none", color=SORGHUM_COLOR,
            label="Sorghum bicolor (BTx623 NCBIv3)", markersize=6,
        ),
        Line2D([0], [0], color="#444444", label="Internal labels: SH-aLRT / UFBoot (%)"),
    ]
    present = {str(row.get("classification", "CONTEXT_ONLY")) for row in source.values()}
    legend_handles.extend(
        Line2D(
            [0], [0], marker="|", linestyle="none", color=CLASS_COLORS[name],
            markeredgewidth=2.2, label=name.replace("_", " ").title(), markersize=8,
        )
        for name in CLASS_COLORS
        if name in present
    )
    axis.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.012),
        frameon=False,
        fontsize=6.4,
        ncol=3,
    )
    axis.set_title(title, loc="left", fontsize=12.2, fontweight="bold", pad=14)
    axis.text(
        0,
        -0.05,
        "UNROOTED maximum-likelihood topology; rectangular orientation is for display only. "
        "All 23 maize and 17 Sorghum source-locked representatives are retained.",
        transform=axis.transAxes,
        fontsize=6.8,
        color="#444444",
        ha="left",
        va="top",
    )
    axis.set_xlim(-span * 0.01, max_x + span * 0.78)
    axis.set_ylim(-0.45, len(terminals) + 1.0)
    axis.axis("off")
    figure.subplots_adjust(left=0.035, right=0.985, top=0.94, bottom=0.09)

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    artifacts = [
        prefix.with_suffix(".svg"),
        prefix.with_suffix(".pdf"),
        prefix.with_suffix(".png"),
    ]
    figure.savefig(artifacts[0], format="svg", bbox_inches="tight")
    figure.savefig(artifacts[1], format="pdf", bbox_inches="tight")
    figure.savefig(artifacts[2], format="png", dpi=600, bbox_inches="tight")
    plt.close(figure)
    return artifacts


def verify_figure_artifacts(
    tree_path: str | Path,
    source_rows: Iterable[Mapping[str, object]],
    artifacts: Iterable[str | Path],
) -> dict[str, object]:
    tree = Phylo.read(str(tree_path), "newick")
    source = _source_map(tree, source_rows)
    paths = [Path(path) for path in artifacts]
    if {path.suffix for path in paths} != {".svg", ".pdf", ".png"}:
        raise ValueError("figure artifacts must include SVG, PDF and PNG")
    if any(not path.is_file() or path.stat().st_size <= 1000 for path in paths):
        raise ValueError("figure artifact is absent or unexpectedly small")
    svg = next(path for path in paths if path.suffix == ".svg").read_text(
        encoding="utf-8"
    )
    all_labels = all(label in svg for label in source)
    if not all_labels:
        raise ValueError("not every tree label is traceable in SVG")
    support_labels = [
        str(clade.name)
        for clade in tree.get_nonterminals()
        if clade.name is not None and "/" in str(clade.name)
    ]
    if any(label not in svg for label in support_labels):
        raise ValueError("tree support label is absent from SVG")
    pdf = next(path for path in paths if path.suffix == ".pdf")
    if not pdf.read_bytes().startswith(b"%PDF"):
        raise ValueError("PDF signature is invalid")
    from PIL import Image

    png = next(path for path in paths if path.suffix == ".png")
    with Image.open(png) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0, 0))
    if width < 3000 or height < 2000:
        raise ValueError("PNG dimensions are below publication-resolution floor")
    return {
        "leaf_count": len(tree.get_terminals()),
        "all_labels_present_in_svg": all_labels,
        "support_labels_reconciled": len(support_labels),
        "png_width_px": width,
        "png_height_px": height,
        "png_reported_dpi": list(dpi),
        "artifacts": [str(path) for path in paths],
    }


def _matrix_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    materialized = [dict(row) for row in rows]
    aliases = [str(row.get("maize_family_alias", "")) for row in materialized]
    expected = {f"ZmNRG2_{index:03d}" for index in range(1, 24)}
    if len(materialized) != 23 or len(set(aliases)) != 23 or set(aliases) != expected:
        raise ValueError("robustness matrix requires exactly 23 unique maize gene aliases")
    return sorted(materialized, key=lambda row: str(row["maize_family_alias"]))


def render_robustness_matrix(
    rows: Iterable[Mapping[str, object]],
    output_prefix: str | Path,
    title: str,
) -> list[Path]:
    """Render the pre-specified 23-gene cross-model robustness matrix."""

    materialized = _matrix_rows(rows)
    candidates = (
        ("m09c_primary_state", "M09C\nprimary"),
        ("m09c_full_state", "M09C\nfull length"),
        ("r1_state", "R1\nC20"),
        ("r2_state", "R2\ndual mask"),
        ("r3_state", "R3\nfull mask"),
        ("r4_state", "R4\nDayhoff-6"),
    )
    columns = [(key, label) for key, label in candidates if key in materialized[0]]
    if len(columns) < 4:
        raise ValueError("robustness matrix is missing required analysis-state columns")

    figure, axis = plt.subplots(figsize=(13.6, 11.2))
    for row_index, row in enumerate(materialized):
        y = len(materialized) - 1 - row_index
        for column_index, (key, _label) in enumerate(columns):
            state = str(row.get(key, "NOT_PRESENT"))
            color = MATRIX_STATE_COLORS.get(state, "#FEE08B")
            rectangle = plt.Rectangle(
                (column_index - 0.47, y - 0.43),
                0.94,
                0.86,
                facecolor=color,
                edgecolor="white",
                linewidth=0.8,
            )
            axis.add_patch(rectangle)
            symbol = "E" if state in {"EXACT_STRONG", "STRONG"} else "C" if state == "COMPONENT_STRONG" else "·"
            axis.text(column_index, y, symbol, ha="center", va="center", fontsize=7.4, color="#202020")
        classification = str(row["final_computational_classification"])
        axis.scatter(
            [len(columns) + 0.13],
            [y],
            s=52,
            marker="s",
            color=CLASS_COLORS.get(classification, CLASS_COLORS["CONTEXT_ONLY"]),
        )
        candidate = str(row.get("sorghum_candidate_ids", "NOT_REPORTED"))
        axis.text(
            len(columns) + 0.42,
            y,
            candidate,
            ha="left",
            va="center",
            fontsize=6.5,
            color="#303030",
        )

    axis.set_xticks(range(len(columns)))
    axis.set_xticklabels([label for _key, label in columns], fontsize=8)
    axis.xaxis.tick_top()
    axis.set_yticks(range(len(materialized)))
    axis.set_yticklabels(
        [str(row["maize_family_alias"]) for row in reversed(materialized)],
        fontsize=7.6,
        color=MAIZE_COLOR,
    )
    axis.set_xlim(-0.55, len(columns) + 3.6)
    axis.set_ylim(-0.65, len(materialized) - 0.35)
    axis.tick_params(axis="both", length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_title(title, loc="left", fontsize=12.2, fontweight="bold", pad=34)
    legend_handles = [
        Line2D([0], [0], marker="s", linestyle="none", color=MATRIX_STATE_COLORS["EXACT_STRONG"], label="E: exact pair strongly supported"),
        Line2D([0], [0], marker="s", linestyle="none", color=MATRIX_STATE_COLORS["COMPONENT_STRONG"], label="C: collinearity component strongly supported"),
        Line2D([0], [0], marker="s", linestyle="none", color=MATRIX_STATE_COLORS["WEAK_OR_ABSENT"], label="·: weak or absent"),
    ]
    legend_handles.extend(
        Line2D([0], [0], marker="s", linestyle="none", color=CLASS_COLORS[name], label=name.replace("_", " ").title())
        for name in CLASS_COLORS
        if name != "CONTEXT_ONLY"
        and any(str(row["final_computational_classification"]) == name for row in materialized)
    )
    axis.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.025),
        frameon=False,
        fontsize=6.6,
        ncol=2,
    )
    axis.text(
        0,
        -0.115,
        "Strong support requires SH-aLRT ≥80 and UFBoot ≥95. Rows with multiple Sorghum partners are summarized conservatively; source TSV retains every pair.",
        transform=axis.transAxes,
        fontsize=6.7,
        color="#444444",
        ha="left",
        va="top",
    )
    figure.subplots_adjust(left=0.14, right=0.98, top=0.86, bottom=0.14)

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    artifacts = [prefix.with_suffix(".svg"), prefix.with_suffix(".pdf"), prefix.with_suffix(".png")]
    figure.savefig(artifacts[0], format="svg", bbox_inches="tight")
    figure.savefig(artifacts[1], format="pdf", bbox_inches="tight")
    figure.savefig(artifacts[2], format="png", dpi=600, bbox_inches="tight")
    plt.close(figure)
    return artifacts


def verify_matrix_artifacts(
    rows: Iterable[Mapping[str, object]], artifacts: Iterable[str | Path]
) -> dict[str, object]:
    materialized = _matrix_rows(rows)
    paths = [Path(path) for path in artifacts]
    if {path.suffix for path in paths} != {".svg", ".pdf", ".png"}:
        raise ValueError("matrix artifacts must include SVG, PDF and PNG")
    if any(not path.is_file() or path.stat().st_size <= 1000 for path in paths):
        raise ValueError("matrix artifact is absent or unexpectedly small")
    svg = next(path for path in paths if path.suffix == ".svg").read_text(encoding="utf-8")
    if any(str(row["maize_family_alias"]) not in svg for row in materialized):
        raise ValueError("not every maize label is traceable in matrix SVG")
    pdf = next(path for path in paths if path.suffix == ".pdf")
    if not pdf.read_bytes().startswith(b"%PDF"):
        raise ValueError("matrix PDF signature is invalid")
    from PIL import Image

    png = next(path for path in paths if path.suffix == ".png")
    with Image.open(png) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0, 0))
    if width < 3000 or height < 3000:
        raise ValueError("matrix PNG dimensions are below publication-resolution floor")
    return {
        "maize_gene_count": len(materialized),
        "all_labels_present_in_svg": True,
        "png_width_px": width,
        "png_height_px": height,
        "png_reported_dpi": list(dpi),
        "artifacts": [str(path) for path in paths],
    }
