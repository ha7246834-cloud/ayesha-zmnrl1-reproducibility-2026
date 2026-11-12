"""Conservative post-synteny coding-sequence gates; no biological inference."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

from Bio.Seq import Seq


def validated_cds(cds: str, protein: str) -> str:
    """Return stop-free CDS only on exact canonical protein translation."""
    cds = cds.upper()
    if not cds or len(cds) % 3 or set(cds) - set("ACGT"):
        raise ValueError("CDS frame/alphabet failure")
    translated = str(Seq(cds).translate())
    if translated.endswith("*"):
        translated, cds = translated[:-1], cds[:-3]
    if "*" in translated or translated != protein:
        raise ValueError("CDS does not exactly translate to retained protein")
    return cds


def alignment_metrics(alignment: dict[str, str]) -> dict:
    if len(alignment) != 2 or len({len(s) for s in alignment.values()}) != 1:
        raise ValueError("Exactly two equal-length aligned proteins required")
    a, b = alignment.values()
    alphabet = set("ACDEFGHIKLMNPQRSTVWY-")
    if set(a + b) - alphabet:
        raise ValueError("Noncanonical protein symbols")
    lengths = [len(s.replace("-", "")) for s in (a, b)]
    paired = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not paired or not min(lengths):
        raise ValueError("No paired canonical residues")
    return {
        "alignment_columns": len(a),
        "paired_codons": len(paired),
        "removed_gap_columns": len(a) - len(paired),
        "protein_identity_fraction": sum(x == y for x, y in paired) / len(paired),
        "retained_fraction_shorter": len(paired) / min(lengths),
        "retained_fraction_longer": len(paired) / max(lengths),
        "maximum_gap_run": max([len(m.group()) for s in (a, b) for m in re.finditer("-+", s)] or [0]),
    }


def independent_paired_codons(alignment: dict[str, str], cds: dict[str, str]) -> tuple[dict[str, str], list[int]]:
    """Independently reconstruct the gap-free PAL2NAL result for comparison."""
    alignment_metrics(alignment)
    if set(alignment) != set(cds):
        raise ValueError("CDS/alignment identifiers differ")
    expanded = {}
    for name, protein in alignment.items():
        sequence = validated_cds(cds[name], protein.replace("-", ""))
        index = 0
        codons = []
        for amino_acid in protein:
            if amino_acid == "-":
                codons.append("---")
            else:
                codons.append(sequence[index:index + 3])
                index += 3
        expanded[name] = codons
    keys = list(alignment)
    sites = [i for i, (a, b) in enumerate(zip(*expanded.values())) if "-" not in a + b]
    return {key: "".join(expanded[key][i] for i in sites) for key in keys}, [i + 1 for i in sites]


def rate_screen(dn: float | None, ds: float | None, ds_se: float | None) -> dict:
    """Screen identifiability without claiming a formal saturation test.

    Ks>=2 is a conservative warning label, not a universal saturation boundary.
    Large sentinel values are unavailable, never evidence of positive selection.
    """
    flags = ["SATURATION_NOT_EXCLUDED"]
    valid = all(v is not None and math.isfinite(v) and 0 <= v < 99 for v in (dn, ds))
    omega = None
    if not valid:
        flags.append("INVALID_OR_SENTINEL_RATE")
    elif ds == 0:
        flags.append("ZERO_KS")
    else:
        omega = dn / ds
        if ds >= 2:
            flags.append("HIGH_KS")
        if ds_se is None or not math.isfinite(ds_se) or ds_se < 0:
            flags.append("KS_UNCERTAINTY_UNAVAILABLE")
        elif ds - 1.96 * ds_se <= 0:
            flags.append("KS_SE_INTERVAL_INCLUDES_ZERO")
    return {"dN": dn if valid else None, "dS": ds if valid else None,
            "omega": omega, "flags": flags, "scientific_status": "REVIEW"}


def require_m09d_completion(release: Path) -> dict:
    """Fail closed on incomplete trees, adjudication, or release hashes."""
    manifest = release / "SHA256SUMS.txt"
    qc_path = release / "M09D_QC.json"
    if not manifest.is_file() or not qc_path.is_file():
        raise ValueError("M09D release completion and manifest unavailable")
    paths = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split(None, 1)
        name = name.lstrip(" *")
        path = (release / name).resolve()
        if not path.is_relative_to(release.resolve()) or not path.is_file():
            raise ValueError("Unsafe or missing M09D manifest entry")
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"M09D source changed: {name}")
        paths.add(name)
    qc = json.loads(qc_path.read_text())
    if qc.get("status") not in {"PASS", "REVIEW"}:
        raise ValueError("M09D scientific FAIL or absent status")
    for analysis in ("R1", "R2", "R3", "R4"):
        if f"CHECKPOINTS/{analysis}_COMPLETE.json" not in paths or analysis not in qc.get("analyses", {}):
            raise ValueError(f"M09D {analysis} incomplete")
    matrix = release / "M09D_ALL_23_FINAL_REVIEW.tsv"
    if matrix.name not in paths:
        raise ValueError("M09D gene adjudication absent from manifest")
    with matrix.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if {r["maize_family_alias"] for r in rows} != {f"ZmNRG2_{i:03d}" for i in range(1, 24)} or len(rows) != 23:
        raise ValueError("M09D gene adjudication incomplete")
    return qc
