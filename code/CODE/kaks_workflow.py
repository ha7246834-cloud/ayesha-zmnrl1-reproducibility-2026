"""Sequential, source-locked post-M09D maize duplicate-pair analysis.

Two explicit stages: align, then estimate after an alignment-review attestation.
No existing release or canonical manuscript is changed.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Phylo.PAML import yn00
import kaks_core as core

BIN = Path('/home/hammadse/miniconda3/envs/wrky_q1/bin')
QA = '11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/QA/'
PAIR_FILE = QA + 'M07B_INTRA_GENOME/RELEASE_20261031T112429Z_219b181d/M07B_FAMILY_UNIQUE_COLLINEAR_PAIRS.tsv'
MAP_FILE = QA + 'M09C_INDEPENDENT_ORTHOLOGY_PHYLOGENY/RELEASE_20261107T133337Z_f5edc166/PHYLOGENY/INPUTS/M09C_SEQUENCE_MAP.tsv'
CDS_FILE = 'outputs/downloads/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cds.fa.gz'
PROTEIN_FILE = 'outputs/downloads/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.gz'
LOCKS = {
    PAIR_FILE: '94ff86d9bcd648393894603cfe51e6cfd4a8c7f798de7044f74760246c9e7e26',
    MAP_FILE: 'd058707aa3c746656371135ab086d726783fa739c530707461216bfba21e8f14',
    CDS_FILE: 'd1a10cf63a626979c1373e5eed2ea504672ff6d8d643c8fa84d5a85eb1958d2c',
    PROTEIN_FILE: '0fc902a74b79af6c7f6d110b9843cf636a30392537852aad6ac31ce3457e82f2',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    if not rows:
        raise ValueError('Cannot silently write an empty evidence table')
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def validate_pairs(pairs, mappings):
    identity = {r['tree_label']: r['gene_id'] for r in mappings}
    seen = set()
    for pair in pairs:
        a, b = pair['alias_a'], pair['alias_b']
        key = frozenset((a, b))
        if len(key) != 2 or key in seen:
            raise ValueError('Self-pair or duplicate unordered pair')
        seen.add(key)
        if identity.get(a) != pair['gene_a'] or identity.get(b) != pair['gene_b']:
            raise ValueError('Pair gene/alias mapping does not match retained representative')
    return pairs


def alignment_gate(metrics):
    # Operational review triage, not a proof of homology or a universal cutoff.
    if (metrics['paired_codons'] < 300 or metrics['retained_fraction_shorter'] < .7
            or metrics['protein_identity_fraction'] < .3):
        return 'REVIEW_ALIGNMENT_INSUFFICIENT_FOR_AUTOMATIC_ELIGIBILITY'
    return 'PASS_TECHNICAL_SCREEN_PENDING_VISUAL_REVIEW'


def command(directory, label, argv, output):
    record = directory / (label + '.command.json')
    done = directory / (label + '.done.json')
    if done.exists():
        saved = json.loads(done.read_text())
        if saved['argv'] != argv or any(sha(directory / p) != digest for p, digest in saved['outputs'].items()):
            raise ValueError(f'Completed {label} command/output changed')
        return
    write_json(record, {'argv': argv, 'cwd': str(directory), 'started_utc': datetime.now(timezone.utc).isoformat()})
    with (directory / output).open('w') as stdout, (directory / (label + '.stderr.txt')).open('w') as stderr:
        result = subprocess.run(argv, cwd=directory, stdout=stdout, stderr=stderr, env={**os.environ, 'OMP_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '1'})
    if result.returncode:
        raise RuntimeError(f'{label} returned {result.returncode}; see stderr')
    outputs = {output: sha(directory / output)}
    if label == 'yn00':
        outputs['yn00.out'] = sha(directory / 'yn00.out')
    write_json(done, {'argv': argv, 'outputs': outputs})


def selected_fasta(path, identifiers):
    records = {}
    with gzip.open(path, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in identifiers:
                if record.id in records:
                    raise ValueError('Duplicate source FASTA identifier')
                records[record.id] = str(record.seq).upper()
    if set(records) != set(identifiers):
        raise ValueError(f'Missing source sequences: {set(identifiers) - set(records)}')
    return records


def fasta_write(path, records):
    Path(path).write_text(''.join(f'>{key}\n{seq}\n' for key, seq in records.items()))


def run(root, m09d, output, phase):
    root, m09d, output = (Path(p).resolve() for p in (root, m09d, output))
    core.require_m09d_completion(m09d)  # Before any downstream write or calculation.
    if not output.is_relative_to(root / QA / 'POST_SYNTENY_KAKS'):
        raise ValueError('Output must be in dedicated POST_SYNTENY_KAKS candidate directory')
    for name, expected in LOCKS.items():
        if sha(root / name) != expected:
            raise ValueError(f'Source-lock mismatch: {name}')
    mappings = read_tsv(root / MAP_FILE)
    pairs = validate_pairs(read_tsv(root / PAIR_FILE), mappings)
    if len(pairs) != 8:
        raise ValueError('Expected exactly eight source-validated maize pairs')
    by_alias = {r['tree_label']: r for r in mappings}
    aliases = {p[k] for p in pairs for k in ('alias_a', 'alias_b')}
    selected = {a: by_alias[a] for a in aliases}
    proteins = selected_fasta(root / PROTEIN_FILE, {r['protein_id'] for r in selected.values()})
    coding = selected_fasta(root / CDS_FILE, {r['protein_id'].replace('_P', '_T') for r in selected.values()})
    sequences, sequence_qc = {}, []
    for alias, row in sorted(selected.items()):
        protein = proteins[row['protein_id']].removesuffix('*')
        if hashlib.sha256(protein.encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError(f'{alias}: current protein differs from M09C retained isoform')
        transcript = row['protein_id'].replace('_P', '_T')
        cds = core.validated_cds(coding[transcript], protein)
        sequences[alias] = (protein, cds)
        sequence_qc.append({'alias': alias, 'gene_id': row['gene_id'], 'protein_id': row['protein_id'], 'transcript_id': transcript,
                            'protein_length': len(protein), 'provider_cds_length': len(coding[transcript]), 'stop_free_cds_length': len(cds),
                            'protein_sha256': row['sequence_sha256'], 'stop_free_cds_sha256': hashlib.sha256(cds.encode()).hexdigest(), 'translation_status': 'PASS_EXACT'})
    caveats = {r['maize_family_alias']: r for r in read_tsv(m09d / 'M09D_ALL_23_FINAL_REVIEW.tsv')}
    current_provenance = {'source_sha256': LOCKS, 'm09d_manifest_sha256': sha(m09d / 'SHA256SUMS.txt'),
                          'executables': {name: sha(BIN / name) for name in ('mafft', 'pal2nal.pl', 'yn00')},
                          'code': {p.name: sha(p) for p in Path(__file__).parent.glob('*.py')}}
    provenance = output / 'INPUT_AND_CODE_LOCK.json'
    if output.exists():
        if not provenance.is_file() or json.loads(provenance.read_text()) != current_provenance:
            raise ValueError('Candidate provenance differs; do not overwrite')
    else:
        if phase != 'align':
            raise ValueError('Alignment phase must precede estimation')
        output.mkdir(parents=True)
        write_json(provenance, current_provenance)
    write_tsv(output / 'MATCHED_ISOFORM_TRANSLATION_QC.tsv', sequence_qc)
    write_tsv(output / 'M07B_EIGHT_VALIDATED_PAIRS.tsv', pairs)
    results = []
    for pair in pairs:
        a, b = pair['alias_a'], pair['alias_b']
        pair_id = f'{a}__{b}'
        directory = output / 'PAIRS' / pair_id
        directory.mkdir(parents=True, exist_ok=True)
        fasta_write(directory / 'proteins.faa', {'seqA': sequences[a][0], 'seqB': sequences[b][0]})
        fasta_write(directory / 'coding.fna', {'seqA': sequences[a][1], 'seqB': sequences[b][1]})
        command(directory, 'mafft', [str(BIN / 'mafft'), '--localpair', '--maxiterate', '1000', '--thread', '2', 'proteins.faa'], 'protein_alignment.faa')
        alignment = {r.id: str(r.seq) for r in SeqIO.parse(directory / 'protein_alignment.faa', 'fasta')}
        metrics = core.alignment_metrics(alignment)
        expected, columns = core.independent_paired_codons(alignment, {'seqA': sequences[a][1], 'seqB': sequences[b][1]})
        command(directory, 'pal2nal', [str(BIN / 'pal2nal.pl'), 'protein_alignment.faa', 'coding.fna', '-output', 'paml', '-nogap', '-codontable', '1'], 'codon_alignment.phy')
        observed = {r.id: str(r.seq) for r in AlignIO.read(directory / 'codon_alignment.phy', 'phylip-sequential')}
        if observed != expected:
            raise ValueError(f'{pair_id}: independent backtranslation differs from PAL2NAL')
        write_tsv(directory / 'RETAINED_PROTEIN_COLUMNS.tsv', [{'retained_codon': i, 'protein_alignment_column_1based': column} for i, column in enumerate(columns, 1)])
        inherited = ';'.join(f'{alias}:{caveats[alias]["final_computational_classification"]}:{caveats[alias]["composition_and_mask_caveats"]}' for alias in (a, b))
        row = {'pair_id': pair_id, **pair, **metrics, 'alignment_screen': alignment_gate(metrics),
               'inherited_m09d_context': inherited, 'estimation_status': 'NOT_RUN_PENDING_ALIGNMENT_REVIEW',
               'Ka': 'NA', 'Ks': 'NA', 'Ka_SE': 'NA', 'Ks_SE': 'NA', 'Ka_Ks': 'NA', 'rate_flags': 'NOT_ESTIMATED', 'scientific_status': 'REVIEW', 'MYA': 'NA'}
        if phase == 'estimate':
            review_path = output / 'ALIGNMENT_REVIEW.json'
            if not review_path.is_file():
                raise ValueError('Independent alignment-review attestation required')
            review = json.loads(review_path.read_text()).get(pair_id, {})
            if review.get('protein_alignment_sha256') != sha(directory / 'protein_alignment.faa'):
                raise ValueError(f'{pair_id}: review is missing or stale')
            if review.get('decision') == 'ELIGIBLE_EXPLORATORY':
                if not row['alignment_screen'].startswith('PASS'):
                    raise ValueError(f'{pair_id}: automatic screen failed; revise justification before estimation')
                (directory / 'yn00.ctl').write_text('seqfile = codon_alignment.phy\noutfile = yn00.out\nverbose = 1\nicode = 0\nweighting = 1\ncommonf3x4 = 0\n')
                command(directory, 'yn00', [str(BIN / 'yn00'), 'yn00.ctl'], 'yn00.stdout.txt')
                try:
                    estimates = yn00.read(str(directory / 'yn00.out'))['seqA']['seqB']['YN00']
                    screen = core.rate_screen(estimates['dN'], estimates['dS'], estimates['dS SE'])
                    row.update(Ka=screen['dN'] if screen['dN'] is not None else 'NA', Ks=screen['dS'] if screen['dS'] is not None else 'NA', Ka_SE=estimates['dN SE'], Ks_SE=estimates['dS SE'],
                               Ka_Ks=screen['omega'] if screen['omega'] is not None else 'NA', rate_flags=';'.join(screen['flags']), estimation_status='COMPUTED_EXPLORATORY_REVIEW')
                except (ValueError, KeyError, IndexError) as error:
                    row.update(estimation_status='REVIEW_UNPARSABLE_ESTIMATE', rate_flags=str(error))
            else:
                row.update(estimation_status='NOT_ESTIMATED_ALIGNMENT_REVIEW', rate_flags=review.get('reason', 'No eligible decision'))
        results.append(row)
        write_tsv(output / 'POST_SYNTENY_KAKS_MASTER.tsv', results)
    write_json(output / 'POST_SYNTENY_KAKS_QC.json', {'status': 'REVIEW', 'phase': phase, 'pairs': len(results),
               'exact_translations': len(sequence_qc), 'focal_gene_pair_count': 0, 'm09d_release': str(m09d),
               'interpretation': 'Descriptive maize collinear-pair context only. No functional orthology, formal selection or divergence-time inference. Saturation not excluded.',
               'completed_utc': datetime.now(timezone.utc).isoformat()})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project_root', type=Path)
    parser.add_argument('m09d_release', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('phase', choices=('align', 'estimate'))
    args = parser.parse_args()
    run(args.project_root, args.m09d_release, args.output, args.phase)


if __name__ == '__main__':
    main()
