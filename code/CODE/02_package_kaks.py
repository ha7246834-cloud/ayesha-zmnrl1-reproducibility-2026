"""Package completed estimates without rerunning alignment or estimation."""
import shutil
import sys
from pathlib import Path
import Bio
from Bio.Phylo.PAML import yn00
from package_core import read_tsv, write_tsv, write_json, sha, workbook, manifest, classify_rates

root, source, dest = map(Path, sys.argv[1:])
if dest.exists():
    raise SystemExit('Refusing existing destination')
rows = read_tsv(source / 'POST_SYNTENY_KAKS_MASTER.tsv')
assert len(rows) == 8
checks = []
for r in rows:
    raw = source / 'PAIRS' / r['pair_id'] / 'yn00.out'
    if r['Ka'] == 'NA':
        assert not raw.exists()
        assert r['alias_a'] == 'ZmNRG2_011' and r['alias_b'] == 'ZmNRG2_013'
        continue
    v = yn00.read(str(raw))['seqA']['seqB']['YN00']
    assert float(r['Ka']) == v['dN'] and float(r['Ks']) == v['dS']
    assert abs(float(r['Ka_Ks']) - v['dN']/v['dS']) < 1e-12
    checks.append({'pair_id':r['pair_id'], 'raw_yn00_sha256':sha(raw),
                   'reported_YN00_omega':v['omega'], 'Ka_Ks_recomputed_from_printed_rates':r['Ka_Ks'],
                   'status':'PASS'})
summary = classify_rates(rows)
assert summary['computed'] == 7 and summary['withheld'] == 1
shutil.copytree(source, dest)
write_tsv(dest / 'INDEPENDENT_RATE_RECONCILIATION.tsv', checks)
workbook(dest / 'POST_SYNTENY_KAKS_AUTHOR_TABLES.xlsx', {
    'Eight_pair_outcomes': rows,
    'Isoform_translation_QC': read_tsv(source / 'MATCHED_ISOFORM_TRANSLATION_QC.tsv'),
    'Rate_reconciliation': checks})
code = root / '11_MANUSCRIPT/UBUNTU_MASTER_TOOLS/M09D_POST_SYNTENY_KAKS_R1'
shutil.copytree(code, dest / 'CODE', ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(__file__, dest / 'CODE/02_package_kaks.py')
shutil.copy2(Path(__file__).with_name('package_core.py'), dest / 'CODE/package_core.py')
write_json(dest / 'PACKAGE_QC.json', {**summary, 'computational_integrity':'PASS',
    'source_candidate':str(source), 'raw_estimates_independently_reconciled':7,
    'workbook_cell_readback':'PASS', 'exact_translation_unique_proteins':13,
    'focal_pair_count':0, 'PAML':'4.10.10 (29 Jan 2026)', 'PAL2NAL':'14.1',
    'Biopython':Bio.__version__, 'MYA':'NA; no clock used',
    'limitations':['Saturation not independently excluded',
                   'Pairwise estimates are not formal selection tests',
                   'Divergent C-terminal alignment precludes ZmNRG2_011/013 estimation',
                   'No eligible pair contains focal ZmNRL1',
                   'Printed PAML rates have four decimal places; derived ratios are not extra precision']})
(dest / 'METHODS_RESULTS_LIMITATIONS.md').write_text('''# Post-synteny Ka/Ks — author-review evidence

Eight unique intra-maize family-family anchors from M07B were screened after
M09D completion. The 13 endpoint representatives were selected by the M09C
identity map and matched to same-isoform official B73 NAM 5.0 CDS records.
All 13 translations matched exactly. Sequential MAFFT L-INS-i searches used
--localpair --maxiterate 1000 --thread 2. PAL2NAL 14.1 used -output paml -nogap
-codontable 1; an independent codon reconstruction matched every retained
protein-alignment column. Controls and commands are retained in PAIRS/.

Prespecified technical triage required >=300 paired codons, >=70% retention
of the shorter protein and >=30% protein identity. This was not a substitute
for alignment inspection. Hash-bound inspection excluded ZmNRG2_011/013:
its divergent, indel-rich C-terminal region lacked convincing positional
correspondence. No local trimming was used to force eligibility.

Seven pairs were estimated using Yang–Nielsen in PAML yn00 4.10.10 (Yang and
Nielsen 2000; Yang 2007), with icode=0, weighting=1 and commonf3x4=0. Ka and Ks
are per nonsynonymous and synonymous site, respectively. Ka ranged from
0.0238 to 0.2575 and Ks from 0.1289 to 1.2656. Ratios recomputed from the
printed rates ranged from 0.0776 to 0.3871; raw program omega is separately
preserved, and final displays avoid spurious extra precision. All seven
ratios were below one, but this descriptive result is not a formal
selection test. No pair includes focal ZmNRL1.

Nonfinite, negative or sentinel rates are unavailable; zero Ks produces
NA rather than infinity. Ks >=2 was prespecified as a warning, not proof
that smaller estimates are unsaturated. None reached that cutoff and no
zero Ks occurred. Saturation has not independently been excluded, notably
for the two Ks estimates above one, so all estimates remain REVIEW.
No divergence time, lineage clock, MYA or functional inference was made.

Sources remain unchanged. INPUT_AND_CODE_LOCK.json pins original data,
programs and workflow. ALIGNMENT_REVIEW.json records actual pre-estimation
inspection. INDEPENDENT_RATE_RECONCILIATION.tsv traces every reported value
to raw yn00 output. Existing .done checkpoints remain reusable at their
original candidate paths. No inference was rerun for this package.
''')
(dest / 'REPRODUCIBILITY.md').write_text(f'''# Reproducibility

Original run directory: `{source}`

Each PAIRS/<pair>/ contains the exact command JSON, alignments, PAL2NAL
output, program stdout/stderr, yn00.ctl and raw results. The original
checkpoint remains intact; no restart is required. To verify/reuse it:

```bash
/home/hammadse/miniconda3/envs/wrky_q1/bin/python \\
  11_MANUSCRIPT/UBUNTU_MASTER_TOOLS/M09D_POST_SYNTENY_KAKS_R1/01_run_post_synteny_kaks.py \\
  {root} \\
  {root}/11_MANUSCRIPT/UBUNTU_SINGLE_MASTER/QA/M09D_COMPOSITION_ROBUSTNESS/RELEASE_20261108T013057Z_a4ab944f \\
  {source} estimate
```

This command validates locks and completed-run hashes and reuses existing
estimates. Do not use it to rerun inference unnecessarily. A changed input,
model or workflow requires a new version, not alteration of this release.
''')
print({'release':str(dest), 'manifest_entries':manifest(dest), **summary})
