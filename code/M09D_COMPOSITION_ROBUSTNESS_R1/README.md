# M09D Composition-Aware Robustness R1

This versioned implementation reuses the immutable M09C 40-protein evidence
and executes four pre-registered, unrooted IQ-TREE 3 robustness analyses. It
does not rerun MAFFT, sequence extraction, DIAMOND, HMMER or MCScanX.

## End-to-end execution

```bash
bash MASTER_RUN_M09D.sh /home/hammadse/AYESHA_ZmNRL1_WORKFLOW
```

The runner prepares an exclusive staging directory, executes R1–R4
sequentially with at most two threads, and builds the author-review package.
It deliberately stops before publication so the five figure sets can be
visually inspected.

## Resume

```bash
bash MASTER_RUN_M09D.sh /home/hammadse/AYESHA_ZmNRL1_WORKFLOW \
  --resume /absolute/path/to/.STAGING_<token>
```

Resume is permitted only when the frozen command, input, executable and core
implementation provenance match. Completed IQ-TREE artifacts are reused.

## Publication gate

After programmatic checks and independent visual review, place a PASS
`PHYLOGENY/FIGURES/M09D_VISUAL_QA_ATTESTATION.json` in staging and run:

```bash
/home/hammadse/miniconda3/envs/wrky_q1/bin/python m09d_finalize.py \
  /home/hammadse/AYESHA_ZmNRL1_WORKFLOW publish --staging /absolute/path/to/.STAGING_<token>
```

Publication reruns source checks, producer-independent validation, complete
manifest verification and the protected manuscript hash gate before an
atomic staging-to-release rename.

## Scientific interpretation

Strong support is SH-aLRT at least 80 and UFBoot at least 95. All trees are
unrooted. Component support is not exact one-to-one orthology, and no output
constitutes experimental or functional proof.
