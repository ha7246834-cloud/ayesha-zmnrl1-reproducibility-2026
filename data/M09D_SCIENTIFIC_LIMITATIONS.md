# M09D Scientific Limitations

M09D status: **REVIEW**. This stage provides composition-aware computational robustness evidence only. It does not establish experimental orthology, biochemical function, functional conservation, or stress-response equivalence.

## Model and composition limitations

- The unchanged dual-domain analysis failed the original IQ-TREE composition chi-square test for 4 taxa: Sb_GeneID_8077197, ZmNRG2_011, ZmNRG2_012, ZmNRG2_013.
- The unchanged full-length analysis failed that test for 18 taxa: Sb_GeneID_8055907, Sb_GeneID_8066701, Sb_GeneID_8069565, Sb_GeneID_8076242, Sb_GeneID_8076865, Sb_GeneID_8077197, Sb_GeneID_8077213, Sb_GeneID_8077392, ZmNRG2_001, ZmNRG2_002, ZmNRG2_006, ZmNRG2_007, ZmNRG2_009, ZmNRG2_011, ZmNRG2_012, ZmNRG2_013, ZmNRG2_016, ZmNRG2_019.
- IQ-TREE symmetry summaries yielded 3 BH-significant alignment/partition-by-test-family results. Individual taxon-pair p-values were not exported by IQ-TREE 3.1.2, so pair-specific symmetry failure cannot be assigned from these outputs.
- R1 uses a bounded C20 profile mixture. It addresses across-site heterogeneity but is not a full nonstationary across-lineage model; its near-zero mixture weights triggered an overfitting warning.
- Dayhoff-6 recoding reduces compositional detail and is a sensitivity analysis, not a preferred biological model.
- AIC/BIC comparisons are valid only for fits to the identical character matrix. Masked and recoded analyses are explicitly non-comparable by information criteria.

## Alignment and sampling limitations

- 32 alignment rows, representing 31 unique taxa, exceed 50% gaps or ambiguity. The row/taxon distinction is retained because Sb_GeneID_8076242 is flagged in both matrices.
- Five source-locked domain-completeness concerns are retained across 5 taxa: Sb_GeneID_8058473, Sb_GeneID_8076242, Sb_GeneID_8077392, ZmNRG2_012, ZmNRG2_021. No sequence was removed or padded.
- The analysis contains exactly 23 maize and 17 Sorghum representatives and no defensible outgroup; every tree is therefore unrooted. Direction, ancestral state, and duplication timing are not inferred.
- Component-level support can reflect duplicated or paralog-associated collinearity and is not equivalent to an exact one-to-one relationship.

## Classification limitations

- Reciprocal rank-1 homology, domain architecture, collinearity and phylogenetic support are jointly interpreted. Bootstrap support alone is never treated as proof of orthology.
- Seven NO_RBH cases are ineligible for either supported-candidate class by the pre-registered rules.
- The release remains AUTHOR_REVIEW. No content was promoted into `MASTER_CURRENT.docx`.

## Reasons for REVIEW

- 4 maize outcomes are CONTRADICTED.
- 3 maize outcomes retain direct composition/masking caveats.
- 3 alignment-by-test-family symmetry summaries have BH q<0.05.
- the R1 profile-mixture fit reported near-zero weights and an overfitting warning.
- IQ-TREE did not export individual taxon-pair symmetry p-values for pair-level caveat assignment.
