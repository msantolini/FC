# Personalized Fold-Change in gene expression

These codes are used in the paper:

Santolini et al., "A personalized, multi-omics approach identifies genes involved in cardiac hypertrophy and heart failure".

The first code, 1_compInitialFiltering.R, reproduces the filtering from Figure S1.
The second code, 2_compFCgenes.R, computes the gene selection process from Figure 2a to select 36 "pFC" genes that best correlate with the individual hypertrophic response to a stressor.

The file data.RData contains the expression fold-change matrix, the hypertrophic response and the list of gene symbols.
