# Personalized Fold-Change in gene expression

These codes are used in the paper:

M. Santolini, C. Rau, S. Ren, J. Saucerman, J. Wang, J. Weiss, Y. Wang, A. Lusis, A. Karma. A personalized, multi-omics approach identifies genes involved in cardiac hyper- trophy and heart failure. bioRxiv 120329; doi:https://doi.org/10.1101/120329, 2017.

## Genes pre-filtering
The first script, 1_compInitialFiltering.R, reproduces the filtering from Figure S1.

## Computation of FC genes
The second script, 2_compFCgenes.R, computes the gene selection process from Figure 2a to select 36 "FC" genes that best correlate with the individual hypertrophic response to a stressor.

## HMDP dataset
The file data.RData contains the expression fold-change matrix, the hypertrophic response and the list of gene symbols.

## Results
The folder output contains the expected output files.
