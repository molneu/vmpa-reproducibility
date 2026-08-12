# Supplementary Figure 9

Supplementary Figure 9 uses the temsirolimus RNA-seq and western-blot datasets from Figure 6.

`reproducibility/figure6/scripts/generate_temsirolimus_vmpa_scores.R` uses the gene-symbol annotation included with the raw-count table, sums duplicate-symbol raw counts, removes all-zero genes, applies TMM normalization and log-CPM transformation, and calculates glioma-context VMPA scores using `molneu/vmpaR` commit `5e1c248` with `n = 250`, `min_conf = 1`, and `unique = FALSE`.

`scripts/run_supplementary_figure9.R` generates the VMPA score files and then runs the western-blot correlation and plotting scripts stored under `reproducibility/figure6/scripts/wb_correlation`.

The annotated raw counts and western-blot source data are provided separately with the source-data archive.
