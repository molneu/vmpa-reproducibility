# VMPA Reproducibility

This repository contains the scripts and lightweight documentation needed to reproduce the analyses and figures from the VMPA manuscript.

Large input data are not stored in GitHub. They should be archived on Figshare and downloaded into the expected `reproducibility/figure*/data/` folders before running the scripts. The `results/` folders are created and populated only when the scripts are run.

## Repository Contents

- `reproducibility/figure1/` to `reproducibility/figure6/`: figure-specific analysis scripts and final figure previews.
- `metadata/figshare_manifest.tsv`: list of data files that should be uploaded to Figshare.
- `scripts/download_figshare_data.R`: helper script for downloading Figshare-hosted files once the Figshare record is available.
- `renv.lock`: R package lockfile for restoring the analysis environment.

## Data Availability

The Figshare entry is currently available as a private share:

https://figshare.com/s/93eb2c2824a3cefdcd92

After the Figshare files are finalized:

1. Upload the files listed in `metadata/figshare_manifest.tsv`.
2. Replace `TODO_FIGSHARE_FILE_URL` values in the manifest with direct Figshare download URLs.
3. Add the Figshare DOI here and in the manuscript data availability statement.

Figshare DOI: `TODO_FIGSHARE_DOI`

## Setup

Install R packages from the lockfile:

```r
install.packages("renv")
renv::restore()
```

Download data after the Figshare URLs are added:

```bash
Rscript scripts/download_figshare_data.R
```

Run a figure script from the repository root, for example:

```bash
Rscript reproducibility/figure1/scripts/'Figure 1a-g.R'
Rscript reproducibility/figure6/scripts/'Fig6a,c,d.R'
Rscript reproducibility/figure6/scripts/'Fig6 dose response curves_final.R'
Rscript reproducibility/figure6/scripts/wb_correlation/analyze_original_wb_mtor_s6_across_cell_lines.R
Rscript reproducibility/figure6/scripts/wb_correlation/make_four_mtor_validation_scatterplots.R
Rscript reproducibility/figure6/scripts/target_selection/prepare_activity_matrix_from_raw_result.R
Rscript reproducibility/figure6/scripts/target_selection/make_vmpa_treatment_response_tables.R
Rscript reproducibility/figure6/scripts/target_selection/make_vmpa_average_temsi_response.R
Rscript reproducibility/figure6/scripts/target_selection/rank_vmpa_candidates_n250_sample_pop_sd.R
Rscript reproducibility/figure6/scripts/target_selection/vmpa_initial_unbiased_targetable_screen_n250_sample_pop_sd.R positive
Rscript reproducibility/figure6/scripts/target_selection/vmpa_three_dotplots_unique_FALSE_n250_sample_pop_sd.R positive 10
Rscript reproducibility/figure6/scripts/apoptosis_pathway/make_compass_intrinsic_apoptosis_pathway.R
Rscript reproducibility/figure6/scripts/caspase_dapi_response/generate_main_2x2_figures.R
Rscript reproducibility/supplementary_figure4/scripts/matched_family_benchmark/paper_context_only_effect_threshold_figures.R
```

## Suggested Citation

Please cite the VMPA manuscript and the Figshare dataset when reusing this code or data.
