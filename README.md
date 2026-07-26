# VMPA Reproducibility

This repository contains the scripts needed to reproduce the analyses and figure panels from the VMPA manuscript.

Large input data are stored separately on Figshare:

https://figshare.com/s/93eb2c2824a3cefdcd92

After downloading the Figshare data, place the extracted `reproducibility/` folder contents into this repository so that each script can find its corresponding `data/` folder, for example:

```text
reproducibility/
  figure1/
    data/
    scripts/
  figure2/
    data/
    scripts/
  ...
```

Generated outputs are written to `results/` folders. These folders are intentionally not tracked in Git.

## Setup

Install R packages from the lockfile:

```r
install.packages("renv")
renv::restore()
```

Run scripts from the repository root. Examples:

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

## Citation

Please cite our VMPA manuscript and the Figshare dataset when reusing this code or data.
