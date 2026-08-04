# Reproducibility repository for the VMPA manuscript

This repository contains the scripts needed to reproduce the main figures from the VMPA manuscript.

When executing scripts, these will produce figure outputs to `results/` folders.

## What You Need

- R/RStudio

## Step 1: Get This Repository

Option A, using Git:

Open a terminal and execute:

```bash
git clone https://github.com/molneu/vmpa-reproducibility.git
cd vmpa-reproducibility
```

Option B, without Git:

1. Open the GitHub repository in a browser.
2. Click `Code`.
3. Click `Download ZIP`.
4. Unzip the repository.
5. Open a terminal in the unzipped `vmpa-reproducibility` folder.

All commands below should be run from the folder that contains this `README.md` file (repository root)

## Step 2: Download the compressed Data from Figshare and extract

Download this file from Figshare:

https://figshare.com/s/f5cc8be3979e6b03b6a5

```text
vmpa_reproducibility_input_data.tar
```
Place the downloaded file into the repository root:

```text
vmpa-reproducibility/
  README.md
  renv.lock
  vmpa_reproducibility_input_data.tar
  reproducibility/
```

Then extract:

e.g.:

```bash
tar -xf vmpa_reproducibility_input_data.tar
```

## Step 3: Restore the R Environment

Open R or RStudio from the repository root and run in the R console:

```r
install.packages("renv")
renv::restore()
```

This installs the R package versions recorded in `renv.lock`.

## Step 4: Run a Figure Script

Run scripts from the repository root. For example, to reproduce Figure 1:

in Terminal: 

```bash
cd /path/to/vmpa-reproducibility
Rscript reproducibility/figure1/scripts/'Figure 1a-g.R'
```

in R/Rstudio

```r
setwd("/path/to/vmpa-reproducibility")
source("reproducibility/figure1/scripts/Figure 1a-g.R")
```

The output files will be written to the corresponding `results/` folder.

## Figure Scripts

Figure 1:

```bash
Rscript reproducibility/figure1/scripts/'Figure 1a-g.R'
```

Figure 2:

```bash
Rscript reproducibility/figure2/scripts/'Figure2 a-f.R'
```

Figure 3:

```bash
Rscript reproducibility/figure3/scripts/figure3b.R
Rscript reproducibility/figure3/scripts/figure3c.R
Rscript reproducibility/figure3/scripts/figure3d.R
Rscript reproducibility/figure3/scripts/figure3e.R
```

Figure 4:

```bash
Rscript reproducibility/figure4/scripts/figure4a-e.R
```

Figure 5:

```bash
Rscript reproducibility/figure5/scripts/Figure5b-g.R
```

Figure 6:

```bash
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
```

Supplementary Figure 5:

```bash
Rscript reproducibility/supplementary_figure5/scripts/matched_family_benchmark/paper_context_only_effect_threshold_figures.R
```

## Notes

Some scripts can take a long time to execute

## Citation

Cima I, et al. Modeling the active proteome by context-matched perturbation signatures. Manuscript in preparation.
