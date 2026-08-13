# Reproducibility repository for the VMPA manuscript

This repository contains the scripts needed to reproduce the main figures from the VMPA manuscript.

Scripts write figures and tables to the corresponding `results/` folders.

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

All commands below should be run from the folder containing this `README.md` file (the repository root).

## Step 2: Download and extract the compressed data from Figshare

Download this file from Figshare:

https://figshare.com/s/f5cc8be3979e6b03b6a5

```text
vmpa_reproducibility_input_data_v1.0.tar
```
Place the downloaded file into the repository root:

```text
vmpa-reproducibility/
  README.md
  renv.lock
  vmpa_reproducibility_input_data_v1.0.tar
  reproducibility/
```

Then extract:

e.g.:

```bash
tar -xf vmpa_reproducibility_input_data_v1.0.tar
```

## Step 3: Restore the R Environment

Open R or RStudio from the repository root and run in the R console:

```r
install.packages("renv")
renv::restore()
```

This installs the R package versions recorded in `renv.lock`, including the revision of `vmpaR` used for the Figure 6 analyses.

## Step 4: Run a Figure Script

Run scripts from the repository root. For example, to reproduce Figure 1:

in Terminal: 

```bash
cd /path/to/vmpa-reproducibility
Rscript reproducibility/figure1/scripts/'Figure 1a-g.R'
```

or in R/RStudio:

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
Rscript reproducibility/figure6/scripts/'figure 6a-d.R'
Rscript reproducibility/figure6/scripts/'Fig6 dose response curves_final.R'
Rscript reproducibility/figure6/scripts/run_figure6h_i_j.R \
  figure6_temsirolimus/Figure6_temsirolimus_raw_counts.csv \
  figure6_temsirolimus/'260725 WB_corrected_data_with Ponceau.csv' \
  .
Rscript reproducibility/figure6/scripts/caspase_dapi_response/generate_main_2x2_figures.R
```

The `run_figure6h_i_j.R` wrapper generates the VMPA scores and runs the western-blot correlation, target-selection, Figure 6i plotting, and apoptosis-pathway scripts in the required order.

Supplementary Figure 2:

```bash
cd reproducibility/supplementary_figure2/scripts/foundation_models
bash run_supplementary_figure2_foundation_models.sh \
  ../../../../supplementary_figure2_foundation_models/data_export \
  /path/to/BulkFormer \
  local_outputs
```

See `reproducibility/supplementary_figure2/scripts/foundation_models/README.md` for the required Python environments and model inputs.

Supplementary Figure 3:

```bash
Rscript reproducibility/supplementary_figure3/scripts/supplementary_figure3_overlap_rate.R
```

Supplementary Figure 5:

```bash
Rscript reproducibility/supplementary_figure5/paper_context_only_effect_threshold_figures.R
```

Supplementary Figure 6:

```bash
Rscript reproducibility/supplementary_figure6/scripts/supplementary_figure6_function_annotation.R
```

Supplementary Figure 7 (run after Figure 5):

```bash
Rscript reproducibility/supplementary_figure7/scripts/supplementary_figure7_random_forest_validation.R
```

Supplementary Figure 9:

```bash
Rscript reproducibility/supplementary_figure9/scripts/run_supplementary_figure9.R \
  figure6_temsirolimus/Figure6_temsirolimus_raw_counts.csv \
  figure6_temsirolimus/'260725 WB_corrected_data_with Ponceau.csv' \
  .
```

## Notes

Some scripts can take a long time to execute.

## Citation

Cima I, et al. Modeling the active proteome by context-matched perturbation signatures. Manuscript in preparation.
