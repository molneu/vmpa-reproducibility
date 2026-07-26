# Figshare Entry

Private share URL: https://figshare.com/s/93eb2c2824a3cefdcd92

Title: VMPA manuscript reproducibility data

Type: Dataset

License: CC BY 4.0

Categories:

- Genomics and transcriptomics
- Translational and applied bioinformatics

Keywords:

- VMPA
- protein activity
- transcriptomics
- single-cell RNA-seq
- perturbation signatures
- reproducibility

Description:

This dataset contains the input data required to reproduce the analyses and figure panels accompanying the VMPA manuscript. The associated GitHub repository contains the executable analysis scripts, documentation, and package lockfile. Large input data are archived here rather than in GitHub; generated `results/` folders and output figures are intentionally excluded from the repository and are recreated by running the scripts.

The file organization follows the GitHub repository layout, with inputs placed under `reproducibility/figure*/data/` and `reproducibility/supplementary_figure*/data/`. The manifest in the GitHub repository maps each Figshare file URL to the relative path expected by the scripts.

Status notes:

- Keep the entry private until all manifest rows have direct Figshare file URLs and a fresh-clone test succeeds.
- Add the public DOI after publication/finalization.
- Add direct file download URLs to `metadata/figshare_manifest.tsv` after all files are uploaded.
