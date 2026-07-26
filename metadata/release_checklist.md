# Release Checklist

## GitHub

1. Create the private GitHub repository `molneu/vmpa-reproducibility`.
2. Initialize this local folder as a Git repository if needed.
3. Commit the code, documentation, `renv.lock`, and final manuscript figure previews.
4. Push to the private GitHub repository.
5. Confirm that no `data/`, `results/`, `.DS_Store`, archive folders, or local editor files are tracked.

## Figshare

1. Upload the input files listed in `metadata/figshare_manifest.tsv`.
2. Use `metadata/figshare_upload_plan.tsv` to locate the local source file for each manifest row.
3. Replace every `TODO_FIGSHARE_FILE_URL` value with the direct Figshare download URL for that file.
4. Reuse the same Figshare file URL for duplicate large files if the same input is needed at multiple relative paths.
5. Keep the Figshare entry private until the manifest is complete.
6. Run a fresh-clone test:

```bash
Rscript scripts/download_figshare_data.R
Rscript reproducibility/figure6/scripts/wb_correlation/make_four_mtor_validation_scatterplots.R
Rscript reproducibility/figure6/scripts/target_selection/vmpa_three_dotplots_unique_FALSE_n250_sample_pop_sd.R positive 10
Rscript reproducibility/figure6/scripts/apoptosis_pathway/make_compass_intrinsic_apoptosis_pathway.R
Rscript reproducibility/figure6/scripts/caspase_dapi_response/generate_main_2x2_figures.R
Rscript reproducibility/supplementary_figure4/scripts/matched_family_benchmark/paper_context_only_effect_threshold_figures.R
```

## Current Figshare Blockers

The current working folder has most manifest files locally, but these paths still need to be resolved before publication:

- `reproducibility/figure2/data/Human.GRCh38.p13.annot.tsv.gz`
- `reproducibility/figure3/data/subsets/`
- `reproducibility/figure3/data/level5_beta_trt_xpr_n142901x12328.gctx`
- `reproducibility/figure3/data/siginfo_beta.txt`
- `reproducibility/figure3/data/LINC gene annotations.csv`
- `reproducibility/figure3/data/250717_original_annotations_Dressler_role.csv`
- `reproducibility/figure3/data/250717_original_annotations_Dressler_status.csv`
- `reproducibility/figure3/data/250717_original_annotations_Kinnersley_role.csv`
- `reproducibility/figure3/data/Human.GRCh38.p13.annot.tsv.gz`
- `reproducibility/figure3/data/geo_raw_counts.zip`
- `reproducibility/figure4/data/Human.GRCh38.p13.annot.tsv.gz`
- `reproducibility/figure4/data/GSE111571_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE169418_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE171163_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE86518_raw_counts_GRCh38.p13_NCBI.tsv.gz`

Some of these are duplicate inputs already present under another figure folder or available elsewhere locally. Resolve them by either uploading once and reusing the same Figshare URL in the manifest, or by removing the duplicate manifest row if the corresponding script no longer needs it.
