# Release Checklist

## GitHub

1. Create the private GitHub repository `molneu/vmpa-reproducibility`.
2. Initialize this local folder as a Git repository if needed.
3. Commit the code, documentation, `renv.lock`, and final manuscript figure previews.
4. Push to the private GitHub repository.
5. Confirm that no `data/`, `results/`, `.DS_Store`, archive folders, or local editor files are tracked.

## Figshare

1. Review `metadata/figshare_ftp_upload_status.tsv`.
2. Replace every `TODO_FIGSHARE_FILE_URL` value in `metadata/figshare_manifest.tsv` with the direct Figshare download URL for that file.
3. For rows marked `reuse_existing_upload_url`, reuse the URL from `reuse_source_relative_path`; do not upload duplicate files.
4. Resolve the rows marked `unresolved_not_uploaded`.
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

The FTP upload to the private Figshare share has been run without duplicate shared files. These paths still need to be resolved before publication:

- `reproducibility/figure3/data/geo_raw_counts.zip`
- `reproducibility/figure4/data/GSE111571_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE169418_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE171163_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `reproducibility/figure4/data/GSE86518_raw_counts_GRCh38.p13_NCBI.tsv.gz`

Shared inputs such as the LINCS GCTX, `siginfo_beta.txt`, `LINC gene annotations.csv`, and `Human.GRCh38.p13.annot.tsv.gz` were uploaded once only. Reuse the same Figshare URL for duplicate manifest rows.
