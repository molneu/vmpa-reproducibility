#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
rank_file <- file.path(out_dir, "vmpa_ranked_candidates_n250_sample_pop_sd_unique_TRUE_FALSE.csv")
if (!file.exists(rank_file)) stop("Missing VMPA ranking file: ", rank_file)

args <- commandArgs(trailingOnly = TRUE)
screen_rule <- if (length(args) > 0) tolower(args[1]) else "both"
if (!screen_rule %in% c("both", "either", "positive")) screen_rule <- "both"
file_tag <- if (screen_rule == "either") "_top5_either" else if (screen_rule == "positive") "_positive_3of4" else ""
top5_column <- if (screen_rule == "either") "n_lines_top5_either_dose" else "n_lines_top5_both_doses"
plot_sort_column <- if (screen_rule == "positive") "n_lines_up_both_doses" else top5_column
display_n <- if (screen_rule == "positive") 10 else 12
display_tag <- if (screen_rule == "positive") "_top10" else ""

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(httr2)
  library(jsonlite)
  library(ggplot2)
})

ranked <- read.csv(rank_file, check.names = FALSE, stringsAsFactors = FALSE)

secondary_rule <- if (screen_rule == "positive") {
  rep(TRUE, nrow(ranked))
} else {
  ranked[[top5_column]] >= 1
}

# These rules are fixed before inspecting individual genes. They define an
# unbiased data-first screen followed by an independent targetability filter.
data_first <- ranked[
  ranked$n_lines_up_both_doses >= 3 &
    secondary_rule,
  ,
  drop = FALSE
]
data_first$targetability_gene <- sub("\\.[0-9]+$", "", toupper(data_first$gene))

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(data_first$targetability_gene),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
)
symbol_map <- symbol_map[
  !is.na(symbol_map$ENSEMBL) & !duplicated(symbol_map$SYMBOL),
  ,
  drop = FALSE
]
data_first$ensembl_id <- symbol_map$ENSEMBL[
  match(data_first$targetability_gene, symbol_map$SYMBOL)
]

query_target <- function(ensembl_id) {
  if (is.na(ensembl_id) || !nzchar(ensembl_id)) {
    return(data.frame(
      ensembl_id = ensembl_id,
      approved_symbol = NA_character_,
      api_status = "unmapped",
      targetability_tier = 1,
      targetability_basis = "unmapped",
      phase1_or_advanced_clinical = FALSE,
      phase1_clinical_flag = FALSE,
      advanced_clinical_flag = FALSE,
      approved_drug_flag = FALSE,
      clinical_tractability_status = "No clinical-phase evidence",
      small_molecule_direct_flag = FALSE,
      antibody_access_flag = FALSE,
      antibody_localization_flag = FALSE,
      protac_support_flag = FALSE,
      literature_flag = FALSE,
      ubiquitination_flag = FALSE,
      half_life_flag = FALSE,
      small_molecule_binder_flag = FALSE,
      pocket_family_flag = FALSE,
      direct_tractability_flag = FALSE,
      supporting_tractability_flag = FALSE,
      positive_tractability_evidence = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  query <- paste0(
    "query($ensemblId: String!){ target(ensemblId: $ensemblId){ ",
    "approvedSymbol tractability { label modality value } } }"
  )

  for (attempt in seq_len(3)) {
    ans <- tryCatch({
      response <- request("https://api.platform.opentargets.org/api/v4/graphql") |>
        req_body_json(list(query = query, variables = list(ensemblId = ensembl_id))) |>
        req_timeout(30) |>
        req_perform()
      resp_body_json(response, simplifyVector = TRUE)
    }, error = function(e) NULL)

    if (!is.null(ans) && !is.null(ans$data) && !is.null(ans$data$target)) {
      tr <- ans$data$target$tractability
      if (is.null(tr) || length(tr) == 0) {
        tr <- data.frame(label = character(), modality = character(), value = logical())
      } else {
        tr <- as.data.frame(tr, stringsAsFactors = FALSE)
      }
      value <- if (nrow(tr)) as.logical(tr$value) else logical()
      label <- if (nrow(tr)) as.character(tr$label) else character()
      modality <- if (nrow(tr)) as.character(tr$modality) else character()

      clinical_labels <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
      small_molecule_labels <- c(
        "Structure with Ligand", "High-Quality Ligand", "Small Molecule Binder",
        "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family"
      )
      antibody_localization_labels <- c(
        "UniProt loc high conf", "GO CC high conf", "UniProt loc med conf",
        "UniProt SigP or TMHMM", "GO CC med conf", "Human Protein Atlas loc"
      )

      approved <- any(value & label == "Approved Drug")
      phase1 <- any(value & label == "Phase 1 Clinical")
      advanced <- any(value & label == "Advanced Clinical")
      clinical <- any(value & label %in% clinical_labels)
      small_molecule_direct <- any(value & modality == "SM" & label %in% small_molecule_labels)
      pocket_or_family <- any(value & modality == "SM" & label %in% c(
        "High-Quality Pocket", "Med-Quality Pocket", "Druggable Family"
      ))
      antibody_access <- any(value & modality == "AB" & label == "UniProt SigP or TMHMM")
      antibody_localization <- any(value & modality == "AB" & label %in% antibody_localization_labels)
      literature <- any(value & modality == "PR" & label == "Literature")
      ubiquitination <- any(value & modality == "PR" & label %in% c(
        "UniProt Ubiquitination", "Database Ubiquitination"
      ))
      half_life <- any(value & modality == "PR" & label == "Half-life Data")
      small_molecule_binder <- any(value & modality == "PR" & label == "Small Molecule Binder")
      protac_support <- any(literature, ubiquitination, half_life, small_molecule_binder)
      direct_tractability <- any(clinical, small_molecule_direct, antibody_access)
      supporting_tractability <- any(protac_support, pocket_or_family, antibody_localization)

      tier <- if (clinical) 4 else if (direct_tractability || protac_support) 3 else if (supporting_tractability) 2 else 1
      clinical_status <- if (approved) {
        "Approved drug"
      } else if (advanced) {
        "Advanced clinical"
      } else if (phase1) {
        "Phase I clinical"
      } else {
        "No clinical-phase evidence"
      }
      basis_bits <- character()
      if (approved) basis_bits <- c(basis_bits, "approved drug")
      if (clinical && !approved) basis_bits <- c(basis_bits, "advanced clinical or phase-I")
      if (small_molecule_direct) basis_bits <- c(basis_bits, "small-molecule evidence")
      if (antibody_access) basis_bits <- c(basis_bits, "antibody-accessible")
      if (protac_support) basis_bits <- c(basis_bits, "PROTAC-supporting evidence")
      if (pocket_or_family) basis_bits <- c(basis_bits, "pocket or druggable family")
      if (antibody_localization && !antibody_access) basis_bits <- c(basis_bits, "antibody localization evidence")
      basis <- if (length(basis_bits)) paste(basis_bits, collapse = "; ") else "no positive tractability flag"
      positive_evidence <- unique(paste(modality[value], label[value], sep = ":"))

      return(data.frame(
        ensembl_id = ensembl_id,
        approved_symbol = if (is.null(ans$data$target$approvedSymbol)) NA_character_ else ans$data$target$approvedSymbol,
        api_status = "ok",
        targetability_tier = tier,
        targetability_basis = basis,
        phase1_or_advanced_clinical = clinical,
        phase1_clinical_flag = phase1,
        advanced_clinical_flag = advanced,
        approved_drug_flag = approved,
        clinical_tractability_status = clinical_status,
        small_molecule_direct_flag = small_molecule_direct,
        antibody_access_flag = antibody_access,
        antibody_localization_flag = antibody_localization,
        protac_support_flag = protac_support,
        literature_flag = literature,
        ubiquitination_flag = ubiquitination,
        half_life_flag = half_life,
        small_molecule_binder_flag = small_molecule_binder,
        pocket_family_flag = pocket_or_family,
        direct_tractability_flag = direct_tractability,
        supporting_tractability_flag = supporting_tractability,
        positive_tractability_evidence = paste(positive_evidence, collapse = "; "),
        stringsAsFactors = FALSE
      ))
    }
    Sys.sleep(attempt)
  }

  data.frame(
    ensembl_id = ensembl_id,
    approved_symbol = NA_character_,
    api_status = "failed",
    targetability_tier = 1,
    targetability_basis = "API query failed",
    phase1_or_advanced_clinical = FALSE,
    phase1_clinical_flag = FALSE,
    advanced_clinical_flag = FALSE,
    approved_drug_flag = FALSE,
    clinical_tractability_status = "No clinical-phase evidence",
    small_molecule_direct_flag = FALSE,
    antibody_access_flag = FALSE,
    antibody_localization_flag = FALSE,
    protac_support_flag = FALSE,
    literature_flag = FALSE,
    ubiquitination_flag = FALSE,
    half_life_flag = FALSE,
    small_molecule_binder_flag = FALSE,
    pocket_family_flag = FALSE,
    direct_tractability_flag = FALSE,
    supporting_tractability_flag = FALSE,
    positive_tractability_evidence = NA_character_,
    stringsAsFactors = FALSE
  )
}

ids <- unique(na.omit(data_first$ensembl_id))
ot_file <- file.path(out_dir, paste0("vmpa_initial_open_targets_tractability", file_tag, ".csv"))
if (file.exists(ot_file)) {
  message("Reading cached Open Targets results: ", ot_file)
  ot <- read.csv(ot_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  message("Querying Open Targets for ", length(ids), " unique genes...")
  ot <- lapply(ids, query_target)
  ot <- do.call(rbind, ot)
  write.csv(ot, ot_file, row.names = FALSE)
}

screen <- merge(data_first, ot, by = "ensembl_id", all.x = TRUE, sort = FALSE)
screen$targetability_tier[is.na(screen$targetability_tier)] <- 1
screen$targetability_basis[is.na(screen$targetability_basis)] <- "unmapped/no positive flag"
screen$phase1_clinical_flag[is.na(screen$phase1_clinical_flag)] <- FALSE
screen$advanced_clinical_flag[is.na(screen$advanced_clinical_flag)] <- FALSE
screen$approved_drug_flag[is.na(screen$approved_drug_flag)] <- FALSE
screen$clinical_tractability_status[is.na(screen$clinical_tractability_status)] <- "No clinical-phase evidence"
screen$targetability_percentile <- (screen$targetability_tier - 1) / 3
screen$data_percentile <- 1 - ((screen$rank - 1) / pmax(max(screen$rank) - 1, 1))
screen$combined_screen_score <- 0.5 * screen$data_percentile +
  0.5 * screen$targetability_percentile
screen$passes_targetability <- screen$targetability_tier >= 3
screen <- screen[order(
  -screen$passes_targetability,
  -screen$targetability_tier,
  -screen[[plot_sort_column]],
  -screen$n_lines_up_both_doses,
  -screen$median_mean_delta,
  screen$rank
), , drop = FALSE]
screen$screen_rank <- seq_len(nrow(screen))
screen$plot_label <- paste0(screen$target_id, " [", screen$signature_mode, "]")
screen$highlight_bcl2l1 <- screen$gene == "BCL2L1"
screen$clinical_plot_group <- ifelse(
  screen$clinical_tractability_status %in% c(
    "Phase I clinical", "Advanced clinical", "Approved drug"
  ),
  "Clinical evidence (Phase I+)",
  "Direct/PROTAC/antibody tractability"
)

write.csv(
  screen,
  file.path(out_dir, paste0("vmpa_initial_unbiased_targetable_screen_n250_sample_pop_sd", file_tag, ".csv")),
  row.names = FALSE
)

targetable <- screen[screen$passes_targetability, , drop = FALSE]
plot_data <- do.call(rbind, lapply(c(FALSE, TRUE), function(mode) {
  z <- targetable[targetable$unique_signatures == mode, , drop = FALSE]
  head(z[order(
    -z[[plot_sort_column]],
    -z$n_lines_up_both_doses,
    -z$median_mean_delta
  ), , drop = FALSE], display_n)
}))

write.csv(
  plot_data,
  file.path(out_dir, paste0("vmpa_initial_unbiased_targetable_plot_data_n250_sample_pop_sd", file_tag, display_tag, ".csv")),
  row.names = FALSE
)

plot <- ggplot(
  plot_data,
  aes(
    x = median_mean_delta,
    y = reorder(plot_label, median_mean_delta),
    color = clinical_plot_group,
    size = .data[[plot_sort_column]],
    shape = highlight_bcl2l1
  )
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ signature_mode, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Direct/PROTAC/antibody tractability" = "#6baed6",
      "Clinical evidence (Phase I+)" = "#b2182b"
    ),
    name = "Targetability"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Other candidate", "TRUE" = "BCL2L1"),
    name = NULL
  ) +
  scale_size_continuous(range = c(2.5, 6), breaks = 1:4) +
  labs(
    x = "Median VMPA activity delta: mean Temsi 50/500 nM vs DMSO",
    y = NULL,
    size = if (screen_rule == "positive") {
      "Cell lines positive at both doses"
    } else {
      paste0(
        "Cell lines top 5% at ",
        ifelse(screen_rule == "either", "either dose", "both doses")
      )
    },
    title = "Candidates emerging from an unbiased VMPA screen",
    subtitle = paste0(
      if (screen_rule == "positive") {
        "Positive at both doses in >=3/4 lines; top-5% status reported descriptively"
      } else {
        paste0(
          "Positive at both doses in >=3/4 lines; top 5% at ",
          ifelse(screen_rule == "either", "either dose", "both doses")
        )
      },
      "; independent Open Targets targetability annotation"
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, paste0("vmpa_initial_unbiased_targetable_dotplot_n250_sample_pop_sd", file_tag, display_tag, ".pdf")),
  plot, width = 12, height = 8
)
ggsave(
  file.path(out_dir, paste0("vmpa_initial_unbiased_targetable_dotplot_n250_sample_pop_sd", file_tag, display_tag, ".png")),
  plot, width = 12, height = 8, dpi = 220
)

cat(
  "Wrote unbiased VMPA/targetability screen (top-5% rule: ",
  screen_rule,
  ") to ",
  normalizePath(out_dir),
  "\n",
  sep = ""
)
