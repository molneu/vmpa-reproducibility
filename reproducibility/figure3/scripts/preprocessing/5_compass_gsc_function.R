#' Extract gene signatures for a given organ-of-origin context
#' 
#' @param context       character; preferred context, one of: 
#'                      "glioma","melanoma","nsclc","gastric","ovarian",
#'                      "crc","breast","prostate","pdac","headneck"
#' @param subset_dir    directory where the per-context RDS files live
#'                      (default: ".../figure3/data/subsets")
#' @param n             integer; how many bottom genes to return per signature
#'                      (default: 200)
#' @param min_conf      numeric; minimum cps_conf_total to include a signature
#'                      (default: 1)
#' @param targets       optional character vector of protein activity targets
#'                      (default: NULL means keep all)
#' @param driver_filter logical; if TRUE, only keep signatures with 
#'                      cancer_driver_summary != "None" (default: FALSE)
#' @param output        one of "list", "df", or "gsc" (default: "list")
#' @return if output="list", a named list of character vectors;
#'         if output="df", a data.frame with one column per signature;
#'         if output="gsc", a GeneSetCollection of GeneSet objects
compass_gsc <- function(context,
                        subset_dir   = here::here("reproducibility", "figure3", "data", "subsets"),
                        n            = 200,
                        min_conf     = 1,
                        targets      = NULL,
                        driver_filter = FALSE,
                        output       = c("list", "df", "gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # 1) Load the context‐specific subsetted GCT object
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) {
    stop("No subset file found for context: ", context)
  }
  gct <- readRDS(subset_file)
  
  # 2) Select high‐confidence signatures by total CPS score
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  
  # 2b) (optional) filter by cancer‐driver annotation
  if (driver_filter) {
    cds <- gct@cdesc$cancer_driver_summary
    keep_idx <- keep_idx[cds[keep_idx] != "None"]
  }
  
  if (length(keep_idx) == 0) {
    if (output == "gsc") return(GeneSetCollection(list()))
    if (output == "df")  return(data.frame())
    return(list())
  }
  
  # 3) Extract matrix & reassign rownames to gene symbols
  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  
  # 4) Build output names: "<id>_c<TOTAL_CONF>"
  ids      <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # 5) (Optional) filter to only specified targets
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[keep_idx] %in% targets
    if (!any(keep_tgt)) {
      warning("No signatures found for targets: ", paste(targets, collapse = ", "))
      if (output == "gsc") return(GeneSetCollection(list()))
      if (output == "df")  return(data.frame())
      return(list())
    }
    mat       <- mat[, keep_tgt,       drop = FALSE]
    out_names <- out_names[keep_tgt]
  }
  
  # 6) Extract bottom N genes per signature, remove "MIA2" duplicate
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    ord  <- suppressWarnings(order(vals, decreasing = FALSE, na.last = "keep"))
    gs   <- head(rownames(mat)[ord], n)
    # Explicitly drop "MIA2" duplicate
    gs   <- gs[gs != "MIA2"]        # drop all MIA2 rows
  })
  names(res_list) <- out_names
  
  # 7) Return as GeneSetCollection
  if (output == "gsc") {
    gs <- lapply(out_names, function(setname) {
      GeneSet(geneIds = res_list[[setname]], setName = setname)
    })
    return(GeneSetCollection(gs))
  }
  
  # 8) Return as data.frame
  if (output == "df") {
    max_len <- max(lengths(res_list))
    df <- as.data.frame(
      do.call(cbind, lapply(res_list, function(vec) {
        length(vec) <- max_len; vec
      })),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(df)
  }
  
  # 9) Default: return as list
  res_list
}
