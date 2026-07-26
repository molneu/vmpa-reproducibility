manifest_path <- file.path("metadata", "figshare_manifest.tsv")

if (!file.exists(manifest_path)) {
  stop("Cannot find ", manifest_path, ". Run this script from the repository root.")
}

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
manifest <- manifest[manifest$figshare_url != "TODO_FIGSHARE_FILE_URL", , drop = FALSE]

if (nrow(manifest) == 0) {
  stop("No Figshare URLs have been added to ", manifest_path, ".")
}

for (i in seq_len(nrow(manifest))) {
  destination <- manifest$relative_path[i]
  url <- manifest$figshare_url[i]

  if (grepl("/$", destination)) {
    message("Skipping directory placeholder: ", destination)
    next
  }

  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(destination)) {
    message("Already exists: ", destination)
    next
  }

  message("Downloading: ", destination)
  utils::download.file(url, destination, mode = "wb", quiet = FALSE)

  if (grepl("\\.zip$", destination, ignore.case = TRUE)) {
    message("Unzipping: ", destination)
    utils::unzip(destination, exdir = dirname(destination))
  }
}
