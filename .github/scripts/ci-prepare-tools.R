arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript .github/scripts/ci-prepare-tools.R OUTPUT_DIRECTORY", call. = FALSE)
}
if (as.character(getRversion()) != "4.1.3") {
  stop("Prepare CI tools only in the pinned R 4.1.3 environment.", call. = FALSE)
}

output <- normalizePath(arguments[[1]], mustWork = FALSE)
library <- file.path(output, "library")
base_libraries <- setdiff(.libPaths(), library)
base_packages <- installed.packages(lib.loc = base_libraries)[,
  c("Package", "Version", "LibPath"), drop = FALSE
]
expected <- c(Seurat = "4.1.1", SeuratObject = "4.1.1", Matrix = "1.5-1", ggplot2 = "3.3.6")
actual <- base_packages[match(names(expected), base_packages[, "Package"]), "Version"]
if (anyNA(actual) || !identical(unname(actual), unname(expected))) {
  stop("The base image does not contain the verified legacy package versions.", call. = FALSE)
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
sources <- jsonlite::fromJSON(file.path(dirname(script), "ci-tool-sources.json"))

dir.create(library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library, base_libraries))

for (index in seq_len(nrow(sources))) {
  package <- sources$Package[[index]]
  version <- sources$Version[[index]]
  if (package %in% base_packages[, "Package"]) {
    message("Reuse existing base package: ", package)
    next
  }
  installed <- installed.packages(lib.loc = library)
  if (!package %in% installed[, "Package"]) {
    archive <- file.path(output, "sources", paste0(package, "_", version, ".tar.gz"))
    if (!file.exists(archive)) stop("Missing pinned source archive: ", package, call. = FALSE)
    if (digest::digest(file = archive, algo = "sha256") != sources$SHA256[[index]]) {
      stop("Source checksum mismatch: ", package, call. = FALSE)
    }
    utils::install.packages(archive, lib = library, repos = NULL, type = "source",
                           dependencies = FALSE)
  }
  if (!requireNamespace(package, quietly = TRUE) ||
      as.character(utils::packageVersion(package)) != version) {
    stop("Pinned CI package cannot load at the expected version: ", package, call. = FALSE)
  }
}

after <- installed.packages(lib.loc = base_libraries)[,
  c("Package", "Version", "LibPath"), drop = FALSE
]
additions <- installed.packages(lib.loc = library)
if (!identical(base_packages, after) ||
    any(additions[, "Package"] %in% base_packages[, "Package"])) {
  stop("CI tool preparation must not replace any base package.", call. = FALSE)
}
utils::write.csv(as.data.frame(base_packages), file.path(output, "base-packages.csv"),
                 row.names = FALSE)
utils::write.csv(sources, file.path(output, "tools-manifest.csv"), row.names = FALSE)
previous <- setwd(output)
tryCatch(utils::tar("library.tar.gz", files = "library", compression = "gzip"),
         finally = setwd(previous))
message("Pinned CI tools prepared; all base packages and versions are unchanged.")