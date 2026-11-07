args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) {
  stop("Usage: Rscript .github/scripts/bundle-test-fixtures.R /external/output/directory [test-fixtures-vN]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
output <- args[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
if (output == root || startsWith(output, paste0(root, "/"))) {
  stop("The release bundle must be written outside the package repository.")
}
tag <- if (length(args) == 2L) args[[2L]] else "test-fixtures-v1"
if (!grepl("^test-fixtures-v[0-9]+$", tag)) {
  stop("Use an explicit fixture tag such as test-fixtures-v2.")
}
asset <- paste0("SCWorkflow-", tag, ".tar.gz")
archive <- file.path(output, asset)
manifest <- file.path(output, "fixture-bundle.json")
checksums <- file.path(output, "SHA256SUMS")
if (any(file.exists(c(archive, manifest, checksums)))) {
  stop("Refusing to overwrite an existing bundle; choose a new output directory.")
}
fixtures <- file.path(root, "tests", "testthat", "fixtures")
files <- sort(list.files(fixtures, recursive = TRUE, all.files = TRUE,
                         no.. = TRUE, full.names = FALSE))
paths <- file.path(fixtures, files)
if (!length(files) || any(dir.exists(paths)) || any(nzchar(Sys.readlink(paths)))) {
  stop("Expected a nonempty fixture directory containing regular files, not links.")
}
records <- lapply(seq_along(files), function(index) {
  list(path = files[[index]], bytes = unname(file.info(paths[[index]])$size),
       sha256 = digest::digest(file = paths[[index]], algo = "sha256"))
})
previous <- setwd(dirname(fixtures))
tryCatch(
  utils::tar(archive, files = file.path("fixtures", files),
             compression = "gzip", tar = "internal"),
  finally = setwd(previous)
)
lock <- list(schema = 1L, repo = "NIDAP-Community/SCWorkflow", tag = tag,
             asset = asset, bytes = unname(file.info(archive)$size),
             sha256 = digest::digest(file = archive, algo = "sha256"),
             files = records)
jsonlite::write_json(lock, manifest, pretty = TRUE, auto_unbox = TRUE)
writeLines(paste(lock$sha256, asset), checksums)
message("Bundled ", length(files), " files (",
        round(sum(vapply(records, `[[`, numeric(1L), "bytes")) / 1024^2, 2),
        " MiB uncompressed) into ", archive)