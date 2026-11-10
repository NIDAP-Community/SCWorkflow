scworkflow_fixture_manifest <- function(root) {
  if (!is.character(root) || length(root) != 1L || is.na(root)) {
    stop("root must be a single directory path.")
  }
  lock <- jsonlite::read_json(file.path(root, "tests", "fixture-bundle.json"),
                            simplifyVector = TRUE)
  required <- c("schema", "repo", "tag", "asset", "bytes", "sha256", "files")
  if (!all(required %in% names(lock))) {
    stop("Fixture manifest is missing required fields.")
  }
  safe_path <- function(paths) {
    is.character(paths) && !anyNA(paths) && all(nzchar(paths)) &&
      !any(grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.{1,2}(/|$)|[[:cntrl:]])", paths))
  }
  if (!identical(lock$schema, 1L) || length(lock$repo) != 1L ||
      !grepl("^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$", lock$repo) ||
      length(lock$tag) != 1L || !safe_path(lock$tag) ||
      length(lock$asset) != 1L || !safe_path(lock$asset) ||
      basename(lock$asset) != lock$asset ||
      length(lock$sha256) != 1L || !grepl("^[a-f0-9]{64}$", lock$sha256) ||
      !is.numeric(lock$bytes) || length(lock$bytes) != 1L ||
      !is.finite(lock$bytes) || lock$bytes <= 0 ||
      !is.data.frame(lock$files) || !nrow(lock$files) ||
      !all(c("path", "bytes", "sha256") %in% names(lock$files))) {
    stop("Invalid fixture manifest.")
  }
  records <- lock$files
  if (!safe_path(records$path) || anyDuplicated(records$path) ||
      !is.numeric(records$bytes) || any(!is.finite(records$bytes)) ||
      any(records$bytes < 0) || !is.character(records$sha256) ||
      anyNA(records$sha256) || any(!grepl("^[a-f0-9]{64}$", records$sha256))) {
    stop("Invalid fixture file records.")
  }
  lock
}

scworkflow_setup_fixtures <- function(
    root = ".",
    cache = Sys.getenv("SCWORKFLOW_FIXTURE_CACHE",
                       unset = tools::R_user_dir("SCWorkflow", "cache")),
    archive = Sys.getenv("SCWORKFLOW_FIXTURE_ARCHIVE", unset = ""),
    download_only = FALSE,
    fixtures = file.path(root, "tests", "testthat", "fixtures")) {
  lock <- scworkflow_fixture_manifest(root)
  if (!is.character(fixtures) || length(fixtures) != 1L || is.na(fixtures) ||
      !nzchar(fixtures)) {
    stop("fixtures must be a single directory path.")
  }
  records <- lock$files
  matches <- function(path, bytes, sha256) {
    file.exists(path) && !dir.exists(path) && file.info(path)$size == bytes &&
      identical(digest::digest(file = path, algo = "sha256"), sha256)
  }
  paths <- file.path(fixtures, records$path)
  present <- file.exists(paths)
  if (!download_only) {
    valid <- vapply(which(present), function(index) {
      matches(paths[[index]], records$bytes[[index]], records$sha256[[index]])
    }, logical(1L))
    if (any(!valid)) {
      stop("Existing fixtures differ from the pinned release: ",
           paste(records$path[which(present)[!valid]], collapse = ", "),
           ". Move those files aside before restoring; no files were overwritten.")
    }
    if (all(present)) {
      return(invisible(fixtures))
    }
  }
  if (!nzchar(archive)) {
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    archive <- file.path(cache, paste0(lock$sha256, "-", lock$asset))
    if (!file.exists(archive)) {
      download_dir <- tempfile("fixture-download-", tmpdir = cache)
      dir.create(download_dir)
      on.exit(unlink(download_dir, recursive = TRUE), add = TRUE)
      piggyback::pb_download(file = lock$asset, dest = download_dir, repo = lock$repo,
                            tag = lock$tag, use_timestamps = FALSE, .token = "")
      download <- file.path(download_dir, lock$asset)
      if (!matches(download, lock$bytes, lock$sha256)) {
        stop("Fixture archive SHA-256 or size mismatch; download rejected.")
      }
      if (!file.rename(download, archive)) {
        stop("Cannot install the verified fixture archive in the cache.")
      }
    }
  }
  if (!matches(archive, lock$bytes, lock$sha256)) {
    stop("Fixture archive SHA-256 or size mismatch: ", archive)
  }
  if (download_only) {
    return(invisible(archive))
  }
  members <- utils::untar(archive, list = TRUE, tar = "internal")
  expected <- paste0("fixtures/", records$path)
  if (!setequal(members, expected) || anyDuplicated(members)) {
    stop("Fixture archive entries do not match the pinned manifest.")
  }
  staging <- tempfile("scworkflow-fixtures-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::untar(archive, exdir = staging, tar = "internal")
  for (index in seq_len(nrow(records))) {
    source <- file.path(staging, expected[[index]])
    if (!matches(source, records$bytes[[index]], records$sha256[[index]])) {
      stop("Extracted fixture checksum mismatch: ", records$path[[index]])
    }
  }
  for (index in which(!present)) {
    dir.create(dirname(paths[[index]]), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(staging, expected[[index]]), paths[[index]],
                   overwrite = FALSE)) {
      stop("Cannot restore fixture: ", paths[[index]])
    }
  }
  message("Restored pinned test fixtures: ", lock$tag)
  invisible(fixtures)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L || (length(args) == 1L && args != "--download-only")) {
    stop("Usage: Rscript tests/setup-fixtures.R [--download-only]")
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  tests <- if (length(script) == 1L) {
    dirname(script)
  } else if (file.exists("fixture-bundle.json")) {
    "."
  } else {
    "tests"
  }
  root <- normalizePath(file.path(tests, ".."), mustWork = TRUE)
  scworkflow_setup_fixtures(root, download_only = "--download-only" %in% args,
                            fixtures = file.path(tests, "testthat", "fixtures"))
}