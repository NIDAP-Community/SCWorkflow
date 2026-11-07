source(testthat::test_path("..", "setup-fixtures.R"), local = TRUE)

make_fixture_bundle <- function() {
  root <- tempfile("fixture-bootstrap-test-")
  dir.create(file.path(root, "tests"), recursive = TRUE)
  payload <- file.path(root, "payload")
  dir.create(file.path(payload, "fixtures", "cohort"), recursive = TRUE)
  fixture <- file.path(payload, "fixtures", "cohort", "input.txt")
  writeLines("small fixture", fixture)
  archive <- file.path(root, "bundle.tar.gz")
  previous <- setwd(payload)
  tryCatch(utils::tar(archive, "fixtures/cohort/input.txt", compression = "gzip",
                     tar = "internal"), finally = setwd(previous))
  lock <- list(schema = 1L, repo = "NIDAP-Community/SCWorkflow", tag = "test-v1",
               asset = "bundle.tar.gz", bytes = unname(file.info(archive)$size),
               sha256 = digest::digest(file = archive, algo = "sha256"),
               files = list(list(path = "cohort/input.txt",
                                 bytes = unname(file.info(fixture)$size),
                                 sha256 = digest::digest(file = fixture,
                                                         algo = "sha256"))))
  manifest <- file.path(root, "tests", "fixture-bundle.json")
  jsonlite::write_json(lock, manifest, auto_unbox = TRUE)
  list(root = root, archive = archive, manifest = manifest, lock = lock)
}

test_that("fixture bootstrap restores exact files and reuses them offline", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  restored <- scworkflow_setup_fixtures(bundle$root, archive = bundle$archive)
  expect_identical(readLines(file.path(restored, "cohort", "input.txt")),
                   "small fixture")
  unlink(bundle$archive)
  expect_identical(scworkflow_setup_fixtures(bundle$root), restored)
})

test_that("fixture bootstrap resolves R CMD BATCH paths without a file argument", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  script <- file.path(bundle$root, "tests", "setup-fixtures.R")
  stopifnot(file.copy(testthat::test_path("..", "setup-fixtures.R"), script))
  fixtures <- file.path(bundle$root, "tests", "testthat", "fixtures")
  for (directory in c(bundle$root, file.path(bundle$root, "tests"))) {
    unlink(fixtures, recursive = TRUE)
    log <- file.path(bundle$root, paste0("batch-", basename(directory), ".Rout"))
    previous <- setwd(directory)
    status <- tryCatch(
      system2(file.path(R.home("bin"), "R"),
              c("CMD", "BATCH", "--vanilla", shQuote(script), shQuote(log)),
              env = paste0("SCWORKFLOW_FIXTURE_ARCHIVE=", shQuote(bundle$archive))),
      finally = setwd(previous)
    )
    expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
    expect_identical(readLines(file.path(fixtures, "cohort", "input.txt")),
                     "small fixture")
  }
})

test_that("fixture bootstrap retains Rscript paths and download-only behavior", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  script <- file.path(bundle$root, "tests", "setup-fixtures.R")
  stopifnot(file.copy(testthat::test_path("..", "setup-fixtures.R"), script))
  fixtures <- file.path(bundle$root, "tests", "testthat", "fixtures")
  log <- file.path(bundle$root, "rscript.log")
  previous <- setwd(file.path(bundle$root, "payload"))
  status <- tryCatch(
    system2(file.path(R.home("bin"), "Rscript"),
            c(shQuote(script), "--download-only"), stdout = log, stderr = log,
            env = paste0("SCWORKFLOW_FIXTURE_ARCHIVE=", shQuote(bundle$archive))),
    finally = setwd(previous)
  )
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  expect_false(dir.exists(fixtures))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                     stdout = log, stderr = log,
                     env = paste0("SCWORKFLOW_FIXTURE_ARCHIVE=", shQuote(bundle$archive)))
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  expect_identical(readLines(file.path(fixtures, "cohort", "input.txt")),
                   "small fixture")
})

test_that("fixture bootstrap rejects missing manifest fields", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  bundle$lock$sha256 <- NULL
  jsonlite::write_json(bundle$lock, bundle$manifest, auto_unbox = TRUE)
  expect_error(scworkflow_setup_fixtures(bundle$root), "missing required fields")
})

test_that("fixture bootstrap validates root argument types", {
  expect_error(scworkflow_setup_fixtures(42), "single directory path")
  expect_error(scworkflow_setup_fixtures(c("one", "two")), "single directory path")
})

test_that("fixture archive can be verified before an offline test run", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  expect_identical(scworkflow_setup_fixtures(bundle$root, archive = bundle$archive,
                                            download_only = TRUE), bundle$archive)
  expect_false(dir.exists(file.path(bundle$root, "tests", "testthat", "fixtures")))
})

test_that("fixture bootstrap verifies extracted files before installing", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  bundle$lock$files[[1L]]$sha256 <- strrep("0", 64L)
  jsonlite::write_json(bundle$lock, bundle$manifest, auto_unbox = TRUE)
  expect_error(scworkflow_setup_fixtures(bundle$root, archive = bundle$archive),
               "Extracted fixture checksum mismatch")
  expect_false(dir.exists(file.path(bundle$root, "tests", "testthat", "fixtures")))
})

test_that("fixture bootstrap rejects bad archives without installing files", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  writeLines("not the pinned archive", bundle$archive)
  expect_error(scworkflow_setup_fixtures(bundle$root, archive = bundle$archive),
               "SHA-256 or size mismatch")
  expect_false(dir.exists(file.path(bundle$root, "tests", "testthat", "fixtures")))
})

test_that("fixture bootstrap rejects unsafe paths and altered existing files", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  restored <- scworkflow_setup_fixtures(bundle$root, archive = bundle$archive)
  writeLines("local change", file.path(restored, "cohort", "input.txt"))
  expect_error(scworkflow_setup_fixtures(bundle$root), "no files were overwritten")
  expect_identical(readLines(file.path(restored, "cohort", "input.txt")), "local change")
  bundle$lock$files[[1L]]$path <- "../escape.txt"
  jsonlite::write_json(bundle$lock, bundle$manifest, auto_unbox = TRUE)
  expect_error(scworkflow_setup_fixtures(bundle$root), "Invalid fixture file records")
})