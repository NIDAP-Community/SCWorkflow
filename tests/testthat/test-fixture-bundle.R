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

make_fixture_package <- function() {
  bundle <- make_fixture_bundle()
  dir.create(file.path(bundle$root, "R"))
  dir.create(file.path(bundle$root, "tests", "testthat"))
  writeLines(c("Package: scworkflowfixtureprobe", "Version: 0.0.0.9000",
               "Title: Fixture Bootstrap Probe",
               "Description: Temporary package for fixture bootstrap runner tests.",
               "Author: Fixture Probe", "Maintainer: Fixture Probe <fixture@example.org>",
               "License: GPL-3", "Imports: digest, jsonlite", "Suggests: testthat"),
             file.path(bundle$root, "DESCRIPTION"))
  writeLines("export(fixture_value)", file.path(bundle$root, "NAMESPACE"))
  writeLines(c("fixture_value <- function() {", "  'small fixture'", "}"),
             file.path(bundle$root, "R", "fixture.R"))
  stopifnot(file.copy(testthat::test_path("..", "setup-fixtures.R"),
                     file.path(bundle$root, "tests", "setup-fixtures.R")))
  stopifnot(file.copy(testthat::test_path("helper-00-fixtures.R"),
                     file.path(bundle$root, "tests", "testthat", "helper-00-fixtures.R")))
  writeLines(c("library(testthat)", "library(scworkflowfixtureprobe)",
               "test_check('scworkflowfixtureprobe')"),
             file.path(bundle$root, "tests", "testthat.R"))
  writeLines(c("test_that('runner reads restored fixture', {",
               "  expect_identical(readLines(test_path('fixtures/cohort', 'input.txt')),",
               "                   fixture_value())", "})"),
             file.path(bundle$root, "tests", "testthat", "test-active-fixtures.R"))
  bundle
}

run_fixture_package <- function(bundle, expression) {
  log <- file.path(bundle$root, "runner.log")
  previous <- setwd(bundle$root)
  on.exit(setwd(previous))
  status <- system2(file.path(R.home("bin"), "Rscript"),
                     c("--vanilla", "-e", shQuote(expression)), stdout = log, stderr = log,
                     env = c(paste0("SCWORKFLOW_FIXTURE_ARCHIVE=", shQuote(bundle$archive)),
                             "TESTTHAT=false"))
  list(status = status, output = paste(readLines(log, warn = FALSE), collapse = "\n"))
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

test_that("fixture bootstrap restores only the active relocated fixture directory", {
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  fixtures <- file.path(bundle$root, "SCWorkflow-tests", "testthat", "fixtures")
  expect_identical(scworkflow_setup_fixtures(bundle$root, archive = bundle$archive,
                                            fixtures = fixtures), fixtures)
  expect_identical(readLines(file.path(fixtures, "cohort", "input.txt")),
                   "small fixture")
  expect_false(dir.exists(file.path(bundle$root, "tests", "testthat", "fixtures")))
  unlink(bundle$archive)
  expect_identical(scworkflow_setup_fixtures(bundle$root, fixtures = fixtures), fixtures)
})

test_that("fixture helpers load before the Gitflow testthat context exists", {
  bundle <- make_fixture_package()
  on.exit(unlink(bundle$root, recursive = TRUE))
  result <- run_fixture_package(bundle, paste(
    "pkgload::load_all('.', quiet = TRUE);",
    "testthat::test_file('tests/testthat/test-active-fixtures.R',",
    "reporter = 'summary', stop_on_failure = TRUE)"
  ))
  expect_equal(result$status, 0L, info = result$output)
  expect_identical(readLines(file.path(bundle$root, "tests", "testthat", "fixtures",
                                       "cohort", "input.txt")), "small fixture")
})

test_that("installed tests run with fixtures in the active copied tree", {
  bundle <- make_fixture_package()
  on.exit(unlink(bundle$root, recursive = TRUE))
  library <- file.path(bundle$root, "library")
  result <- run_fixture_package(bundle, paste0(
    "library <- ", encodeString(library, quote = '"'), ";",
    "dir.create(library);",
    "utils::install.packages('.', lib = library, repos = NULL, type = 'source',",
    "INSTALL_opts = '--install-tests');",
    ".libPaths(c(library, .libPaths()));",
    "Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep));",
    "status <- tools::testInstalledPackage('scworkflowfixtureprobe',",
    "lib.loc = library, outDir = file.path(library, 'scworkflowfixtureprobe'),",
    "types = 'tests');",
    "stopifnot(status == 0L)"
  ))
  expect_equal(result$status, 0L, info = result$output)
  package <- file.path(library, "scworkflowfixtureprobe")
  active <- file.path(package, "scworkflowfixtureprobe-tests", "testthat", "fixtures")
  expect_identical(readLines(file.path(active, "cohort", "input.txt")), "small fixture")
  expect_false(dir.exists(file.path(package, "tests", "testthat", "fixtures")))
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

test_that("fixture bootstrap validates directory argument types", {
  expect_error(scworkflow_setup_fixtures(42), "single directory path")
  expect_error(scworkflow_setup_fixtures(c("one", "two")), "single directory path")
  bundle <- make_fixture_bundle()
  on.exit(unlink(bundle$root, recursive = TRUE))
  expect_error(scworkflow_setup_fixtures(bundle$root, fixtures = 42), "single directory path")
  expect_error(scworkflow_setup_fixtures(bundle$root, fixtures = character()),
               "single directory path")
  expect_error(scworkflow_setup_fixtures(bundle$root, fixtures = ""), "single directory path")
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