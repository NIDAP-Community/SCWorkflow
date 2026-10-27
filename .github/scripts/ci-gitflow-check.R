# CI-only adapter for the frozen R runtime and restored CI-tool library.
# Shared Gitflow checks assume Conda and an old renv path, with no setup hook
# for restoring this library, so changing their image alone is insufficient.
# This runner retains push, main, and release validation; deployment stays on
# the existing shared workflow. It does not change SCWorkflow functions.
# This adapter is an integration choice, not a frozen-container requirement;
# upstream Gitflow support for the same setup could replace it.

gitflow_git <- function(arguments) {
  output <- system2("git", shQuote(arguments), stdout = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Gitflow could not determine the changed files.", call. = FALSE)
  }
  output
}

gitflow_base <- function(base_sha, base_ref) {
  if (nzchar(base_sha) && !grepl("^0+$", base_sha)) {
    return(gitflow_git(c("rev-parse", "--verify", paste0(base_sha, "^{commit}")))[[1]])
  }
  if (!nzchar(base_ref)) {
    stop("Gitflow requires a comparison commit or default branch.", call. = FALSE)
  }
  gitflow_git(c("merge-base", "HEAD", paste0("origin/", base_ref)))[[1]]
}

gitflow_scope <- function(event, base_ref) {
  if (event == "push") return("push")
  if (event == "pull_request" && base_ref %in% c("main", "master")) return("main")
  if (event == "pull_request" && startsWith(base_ref, "release")) return("release")
  stop("Unsupported Gitflow validation event or target branch.", call. = FALSE)
}

gitflow_select_files <- function(changed_files, test_files, release = FALSE) {
  excluded <- "(^|/)(fixtures|_snaps)/|(^|/)helper-|(^|/)test-Integration_Test[.]R$"
  test_files <- test_files[!grepl(excluded, test_files, ignore.case = TRUE)]
  source_files <- grep("^R/.*[.][rR]$", changed_files, value = TRUE)
  lint_files <- grep("^(R/|tests/testthat/).*[.][rR]$", changed_files, value = TRUE)
  lint_files <- lint_files[!grepl("(^|/)(fixtures|_snaps)/", lint_files)]
  selected_tests <- intersect(changed_files, test_files)
  for (source_file in source_files) {
    source_name <- tolower(basename(source_file))
    matching <- if (release) {
      grepl(sub("[.][rR]$", "", source_name), tolower(basename(test_files)), fixed = TRUE)
    } else {
      endsWith(tolower(basename(test_files)), source_name)
    }
    selected_tests <- union(selected_tests, test_files[matching])
  }
  list(lint = sort(unique(lint_files)), tests = sort(unique(selected_tests)))
}

gitflow_run_test <- function(path) {
  expression <- paste(
    "pkgload::load_all('.', quiet = TRUE);",
    "testthat::test_file(commandArgs(TRUE)[[1]], reporter = 'summary',",
    "stop_on_failure = TRUE)"
  )
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote(expression), shQuote(path)))
  if (status != 0L) {
    stop("Gitflow unit test failed: ", path, call. = FALSE)
  }
}

run_gitflow_checks <- function() {
  if (as.character(getRversion()) != "4.1.3") {
    stop("Run Gitflow checks only in the pinned R 4.1.3 environment.", call. = FALSE)
  }
  scope <- gitflow_scope(Sys.getenv("GITFLOW_EVENT", unset = "push"),
                         Sys.getenv("GITFLOW_BASE_REF"))
  required <- if (scope == "main") {
    c("devtools", "roxygen2", "rcmdcheck")
  } else {
    c("pkgload", "testthat", if (scope == "push") "lintr")
  }
  for (package in required) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing required Gitflow tool: ", package, call. = FALSE)
    }
  }
  print(.libPaths())
  print(utils::sessionInfo())
  message("Gitflow validation scope: ", scope)
  if (scope == "main") {
    devtools::document()
    devtools::check(error_on = "warning")
    return(invisible(NULL))
  }
  available <- list.files("tests/testthat", pattern = "^test-.*[.][rR]$",
                          recursive = TRUE, full.names = TRUE)
  if (scope == "release") {
    sources <- list.files("R", pattern = "[.][rR]$", full.names = TRUE)
    selected <- gitflow_select_files(sources, available, release = TRUE)
    selected$lint <- character()
  } else {
    base <- gitflow_base(Sys.getenv("GITFLOW_BASE_SHA"), Sys.getenv("GITFLOW_BASE_REF"))
    changed <- gitflow_git(c("-c", "core.quotepath=false", "diff", "--name-only",
                            "--diff-filter=ACMRT", base, "HEAD", "--"))
    selected <- gitflow_select_files(changed, available)
    message("Comparison commit: ", base)
  }
  message("Changed-file lint inputs: ", paste(selected$lint, collapse = ", "))
  for (path in selected$lint) print(lintr::lint(path))
  message("Selected unit tests: ", paste(selected$tests, collapse = ", "))
  if (length(selected$tests) == 0L) {
    message("No changed R sources or unit tests selected a runnable test file.")
    return(invisible(NULL))
  }
  for (path in selected$tests) gitflow_run_test(path)
}

if (sys.nframe() == 0L) run_gitflow_checks()