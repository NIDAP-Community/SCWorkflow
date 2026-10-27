arguments <- commandArgs(trailingOnly = TRUE)
tool_sets <- list(
  check = c("rcmdcheck", "testthat"),
  lint = c("goodpractice", "lintr"),
  coverage = c("covr", "testthat")
)

if (length(arguments) != 1L || !arguments[[1]] %in% names(tool_sets)) {
  stop("Usage: Rscript .github/scripts/ci-preflight.R check|lint|coverage", call. = FALSE)
}
mode <- arguments[[1]]
metadata <- read.dcf("DESCRIPTION")
if (!identical(unname(metadata[1, "Package"]), "SCWorkflow")) {
  stop("Run the preflight from the checked-out SCWorkflow package root.", call. = FALSE)
}

dir.create("ci-artifacts", showWarnings = FALSE)
session <- capture.output(sessionInfo())
writeLines(c(
  paste("Source commit:", Sys.getenv("GITHUB_SHA", unset = "not-set")),
  paste("GitHub run:", Sys.getenv("GITHUB_RUN_ID", unset = "local")),
  paste("Validation mode:", mode),
  session
), "ci-artifacts/session-info.txt")
inventory <- as.data.frame(
  utils::installed.packages()[, c("Package", "Version", "LibPath", "Built"), drop = FALSE],
  stringsAsFactors = FALSE
)
utils::write.csv(inventory, "ci-artifacts/installed-packages.csv", row.names = FALSE)
cat(paste(session, collapse = "\n"), "\n")

namespace <- parse("NAMESPACE", keep.source = FALSE)
namespace_packages <- unique(unlist(lapply(namespace, function(entry) {
  if (is.call(entry) && as.character(entry[[1]]) %in% c(
    "import", "importFrom", "importClassesFrom", "importMethodsFrom"
  )) {
    return(as.character(entry[[2]]))
  }
  character()
}), use.names = FALSE))
namespace_dependencies <- data.frame(
  Package = namespace_packages,
  Loadable = vapply(namespace_packages, requireNamespace, logical(1), quietly = TRUE)
)
utils::write.csv(
  namespace_dependencies, "ci-artifacts/namespace-dependencies.csv", row.names = FALSE
)

required <- tool_sets[[mode]]
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
missing_imports <- namespace_dependencies$Package[!namespace_dependencies$Loadable]
issues <- c(
  if (length(missing)) paste("preinstalled tools:", paste(missing, collapse = ", ")),
  if (length(missing_imports)) paste("namespace imports:", paste(missing_imports, collapse = ", "))
)
if (length(issues)) {
  stop(
    "NIDAP CI image cannot validate the checked-out SCWorkflow package; unavailable ",
    paste(issues, collapse = "; "),
    ". Missing namespace imports can prevent package loading before tests run",
    ". No dependency installation was attempted; review the image separately.",
    call. = FALSE
  )
}
cat("NIDAP-image preflight passed for", mode, "\n")
