#' Launch the ModuleScore Shiny App
#'
#' Opens the interactive app to explore ModuleScores and adjust thresholds.
#'
#' @return The result of `shiny::runApp()`
#' @export
launch_module_score_app <- function() {
  appDir <- system.file("shiny/ModuleScoreApp", package = "SCWorkflow")
  if (appDir == "") stop("Shiny app directory not found in SCWorkflow.")
  shiny::runApp(appDir, display.mode = "normal")
}
