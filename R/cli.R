# These functions were first inspired by and adapted from renv:
#   https://github.com/rstudio/renv/blob/d0eb86349d35679eb6920ca59072bd7369fe620f/R/cli.R
# and implemented in MOSuite:
#   https://github.com/CCBR/MOSuite/blob/671da78ab50f8cb4b5738bda4f62142e7d9eca36/R/cli.R

#' Execute SCWorkflow from the CLI
#'
#' @export
#' @keywords internal
cli_exec <- function(clargs = commandArgs(trailingOnly = TRUE)) {
  invisible(cli_exec_impl(clargs))
}

cli_exec_impl <- function(clargs) {
  # check for tool called without arguments, or called with '--help'
  usage <-
    length(clargs) == 0 ||
    clargs[1L] %in% c("help", "--help")

  if (usage) {
    return(cli_usage())
  }

  # extract method
  method <- clargs[1L]

  # check request for help on requested method
  help <-
    clargs[2L] %in% c("help", "--help")

  if (help) {
    return(cli_help(method))
  }

  # check for known function in SCWorkflow
  exports <- getNamespaceExports("SCWorkflow")
  if (!method %in% exports) {
    return(cli_unknown(method, exports))
  }

  # begin building call
  # if --json in arguments, call cli_from_json()
  if (any(stringr::str_detect(clargs, "^--json"))) {
    args <- list(call(
      ":::",
      as.symbol("SCWorkflow"),
      as.symbol("cli_from_json")
    ), method = method)
  } else {
    # otherwise call the method directly
    args <- list(call("::", as.symbol("SCWorkflow"), as.symbol(method)))
  }

  for (clarg in clargs[-1L]) {
    # convert '--no-<flag>' into a FALSE parameter
    if (grepl("^--no-", clarg)) {
      key <- substring(clarg, 6L)
      args[[key]] <- FALSE
    }

    # convert '--param=value' flags
    else if (grepl("^--[^=]+=", clarg)) {
      index <- regexpr("=", clarg, fixed = TRUE)
      key <- substring(clarg, 3L, index - 1L)
      val <- substring(clarg, index + 1L)
      args[[key]] <- cli_parse(val)
    }

    # convert '--flag' into a TRUE parameter
    else if (grepl("^--", clarg)) {
      key <- substring(clarg, 3L)
      args[[key]] <- TRUE
    }

    # convert 'param=value' flags
    else if (grepl("=", clarg, fixed = TRUE)) {
      index <- regexpr("=", clarg, fixed = TRUE)
      key <- substring(clarg, 1L, index - 1L)
      val <- substring(clarg, index + 1L)
      args[[key]] <- cli_parse(val)
    }

    # take other parameters as-is
    else {
      args[[length(args) + 1L]] <- cli_parse(clarg)
    }
  }

  # invoke method with parsed arguments
  expr <- as.call(args)
  eval(expr = expr, envir = globalenv())
}

cli_usage <- function(con = stderr()) {
  usage <- "
Usage: SCWorkflow [function] [--json=path/to/args.json]

[function] should be the name of a function exported from SCWorkflow.
[--json] should specify the path to a JSON file with arguments accepted by that function. The equals sign (=) is required to separate --json from the path.

Additionally, the JSON file can contain the following keys:
  - object_input_rds: file path to an existing MultiOmicsDataset object in RDS format. This is required if `method` has `object` as an argument.
  - object_output_rds: file path to write the result to.

Use `scworkflow [function] --help` for more information about the associated function.

Main functions:
  scworkflow tSNE3D
  scworkflow annotateCellTypes
  scworkflow appendMetadataToSeuratObject
  scworkflow colorByGene
  scworkflow colorByMarkerTable
  scworkflow diff_countscombineNormalize
  scworkflow combineNormalize
  scworkflow dotPlotMet
  scworkflow dualLabeling
  scworkflow filterQC
  scworkflow filterSeuratObjectByMetadata
  scworkflow harmonyBatchCorrect
  scworkflow heatmapSC
  scworkflow modScore
  scworkflow nameClusters
  scworkflow plotMetadata
  scworkflow processRawData
  scworkflow reclusterFilteredSeuratObject
  scworkflow reclusterSeuratObject
  scworkflow violinPlot_mod
"
  writeLines(usage, con = con)
}

cli_help <- function(method) {
  print(utils::help(method, package = "SCWorkflow"))
}

cli_unknown <- function(method, exports) {
  # report unknown command
  warning(glue::glue("SCWorkflow: {method} is not a known function."))

  # check for similar commands
  distance <- c(utils::adist(method, exports))
  names(distance) <- exports
  n <- min(distance)
  if (n > 2) {
    return(1L)
  }

  candidates <- names(distance)[distance == n]
  fmt <- "did you mean %s?"
  warning(fmt, paste(shQuote(candidates), collapse = " or "))
  return(1L)
}

cli_parse <- function(text) {
  # handle logical-like values up-front
  if (text %in% c("true", "True", "TRUE")) {
    return(TRUE)
  } else if (text %in% c("false", "False", "FALSE")) {
    return(FALSE)
  }

  # parse the expression
  value <- parse(text = text)[[1L]]
  if (is.language(value))
    text
  else
    value
}

#' Call an SCWorkflow function with arguments specified in a json file
#'
#' @param method function in SCWorkflow to call
#' @param json path to a JSON file containing arguments for the function.
#'  Additionally, the JSON can contain the following keys:
#'    - `object_input_rds` - filepath to an existing MultiOmicsDataset object in RDS format. This is required if the SCWorkflow function contains `object` as an argument.
#'    - `object_output_rds` - filepath to write the result to.
#' @param debug when TRUE, do not call the command, just return the expression.
#'
#' @returns invisible returns the function call
#' @export
#' @keywords internal
#'
cli_from_json <- function(method, json, debug = FALSE) {
  # begin building function call
  fcn_args <- list(call("::", as.symbol("SCWorkflow"), as.symbol(method)))
  # get function arguments from json
  json_args <- jsonlite::read_json(json)

  # if needed, get object from object_input_rds
  accepted_args <- formals(method, envir = getNamespace("SCWorkflow"))
  first_arg <- names(formals(method, envir = getNamespace("SCWorkflow")))[1]
  if (stringr::str_detect(first_arg, glue::glue("^object$"))) {
    assertthat::assert_that(
      "object_input_rds" %in% names(json_args),
      msg = glue::glue(
        "object_input_rds must be included in the JSON because `object` is required for {method}()"
      )
    )
    # most SCWorkflow functions return a list containing an "object" element plus a list of plots or other output.
    # here we extract only the "object" element to pass as the first argument to this function.
    object_list <- readr::read_rds(json_args[["object_input_rds"]])
    if (!("object" %in% names(object_list))) {
      stop('Expected `object_input_rds` to contain an element with the name "object".')
    }
    fcn_args[[first_arg]] <- object_list[["object"]]
  }
  # all other json keys should be arguments for the method
  fcn_args <- c(fcn_args, json_args[!stringr::str_detect(names(json_args), "object_.*_rds")])

  # invoke method with parsed arguments from json
  expr <- as.call(fcn_args)
  if (isFALSE(debug)) {
    result <- eval(expr = expr, envir = globalenv())

    # save result to output_rds
    if ("object_output_rds" %in% names(json_args)) {
      readr::write_rds(result, json_args[["object_output_rds"]])
    }
  }

  invisible(expr)
}
