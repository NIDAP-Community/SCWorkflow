#!/usr/bin/env Rscript --vanilla
# Usage:
#   Rscript .github/install-pak.R [path/to/versions_file]

clargs <- commandArgs(trailingOnly = TRUE)
versions_file <- if (!is.na(clargs[1])) {
  clargs[1]
} else{
  '.github/package-versions.txt'
}
message(cat("versions_file:", versions_file, "\n"))

if (!require('pak', quietly = TRUE)) {
  install.packages('pak', repos = 'https://packagemanager.posit.co/cran/latest')
}


pkgs <- strsplit(readChar(versions_file, file.info(versions_file)$size) |> trimws(),
                 "[[:space:],]+")[[1]]

pak::pkg_install(c(pkgs, "deps::.", "local::."),
                 upgrade = FALSE,
                 dependencies = "all")