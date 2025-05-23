% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Process_Raw_Data.R
\name{processRawData}
\alias{processRawData}
\title{Process Raw Data}
\usage{
processRawData(
  input,
  sample.metadata.table = NULL,
  sample.name.column = NULL,
  organism,
  rename.col = NULL,
  keep = T,
  file.filter.regex = c(),
  split.h5 = F,
  cell.hash = F,
  do.normalize.data = T
)
}
\arguments{
\item{input}{Input can be a vector of scRNA-Seq .h5 files, or a list of
seurat objects for each sample. TCRseq Metadata .csv files can also be
included and will be added to the corrisponding sample seurat object.
Vector of files should include entire file path.}

\item{sample.metadata.table}{A table of sample metadata that you want to
append to the already-existing metadata within the input Seurat Object(s).
(optional)}

\item{sample.name.column}{The column of the input metadata.to.append table
that contains sample names matching the orig.idents in the input object(s).
(optional)}

\item{organism}{Please select species. Choices are Human or Mouse.
(Default: Human).}

\item{rename.col}{Select column name from metadata table that contains new
samples name (optional).}

\item{keep}{If TRUE, keep files when pattern is found in sample name.
If FALSE, remove files when pattern is found in sample name.
The pattern is set in the file.filter.regex parameter (below).}

\item{file.filter.regex}{Pattern or regular expression in sample
name. Use the 'keep' parameter to keep or remove files that contain pattern.
If samples have been renamed set regular expression based on new names}

\item{split.h5}{If TRUE, split H5 into individual files. (Default: FALSE)}

\item{cell.hash}{If TRUE, dataset contains cell hashtags. (Default: FALSE)}

\item{do.normalize.data}{If TRUE counts table will be log2 normalized. If
input contains counts that are already normalzed set to FALSE.
(Default: TRUE)}
}
\value{
Seurat Object and QC plots
}
\description{
Creates Seurat Objects from h5 files for individual or
merged samples. Will log normalize and produce QC figures for
individual samples
}
\details{
This is Step 1 in the basic Single-Cell RNA-seq workflow.
Returns data as a Seurat Object, the basic data structure for
Seurat Single Cell analysis.
}
