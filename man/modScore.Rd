% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModuleScore.R
\name{modScore}
\alias{modScore}
\title{Compute ModScore}
\usage{
modScore(
  object,
  marker.table,
  ms.threshold,
  general.class,
  lvl.vec = c(),
  reduction = "tsne",
  nbins = 10,
  gradient.ft.size = 6,
  violin.ft.size = 6,
  step.size = 0.1
)
}
\arguments{
\item{object}{Seurat-class object}

\item{marker.table}{Table of marker genes for each celltype
(column names of the table), append "_prot" or
"_neg" for proteins or negative markers}

\item{general.class}{Base population of cells to classify}

\item{reduction}{Choose among tsne, umap, and pca (Default: tsne)}

\item{nbins}{Number of bins for storing control features and analyzing
average expression (Default: 10)}

\item{gradient.ft.size}{Set size of axis labels on gradient
density plot of ModuleScore distribution
(Default: 6)}

\item{violin.ft.size}{Set size of axis labels on violin plot of
ModuleScore distribution (Default: 6)}

\item{step.size}{Set step size of distribution plots (Default: 0.1)}

\item{cite.seq}{Set to TRUE if there are CITE-seq markers in
marker.table (Default: FALSE)}

\item{celltypes}{Vector of celltypes from marker.table to
screen for}

\item{threshold}{Specify bimodal thresholds for cell classification,
should be of the same length as celltypes
vector}

\item{multi.lvl}{Toggle to TRUE if there are subpopulations of cells
you want to screen for (Default: FALSE)}

\item{lvl.df}{Table of subpopulation levels and parent-child
information (e.g. Tcells-CD4, Tcells-CD8)}
}
\value{
List containing annotated dimension plot with ModuleScore
distribution of cell marker gene, Seurat Object with cell
classification metadata
}
\description{
Returns Seurat-class object with metadata containing
ModuleScores and Likely_CellType calls
}
\details{
Analyzed features are binned based on averaged expression;
control features are randomly selected from each bin.
}
