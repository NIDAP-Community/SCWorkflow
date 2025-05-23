% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap.R
\name{heatmapSC}
\alias{heatmapSC}
\title{Heatmap of transcript and/or protein expression values in single cells}
\usage{
heatmapSC(
  object,
  sample.names,
  metadata,
  transcripts,
  proteins = NULL,
  heatmap.color = "Bu Yl Rd",
  plot.title = "Heatmap",
  add.gene.or.protein = FALSE,
  protein.annotations = NULL,
  rna.annotations = NULL,
  arrange.by.metadata = TRUE,
  add.row.names = TRUE,
  add.column.names = FALSE,
  row.font = 5,
  col.font = 5,
  legend.font = 5,
  row.height = 15,
  set.seed = 6,
  scale.data = TRUE,
  trim.outliers = TRUE,
  trim.outliers.percentage = 0.01,
  order.heatmap.rows = FALSE,
  row.order = c()
)
}
\arguments{
\item{object}{Seurat-class object}

\item{sample.names}{Sample names}

\item{metadata}{Metadata column to plot}

\item{transcripts}{Transcripts to plot}

\item{proteins}{Proteins to plot (default is NULL)}

\item{heatmap.color}{Color for heatmap. Choices are "Cyan to Mustard",
"Blue to Red", "Red to Vanilla", "Violet to Pink", "Bu Yl Rd",
"Bu Wt Rd" (default is "Bu Yl Rd")}

\item{plot.title}{Title of plot (default is "Heatmap")}

\item{add.gene.or.protein}{Add Gene or protein annotations (default is FALSE)}

\item{protein.annotations}{Protein annotations to add (defulat is NULL)}

\item{rna.annotations}{Gene annotations to add (default is NULL)}

\item{arrange.by.metadata}{Arrange by metadata (default is TRUE)}

\item{add.row.names}{Add row names (default is TRUE)}

\item{add.column.names}{Add column names (default is FALSE)}

\item{row.font}{Font size for rows (default is 5)}

\item{col.font}{Font size for columns (default is 5)}

\item{legend.font}{Font size for legend (default is 5)}

\item{row.height}{Height of row. If NA, adjust to plot size (default is 15)}

\item{set.seed}{Seed for colors (default is 6)}

\item{scale.data}{Perform z-scaling on rows (default is TRUE)}

\item{trim.outliers}{Remove outlier data (default is TRUE)}

\item{trim.outliers.percentage}{Set outlier percentage (default is 0.01)}

\item{order.heatmap.rows}{Order heatmap rows (default is FALSE)}

\item{row.order}{Gene vector to set row order. If NULL, use cluster order
(default is NULL)}
}
\value{
This function returns a heatmap plot and the data underlying the
heatmap.
}
\description{
This method provides a heatmap of single cell data from a Seurat
object given a set of genes and optionally orders by various metadata and/or
gene or protein expression levels. Method is based on ComplexHeatmap::pheatmap
}
