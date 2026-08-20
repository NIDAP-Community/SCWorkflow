skip_if_packages_missing <- function(pkgs) {
  for (pkg in pkgs) {
    testthat::skip_if_not_installed(pkg)
  }
}

source(testthat::test_path("..", "..", "R", "Trajectory_CellTypes_v12.R"))

build_test_seurat <- function() {
  set.seed(42)
  counts <- matrix(
    rpois(200, lambda = 5),
    nrow = 20,
    dimnames = list(paste0("g", seq_len(20)), paste0("c", seq_len(10)))
  )
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- suppressWarnings(
    Seurat::FindVariableFeatures(
      seurat_obj,
      selection.method = "vst",
      nfeatures = 10,
      verbose = FALSE
    )
  )
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 5, verbose = FALSE)
  seurat_obj <- Seurat::RunTSNE(
    seurat_obj,
    dims = 1:5,
    perplexity = 2,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    dims = 1:5,
    n.neighbors = 5,
    verbose = FALSE
  )
  clusters <- sample(letters[1:2], ncol(seurat_obj), replace = TRUE)
  names(clusters) <- colnames(seurat_obj)
  Seurat::AddMetaData(
    object = seurat_obj,
    metadata = data.frame(cluster = clusters, row.names = names(clusters))
  )
}

required_pkgs <- c(
  "Seurat", "scater", "scran", "scuttle", "TrajectoryUtils", "TSCAN",
  "ggplot2", "gridExtra", "mclust",
  "SingleCellExperiment", "SummarizedExperiment"
)


 testthat::test_that("Trajectory_CellTypes errors when clustering column missing", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(other = rep("a", ncol(seurat_obj)), row.names = colnames(seurat_obj))
  testthat::expect_error(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = "cluster"
    ),
    "Column 'cluster' not found"
  )
})

 testthat::test_that("Trajectory_CellTypes validates argument types", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(cluster = seurat_obj$cluster, row.names = colnames(seurat_obj))
  testthat::expect_error(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = 123
    ),
    "non-missing character scalar"
  )
})

 testthat::test_that("Trajectory_CellTypes runs on minimal input", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(cluster = seurat_obj$cluster, row.names = colnames(seurat_obj))
  set.seed(7)
  testthat::expect_invisible(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = "cluster",
      Custom_Gene_List = "g1,g2"
    )
  )
})
