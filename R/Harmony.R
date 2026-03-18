#' @title Harmony Batch Correction from Singular Value Decomposed PCA
#' @description Adjusts cell embeddings and gene expression data to account for 
#'              variations due to user specified variable
#' @details Runs singular value decomposition on pearson residuals 
#'          (SCT scale.data) to obtain PCA embeddings. Performs harmony on 
#'          decomposed embedding and adjusts decomposed gene expression values 
#'          by harmonized embedding. 
#' @param seurat_object Seurat-class object
#' @param nvar Number of variable genes to subset the gene expression data by
#'             (Default: 2000)
#' @param genes.to.add Add genes that might not be found among variably 
#'                     expressed genes
#' @param group.by.var Which variable should be accounted for when running 
#'                     batch correction
#' @param npc Number of principal components to use when running Harmony
#'            (Default: 20)

#' @import Seurat 
#' @import harmony
#' @import gridExtra
#' @import RColorBrewer
#' @import ggplot2
#'   
#' @export
#' @examples
#' \dontrun{
#' harmonyBatchCorrect(
#'   object = seurat,
#'   nvar = 2000,
#'   genes.to.add = c("Cd4", "Cd8a"),
#'   group.by.var = "Mouse_Origin",
#'   npc = 20
#' )
#' }

#' @return A list: adj.object with harmony-adjusted gene expression (SCT slot) 
#'                 adj.tsne: harmonized tSNE plot

harmonyBatchCorrect <- function(object, 
                                nvar = 2000, 
                                genes.to.add = c(),
                                group.by.var,
                                return_lognorm = TRUE,
                                npc = 30) {
  
library(patchwork)  
library(harmony)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

  ensure_single_layer_assay <- function(seurat_obj, assay_name) {
    assay_obj <- seurat_obj[[assay_name]]
    if (inherits(assay_obj, "Assay5")) {
      seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = assay_name)
    }
    seurat_obj
  }

  # Error and Warning Messages
  if (missing(group.by.var) || length(group.by.var) != 1 || !(group.by.var %in% colnames(object@meta.data))) {
    stop("group.by.var must be a single metadata column present in object@meta.data")
  }

  if(is.null(genes.to.add)){
    print("no genes will be added")
  } else if (all(!genes.to.add %in% rownames(object))){
    warning("specified genes were not found and therefore cannot be added")
  }

  assay_use <- if ("SCT" %in% names(object@assays)) "SCT" else DefaultAssay(object)
  object <- ensure_single_layer_assay(object, assay_use)

  var_features <- Seurat::VariableFeatures(object = object, assay = assay_use)
  if (length(var_features) == 0) {
    stop(paste0("No variable features found in assay '", assay_use, "'."))
  }
  
  if (nvar > length(var_features)){
    stop("nvar exceed total number of variable genes in the data")
  }
  
  ## Main Code Block
  # Save original tSNE and UMAP for later comparison
  n <- 60
  qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual.col.pals = qual.col.pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual.col.pals$maxcolors, 
                       rownames(qual.col.pals)))
  
  tsne_embed_orig <- Embeddings(object, reduction = "tsne")
  sdat.tsne.orig <- data.frame(as.vector(tsne_embed_orig[,1]),
                               as.vector(tsne_embed_orig[,2]),
                               object[[group.by.var, drop = TRUE]])
  names(sdat.tsne.orig) <- c("TSNE1","TSNE2","ident")
  
  umap_embed_orig <- Embeddings(object, reduction = "umap")
  sdat.umap.orig <- data.frame(as.vector(umap_embed_orig[,1]),
                               as.vector(umap_embed_orig[,2]),
                               object[[group.by.var, drop = TRUE]])
  names(sdat.umap.orig) <- c("UMAP1","UMAP2","ident")
  
  orig.tsne <- ggplot(sdat.tsne.orig, aes(x=TSNE1, y=TSNE2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=ident),size=1) +
    scale_color_manual(values=cols) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="none",
          panel.background = element_blank(),
          legend.text=element_text(size=rel(0.5))) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    annotate("text", x = Inf, y = -Inf, label = "Original tSNE", hjust = 1.1, vjust = -1, size = 5) 
  
  orig.umap <- ggplot(sdat.umap.orig, aes(x=UMAP1, y=UMAP2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=ident),size=1) +
    scale_color_manual(values=cols) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="none",
          panel.background = element_blank(),
          legend.text=element_text(size=rel(0.5))) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    annotate("text", x = Inf, y = -Inf, label = "Original UMAP", hjust = 1.1, vjust = -1, size = 5)
  
  # Adjusts SCT scale.data based on harmonized embedding
  seur.SCT <- tryCatch(
    Seurat::GetAssayData(object, assay = assay_use, layer = "scale.data"),
    error = function(e) NULL
  )
  if (is.null(seur.SCT) || nrow(seur.SCT) == 0) {
    warning("scale.data layer not available; using data layer for SVD input")
    seur.SCT <- Seurat::GetAssayData(object, assay = assay_use, layer = "data")
  }
  
  # Add more genes to analyze - must be present in SCT scale.data
  genes.to.add <- genes.to.add[genes.to.add %in% rownames(seur.SCT)]
  
  # Most variable features (mvf) with user added genes
  mvf <- unique(c(var_features[1:nvar], genes.to.add))
  
  # Double check that mvf genes are found in scale data genes
  mvf <- mvf[mvf %in% rownames(seur.SCT)]
  
  # Subset SCT scale.data matrix by most variable features
  seur.SCT <- seur.SCT[mvf,]
  seur.SCT <- t(seur.SCT) 
  
  # Singular Value Decomposition (SVD) on scaled counts
  pppca <- svd(seur.SCT) 
  
  # pppca$u: matrix that contains the left singular values
  # pppca$d: vector that contains singular values, sorted in descending order
  ppembed <- pppca$u %*% diag(pppca$d)
  pcnames <- vector(mode = "character")
  for (i in 1:dim(ppembed)[2])pcnames[i] <- paste("PC", i, sep = "_")
  colnames(ppembed) <- pcnames
  rownames(ppembed) <- rownames(seur.SCT)
  
  # pppca$v: matrix that contains the right singular values
  ppldngs <- pppca$v
  colnames(ppldngs) <- pcnames
  rownames(ppldngs) <- mvf 
  
  # Redo pca embeddings based on SVD of gene expression values
  object[["pca"]] <- CreateDimReducObject(
    embeddings = ppembed,
    loadings = ppldngs,
    stdev = pppca$d,
    assay = DefaultAssay(object),
    key = "PC_"
  )

   # Store original log-normalized data and scaling parameters for back-calculation
    if (return_lognorm) {
      library(Matrix)
      # Get log-normalized data for the variable features
      data_layer <- tryCatch(
        Seurat::GetAssayData(object, assay = assay_use, layer = "data"),
        error = function(e) Seurat::GetAssayData(object, assay = assay_use, layer = "counts")
      )
      lognorm_data <- data_layer[mvf, , drop = FALSE]
      
      # Calculate scaling parameters from the original scaled data
      #scale_center <- Matrix::rowMeans(lognorm_data)
      scale_center <- Matrix::rowMeans(as.matrix(lognorm_data))
      scale_scale <- apply(lognorm_data, 1, sd)
      
      # Store these for later reconstruction
      scaling_params <- list(
        center = scale_center,
        scale = scale_scale,
        genes = mvf
      )
    }

  # By default, Harmony corrects pca embeddings. 
  # Set do_pca to FALSE to use your own pca embeddings. 
  # Stores adjusted embeddings in harmony reduction slot
  object <- harmony::RunHarmony(
    object = object,
    group.by.vars = group.by.var,
    reduction.use = "pca",
    dims.use = seq_len(min(npc, ncol(Embeddings(object, reduction = "pca")))),
    plot_convergence = FALSE
  )

  harm_dims <- seq_len(min(npc, ncol(Embeddings(object, reduction = "harmony"))))
  object <- RunUMAP(object, reduction = "harmony", dims = harm_dims)
  object <- RunTSNE(object, reduction = "harmony", dims = harm_dims)
  
  # Plot harmony embeddings annotated by variable to batch correct for
  tsne_embed <- Embeddings(object, reduction = "tsne")
  sdat.tsne <- data.frame(as.vector(tsne_embed[,1]),
                          as.vector(tsne_embed[,2]),
                          object[[group.by.var, drop = TRUE]])
  names(sdat.tsne) <- c("TSNE1","TSNE2","ident")
  
  umap_embed <- Embeddings(object, reduction = "umap")
  sdat.umap <- data.frame(as.vector(umap_embed[,1]),
                          as.vector(umap_embed[,2]),
                          object[[group.by.var, drop = TRUE]])
  names(sdat.umap) <- c("UMAP1","UMAP2","ident")
  
  harm.tsne <- ggplot(sdat.tsne, aes(x=TSNE1, y=TSNE2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=ident),size=1) +
    scale_color_manual(values=cols) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="right",
          panel.background = element_blank(),
          legend.text=element_text(size=rel(1.5))) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    annotate("text", x = Inf, y = -Inf, label = "Harmonized tSNE", hjust = 1.1, vjust = -1, size = 5)
  
  harm.umap <- ggplot(sdat.umap, aes(x=UMAP1, y=UMAP2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=ident),size=1) +
    scale_color_manual(values=cols) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="right",
          panel.background = element_blank(),
          legend.text=element_text(size=rel(1.5))) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    annotate("text", x = Inf, y = -Inf, label = "Harmonized UMAP", hjust = 1.1, vjust = -1, size = 5)
  
  print((orig.tsne + harm.tsne) + plot_layout(ncol = 2))
  print((orig.umap + harm.umap) + plot_layout(ncol = 2))
  
  # Calculate adjusted gene expression from embeddings
  harm.embeds <- Embeddings(object, reduction = "harmony")
  n_common_pcs <- min(ncol(harm.embeds), ncol(ppldngs))
  if (n_common_pcs < 2) {
    stop("Insufficient overlapping PCs between Harmony embeddings and loadings for back-calculation")
  }
  harm.embeds.use <- harm.embeds[, seq_len(n_common_pcs), drop = FALSE]
  ppldngs.use <- ppldngs[, seq_len(n_common_pcs), drop = FALSE]
  harm.lvl.backcalc.scaled <- harm.embeds.use %*% t(ppldngs.use)
  
  # Store batch-corrected scaled data in Harmony assay
  
  if (return_lognorm) {
      # Fast conversion back to log-normalized space
      common_genes <- intersect(colnames(harm.lvl.backcalc.scaled), scaling_params$genes)
      if (length(common_genes) == 0) {
        stop("No overlapping genes available to reconstruct log-normalized values")
      }
      harm_scaled_t <- t(harm.lvl.backcalc.scaled[, common_genes, drop = FALSE])
      harm.lvl.backcalc.lognorm <- sweep(harm_scaled_t, 1, scaling_params$scale[common_genes], `*`)
      harm.lvl.backcalc.lognorm <- sweep(harm.lvl.backcalc.lognorm, 1, scaling_params$center[common_genes], `+`)
      
      print("Batch-corrected data stored in 'Harmony' assay:")
      print("- Log-normalized data: Harmony assay data layer")
      print("- Scaled data: Harmony assay scale.data layer")
    } else {
      harm.lvl.backcalc.lognorm <- t(harm.lvl.backcalc.scaled)
      print("Batch-corrected scaled data was also copied into the Harmony data layer")
    }
  
   # Insert back-calculated data into seurat
   object[["Harmony"]] <- CreateAssayObject(data = harm.lvl.backcalc.lognorm)
   object <- SetAssayData(
     object = object,
     assay = "Harmony",
     layer = "scale.data",
     new.data = t(harm.lvl.backcalc.scaled)
   )

   object <- FindNeighbors(object, reduction = "harmony", dims = seq_len(min(10, length(harm_dims))))
  
  return(
    list("object"=object)
    )
}
