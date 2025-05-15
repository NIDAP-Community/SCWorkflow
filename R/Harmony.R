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
#' @example Do not run: harmonyBatchCorrect(object = seurat,
#'                                          nvar = 2000,
#'                                          genes.to.add = c("Cd4","Cd8a"),
#'                                          group.by.var = "Mouse_Origin",
#'                                          npc = 20)

#' @return A list: adj.object with harmony-adjusted gene expression (SCT slot) 
#'                 adj.tsne: harmonized tSNE plot

harmonyBatchCorrect <- function(object, 
                                nvar = 2000, 
                                genes.to.add = c(),
                                group.by.var,
                                npc = 20) {
  
  # Error and Warning Messages
  if(is.null(genes.to.add)){
    print("no genes will be added")
  } else if (all(!genes.to.add %in% rownames(object))){
    warning("specified genes were not found and therefore cannot be added")
  }
  
  if (nvar > length(VariableFeatures(object))){
    stop("nvar exceed total number of variable genes in the data")
  }
  
  ## Main Code Block
  # Save original tSNE and UMAP for later comparison
  n <- 60
  qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual.col.pals = qual.col.pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual.col.pals$maxcolors, 
                       rownames(qual.col.pals)))
  
  sdat.tsne.orig <- data.frame(as.vector(object@reductions$tsne@cell.embeddings[,1]),
                               as.vector(object@reductions$tsne@cell.embeddings[,2]),
                               object@meta.data[eval(parse(text = "group.by.var"))])
  names(sdat.tsne.orig) <- c("TSNE1","TSNE2","ident")
  
  sdat.umap.orig <- data.frame(as.vector(object@reductions$umap@cell.embeddings[,1]),
                               as.vector(object@reductions$umap@cell.embeddings[,2]),
                               object@meta.data[eval(parse(text = "group.by.var"))])
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
  seur.SCT <- object@assays$SCT@scale.data
  
  # Add more genes to analyze - must be present in SCT scale.data
  genes.to.add <- genes.to.add[
    genes.to.add %in% rownames(object@assays$SCT@scale.data)]
  
  # Most variable features (mvf) with user added genes
  mvf <- unique(c(VariableFeatures(object)[1:nvar], genes.to.add))
  
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
  object@reductions$pca@cell.embeddings <- ppembed
  object@reductions$pca@feature.loadings <- ppldngs
  object@reductions$pca@stdev <- pppca$d
  
  # By default, Harmony corrects pca embeddings. 
  # Set do_pca to FALSE to use your own pca embeddings. 
  # Stores adjusted embeddings in harmony reduction slot
  
  object <- RunHarmony(object, 
                       group.by.var,
                       do_pca=FALSE,
                       assay.use = "SCT",
                       plot_convergence = FALSE)
  
  object <- RunUMAP(object, reduction = "harmony", dims=1:npc)
  object <- RunTSNE(object, reduction = "harmony", dims=1:npc)
  
  # Plot harmony embeddings annotated by variable to batch correct for
  sdat.tsne <- data.frame(as.vector(object@reductions$tsne@cell.embeddings[,1]),
                          as.vector(object@reductions$tsne@cell.embeddings[,2]),
                          object@meta.data[eval(parse(text = "group.by.var"))])
  names(sdat.tsne) <- c("TSNE1","TSNE2","ident")
  
  sdat.umap <- data.frame(as.vector(object@reductions$umap@cell.embeddings[,1]),
                          as.vector(object@reductions$umap@cell.embeddings[,2]),
                          object@meta.data[eval(parse(text = "group.by.var"))])
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
  harm.embeds <- object@reductions$harmony@cell.embeddings
  harm.lvl.backcalc <- harm.embeds %*% t(ppldngs)
  
  
  # Insert back-calculated data into seurat
  object[["Harmony"]] <- CreateAssayObject(data = t(harm.lvl.backcalc))
  
  return(object)
}
