#' @title Annotating cell types using SingleR module
#' @description SingleR is an automatic annotation method for single-cell
#' RNA sequencing (scRNAseq) data (Aran et al. 2019). Given a reference dataset
#' of samples (single-cell or bulk) with known labels, it labels new cells
#' from a test dataset based on similarity to the reference.
#' @details This function is Step 5 of the basic Single-Cell RNA-seq workflow.
#' It is the starting point for downstream visualization, subsetting, and
#' analysis. It takes a combined seurat object as input, such as the one created
#' by the Combined&Renormalized function at the end of the Filter&QC Path
#'
#' @param object Object of class Seurat (your combined Seurat Object
#' after PC reduction has been performed)
#' @param species The species of your samples ("Human" or "Mouse").
#' Default is "Mouse"
#' @param reduction.type  Select the kind of clustering visualization you would
#' like to use to visualize the cell type results ("umap", "tsne", "pca")
#' @param legend.dot.size the size of the colored dots on your chart legend.
#' Default is 2
#' @param do.finetuning Performs the SingleR fine-tuning function.
#' Default is FALSE
#' @param local.celldex Provide a local copy of CellDex library as either an
#' in-memory list object or a path to an .rds file containing references.
#' Default is NULL
#' @param use.clusters Provide cluster identities for each cell.
#' Default is NULL


#'
#' @import Seurat
#' @import cowplot
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom SingleR SingleR
#' @import celldex
#'
#'
#' @export
#'
#' @return a Seurat object with additional metadata



annotateCellTypes <- function(object,
                              species = "Mouse",
                              reduction.type = "umap",
                              legend.dot.size = 2,
                              do.finetuning = FALSE,
                              local.celldex = NULL,
                              use.clusters = NULL) {
  if (is.character(local.celldex)) {
    if (length(local.celldex) != 1 || is.na(local.celldex) || !nzchar(local.celldex)) {
      stop("local.celldex must be NULL, an in-memory CellDex list, or a single valid file path.")
    }
    if (!file.exists(local.celldex)) {
      stop(paste0("local.celldex file not found: ", local.celldex))
    }
    local.celldex <- readRDS(local.celldex)
  }

  .getCounts <- function(so) {
    assay_names <- names(Assays(so))
    if (is.null(assay_names) || length(assay_names) == 0) {
      assay_names <- tryCatch(as.character(Assays(so)), error = function(e) character(0))
    }
    if ("RNA" %in% assay_names) {
      return(GetAssayData(object = so, assay = "RNA", layer = "counts")[, colnames(x = so)])
    }
    return(GetAssayData(object = so, layer = "counts")[, colnames(x = so)])
  }

  .getCelldexRef <- function(local.db, index, name, fetch.fun) {
    if (is.null(local.db)) {
      return(fetch.fun())
    }

    if (is.list(local.db) && !is.null(names(local.db)) && name %in% names(local.db)) {
      return(local.db[[name]])
    }

    if (is.list(local.db) && length(local.db) >= index) {
      return(local.db[[index]])
    }

    stop(paste0("local.celldex does not contain reference '", name, "' (index ", index, ")."))
  }

  ## -------------------------------- ##
  ## Functions                        ##
  ## -------------------------------- ##
  
  
  .annotations <- function(so) {
    so.counts = .getCounts(so)
    if (species == "Human") {

      #HPCA block
      HPCA <- .getCelldexRef(
        local.db = local.celldex,
        index = 1,
        name = "HumanPrimaryCellAtlasData",
        fetch.fun = function() celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
      )
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.main
      )
      if(is.null(use.clusters))
        {
          so[["HPCA_main"]] <- 
            singler$labels[match(rownames(so[[]]), rownames(singler))]
        } else {
          so[["HPCA_main"]] <- 
            singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
        }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["HPCA"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["HPCA"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }
      
      #BP_encode block
      BP <- .getCelldexRef(
        local.db = local.celldex,
        index = 2,
        name = "BlueprintEncodeData",
        fetch.fun = function() celldex::BlueprintEncodeData(ensembl = FALSE)
      )
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.main
      )
      if(is.null(use.clusters))
      {
        so[["BP_encode_main"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["BP_encode_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["BP_encode"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["BP_encode"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }      
    }
    if (species == "Mouse") {
      
      #mouseRNAseq block
      mousernaseq <- .getCelldexRef(
        local.db = local.celldex,
        index = 3,
        name = "MouseRNAseqData",
        fetch.fun = function() celldex::MouseRNAseqData(ensembl = FALSE)
      )
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.main
      )
      if(is.null(use.clusters))
        {
          so[["mouseRNAseq_main"]] <-
            singler$labels[match(rownames(so[[]]), rownames(singler))]
        } else {
          so[["mouseRNAseq_main"]] <- 
            singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
        }   
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.fine
      )
      if(is.null(use.clusters))
        {
          so[["mouseRNAseq"]] <-
            singler$labels[match(rownames(so[[]]), rownames(singler))]
        } else {
          so[["mouseRNAseq"]] <- 
            singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
        }       
      
      
      #ImmGen block
      immgen <- .getCelldexRef(
        local.db = local.celldex,
        index = 4,
        name = "ImmGenData",
        fetch.fun = function() celldex::ImmGenData(ensembl = FALSE)
      )
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.main
      )
      if(is.null(use.clusters))
      {
        so[["immgen_main"]] <- 
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["immgen_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }     
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.fine
      )
      if(is.null(use.clusters))
        {
          so[["immgen"]] <- 
            singler$labels[match(rownames(so[[]]), rownames(singler))]
        } else {
          so[["immgen"]] <- 
            singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
        }     
      
    }
    return(so)
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # Getting Cluster List from MetaData of Seurat Object:
  if(!is.null(use.clusters))
    {
      clusterList <- object@meta.data[[use.clusters]]
    } else {
      clusterList <- NULL
    }
  
  # Running Annotation:
  object <- .annotations(object)
  print("done")
  
  # Assigning Colors (choice of datasets depends on species):
  if (species == "Human") {
    numColors = max(length(unique(object@meta.data$BP_encode_main)), length(unique(object@meta.data$HPCA_main)))
  } else {
    numColors = max(length(unique(object@meta.data$mouseRNAseq_main)), length(unique(object@meta.data$immgen_main)))
  }
  colpaired = colorRampPalette(brewer.pal(12, "Paired"))
  cols = c(
    "#e6194B",
    "#3cb44b",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#fabebe",
    "#469990",
    "#e6beff",
    "#9A6324",
    "#800000",
    "#aaffc3",
    "#808000",
    "#000075",
    colpaired(numColors)
  )
  
  # Creating plots (choice of datasets depends on species):
  if (species == "Human") {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "HPCA_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                     colour = guide_legend(ncol = 4)) + ggtitle("HPCA Main Cell Type Annotations")
    p2 = DimPlot(object, reduction = reduction.type, group.by="BP_encode_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                          colour = guide_legend(ncol = 4)) + ggtitle("BP Encode Main Cell Type Annotations")
  } else {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "immgen_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                       colour = guide_legend(ncol = 4)) + ggtitle("Immgen Main Cell Type Annotations")
    p2 = DimPlot(object, reduction=reduction.type, group.by="mouseRNAseq_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                            colour = guide_legend(ncol = 4)) + ggtitle("Mouse RNAseq Main Cell Type Annotations")
  }
  
  slot(object, "commands") <- list()
  #  cat("\nSingleR Object Checksum:\n")
  #  print(digest::digest(object))
  
  # Returning Seurat Object and 2 plots:
  return(list(
    "object" = object,
    "plots"=list("p1" = p1,
                 "p2" = p2)  
  ))
}
