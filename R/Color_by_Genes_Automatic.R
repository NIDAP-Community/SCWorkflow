#' @title Color by Gene List
#' @description Returns a panel of reduction plots colored by marker expression
#' @details Takes in a gene table inputted by the user, displays a panel of
#'          tsne, umap, or pca colored by marker expression. The panel will be 
#'          organized in a similar format as the gene table, but with the 
#'          omission of genes not found in the data
#'          
#' @param object Seurat-class object
#' @param samples.subset List of samples to subset the data by
#' @param samples.to.display List of samples to depict on dimension plot, 
#'                           samples not in the list would be colored gray in 
#'                           the background
#' @param marker.table Table of marker genes for each celltype 
#'                    (column names of the table), append "_prot" or "_neg" for 
#'                    proteins or negative markers
#' @param cells.of.interest Celltypes from geneset_dataframe to screen for
#' @param protein.presence Set to TRUE if protein markers are used
#' @param assay Assay to extract gene expression data from (Default: "SCT")
#' @param reduction.type Choose among tsne, umap, and pca (Default: "umap")
#' @param point.transparency Set to lower values for more see through points on 
#'                           dimension plot (Default: 0.5)
#' @param point.shape Change the shape of points for between visualization
#'                    (Default: 16)
#' @param cite.seq Set to TRUE to use CITE-seq embedding for dimension reduction
#'
#' @importFrom Seurat Assays Reductions DimPlot subset
#' @importFrom tidyverse mutate arrange
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggpubr text_grob
#' @importFrom ggplot2 ggplot theme_void theme_bw theme legend.title 
#'  element_blank ggtitle geom_point theme legend.text element_text 
#'  scale_color_gradient xlab ylab aes
#'
#' @export
#' @example Do not run: colorByMarkerTable(object = seurat,
#'                                         samples.subset = c("mouse1","mouse2),
#'                                         samples.to.display = c("mouse1"),
#'                                         marker.table = immuneCellMarkers,
#'                                         cells.of.interest = c("CD4","Treg")
#'                                         )

#' @return arranged grob of dimension reduction plots colored by individual 
#'         marker expression

colorByMarkerTable <- function(object,
                               samples.subset,
                               samples.to.display,
                               marker.table,
                               cells.of.interest,
                               protein.presence = FALSE,
                               assay = "SCT",
                               reduction.type = "umap",
                               point.transparency = 0.5,
                               point.shape = 16,
                               cite.seq = FALSE) {

  # Error Messages
  if (!assay %in% Assays(object)) {
    stop("assay type not found in seurat")
  } else if (!reduction.type %in% Reductions(object)) {
    stop("reduction type not found in seurat")
  }
  
  # Functions
  .plotMarkers <- function(markers) {
    if (is.na(markers) == TRUE) {
      g <- ggplot() + theme_void()
      return(g)
    } else {
      markers.mat = object.sub[[assay]]@scale.data[markers, 
      ]
      markers.quant = quantile(markers.mat[markers.mat > 
                                             1], probs = c(0.1, 0.5, 0.9))
      markers.mat[markers.mat > markers.quant[3]] = markers.quant[3]
      markers.mat[markers.mat < markers.quant[1]] = 0
      if (!(cite.seq)) {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object.sub, reduction = "tsne", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$tSNE_1, 
                               umap2 = p1$data$tSNE_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else if (reduction.type == "umap") {
          p1 <- DimPlot(object.sub, reduction = "umap", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$UMAP_1, 
                               umap2 = p1$data$UMAP_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else {
          p1 <- DimPlot(object.sub, reduction = "pca", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$PC_1, 
                               umap2 = p1$data$PC_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        }
      } else {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object.sub, reduction = "protein_tsne", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_tsne_1, 
                               umap2 = p1$data$protein_tsne_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else if (reduction.type == "umap") {
          p1 <- DimPlot(object.sub, reduction = "protein_umap", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_umap_1, 
                               umap2 = p1$data$protein_umap_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else {
          p1 <- DimPlot(object.sub, reduction = "protein_pca", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_pca_1, 
                               umap2 = p1$data$protein_pca_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        }
      }
      samples.caption <- paste(samples.to.display, sep = "", 
                               collapse = "\n")
      final_caption <- paste("Samples Displayed: ", samples.caption, 
                             sep = "", collapse = "\n")
      clusmat <- dplyr::mutate(clusmat, sample.markers = clusmat$markers * 
                                 grepl(paste(samples.to.display, collapse = "|"), 
                                       clusmat$ident))
      clusmat <- dplyr::arrange(clusmat,sample.markers)
      if (grepl("_neg", markers) == TRUE) {
        clusmat <- dplyr::arrange(clusmat, desc(sample.markers))
        g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                                 group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
          ggtitle(markers) + geom_point(aes(color = sample.markers, 
                                            shape = ident), alpha = point.transparency, 
                                        shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                                                                               panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                               legend.text = element_text(size = rel(0.5))) + 
          scale_color_gradient(limits = c(0, markers.quant[3]), 
                               low = "lightgrey", high = "red") + xlab("umap-1") + 
          ylab("umap-2")
        return(g)
      } else {
        clusmat <- dplyr::arrange(clusmat, sample.markers)
        g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                                 group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
          ggtitle(markers) + geom_point(aes(color = sample.markers, 
                                            shape = ident), alpha = point.transparency, 
                                        shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                                                                               panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                               legend.text = element_text(size = rel(0.5))) + 
          scale_color_gradient(limits = c(0, markers.quant[3]), 
                               low = "lightgrey", high = "red") + xlab("umap-1") + 
          ylab("umap-2")
        return(g)
      }
    }
  }
  if (length(samples.subset) == 0) {
    samples.subset = unique(object@meta.data$sample.name)
  }
  if ("active.ident" %in% slotNames(object)) {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  } else {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  }
  marker.table <- marker.table[cells.of.interest]
  present.marker.ls <- list()
  for (celltype in colnames(marker.table)) {
    print(names(marker.table[celltype]))
    present = lapply(marker.table[[celltype]], function(x) x %in% 
                       rownames(object.sub$SCT@scale.data))
    absent.genes = unlist(marker.table[[celltype]])[present == 
                                                      FALSE]
    present.genes = unlist(marker.table[[celltype]])[present == 
                                                       TRUE]
    print(paste0("Genes not present: ", paste0(absent.genes, 
                                               collapse = ",")))
    print(paste0("Genes present: ", paste0(present.genes, 
                                           collapse = ",")))
    if (length(present.genes) == 0) {
      print(paste0(names(marker.table[celltype]), " genes were not found in object and will not be analyzed"))
    } else {
      present.marker.ls[[celltype]] <- present.genes
    }
  }
  
  padded.ls <- lapply(present.marker.ls, `length<-`, max(lengths(present.marker.ls)))
  markers.from.list <- do.call(cbind, padded.ls)
  markers.present = unlist(markers.from.list)
  
  if (!length(markers.present) > 0) {
    print("No markers found in dataset")
    return(NULL)
  }
  
  gg.storage <- list()
  
  for (cell in colnames(markers.from.list)) {
    title <- cell
    markers.to.analyze <- as.character(markers.from.list[, 
                                                         cell])
    grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))
    gg.storage[[cell]] <- gridExtra::arrangeGrob(grobs = grob, newpage = F, 
                                as.table = F, 
                                top = ggpubr::text_grob(title,size = 15, 
                                                        face = "bold"))
  }
  
  return(gg.storage)
}
