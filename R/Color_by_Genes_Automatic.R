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
#' @import Seurat
#' @import tidyverse
#' @import gridExtra
#' @import ggpubr
#' @import ggplot2
#'
#' @export
#' @example Do not run: colorByMarkerTable(object = seurat,
  #'                                       samples.subset = c("mouse1","mouse2),
  #'                                       samples.to.display = c("mouse1"),
  #'                                       marker.table = immuneCellMarkers,
  #'                                       cells.of.interest = c("CD4","Treg")
  #'                                       )

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

  .plotMarkers <- function(markers) {
    if (is.na(markers) == TRUE) {
      g <- ggplot() + theme_void()
      return(g)
    } else {
      markers.mat = GetAssayData(object, assay = assay, slot = slot)[markers, ]
      markers.quant = quantile(markers.mat[markers.mat > 
                                             1], probs = c(0.1, 0.5, 0.9))
      markers.mat[markers.mat > markers.quant[3]] = markers.quant[3]
      markers.mat[markers.mat < markers.quant[1]] = 0
      if (!(cite.seq)) {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object, reduction = "tsne", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$tSNE_1, 
                               umap2 = p1$data$tSNE_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else if (reduction.type == "umap") {
          p1 <- DimPlot(object, reduction = "umap", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$UMAP_1, 
                               umap2 = p1$data$UMAP_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else {
          p1 <- DimPlot(object, reduction = "pca", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$PC_1, 
                               umap2 = p1$data$PC_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        }
      } else {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object, reduction = "protein_tsne", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_tsne_1, 
                               umap2 = p1$data$protein_tsne_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else if (reduction.type == "umap") {
          p1 <- DimPlot(object, reduction = "protein_umap", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_umap_1, 
                               umap2 = p1$data$protein_umap_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        } else {
          p1 <- DimPlot(object, reduction = "protein_pca", 
                        group.by = "ident")
          clusmat = data.frame(umap1 = p1$data$protein_pca_1, 
                               umap2 = p1$data$protein_pca_2, markers = markers.mat, 
                               ident = as.factor(p1$data$ident))
        }
      }
      
      g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                               group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
        ggtitle(markers) + geom_point(aes(color = markers, 
                                          shape = ident), alpha = point.transparency, 
                                      shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                                                                             panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                             legend.text = element_text(size = rel(0.5)), plot.background=element_rect(fill="white", colour="white"), rect = element_rect(fill = 'white')) +
        scale_color_gradient(limits = c(0, markers.quant[3]), 
                             low = "lightgrey", high = "red") + xlab("umap-1") + 
        ylab("umap-2")
      return(g)
      
    }
  }
  
  marker.table <- marker.table[cells.of.interest]
  present.marker.ls <- list()
  for (celltype in colnames(marker.table)) {
    print(names(marker.table[celltype]))
    present = lapply(marker.table[[celltype]], function(x) x %in% 
                       rownames(GetAssayData(object, assay = assay, slot = slot)))
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
  
  # Plot manual genes if not empty
  if(!is.null(manual.genes)){
    
    # Str-spit and use only present genes
    manual.genes.processed <- str_split(manual.genes, pattern = "[^a-zA-Z0-9]+")[[1]]
    manual.genes.processed <- manual.genes.processed[manual.genes.processed %in% rownames(GetAssayData(object, assay = assay, slot = slot))]
    
    tryCatch({
      manual.grob <- lapply(manual.genes.processed, function(x) .plotMarkers(x))
      manual.arranged <- gridExtra::arrangeGrob(grobs = manual.grob, newpage = F, as.table = F, top = ggpubr::text_grob("Manual Genes", size = 15, face = "bold"))
      plot(manual.arranged)
      grid.newpage()}, error = function(e) {
        cat(e$message, "\n", "Possible Reason: No manual genes were found in expression matrix")
      })
  }
  
  # Store images for making consolidated figure
  cons.gg.storage <- list()
  
  # Plot arranged figures with padding
  for (cell in colnames(markers.from.list)) {
    title <- cell
    markers.to.analyze <- as.character(markers.from.list[, 
                                                         cell])
    grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))
    arranged <- gridExtra::arrangeGrob(grobs = grob, newpage = F, as.table = F, top = ggpubr::text_grob(title, size = 15, face = "bold"))
    
    cons.gg.storage[[cell]] <- arranged
    
    # Manually specify background 
    background <- rectGrob(gp = gpar(fill = "white", col = NA)) 
    # Combine the background with the arranged plots
    combined <- grobTree(background, arranged)
    # Draw the combined grob with background
    grid.draw(combined)
    grid.newpage()
  }
  
  if(consolidate.marker.fig == TRUE){
    cons.fig <- do.call(arrangeGrob, c(cons.gg.storage, ncol = ncol(markers.from.list)))
    cons.fig.bkgrd <- grobTree(background, cons.fig)
    grid.draw(cons.fig.bkgrd)
  }
  
  return(NULL)
}
