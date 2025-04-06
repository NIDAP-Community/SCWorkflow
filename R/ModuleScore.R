#' @title Compute ModScore
#' @description Returns Seurat-class object with metadata containing 
#'              ModuleScores and Likely_CellType calls
#' @details Analyzed features are binned based on averaged expression; 
#'          control features are randomly selected from each bin.
#' 
#' @param object Seurat-class object
#' @param marker.table Table of marker genes for each celltype 
#'                          (column names of the table), append "_prot" or 
#'                          "_neg" for proteins or negative markers
#' @param cite.seq Set to TRUE if there are CITE-seq markers in 
#'                          marker.table (Default: FALSE)
#' @param celltypes Vector of celltypes from marker.table to 
#'                             screen for
#' @param threshold Specify bimodal thresholds for cell classification, 
#'                         should be of the same length as celltypes 
#'                         vector
#' @param general.class Base population of cells to classify
#' @param multi.lvl Toggle to TRUE if there are subpopulations of cells 
#'                          you want to screen for (Default: FALSE)
#' @param lvl.df Table of subpopulation levels and parent-child 
#'                         information (e.g. Tcells-CD4, Tcells-CD8)
#' @param reduction Choose among tsne, umap, and pca (Default: tsne)
#' @param nbins Number of bins for storing control features and analyzing 
#'              average expression (Default: 10)
#' @param gradient.ft.size Set size of axis labels on gradient 
#'                                   density plot of ModuleScore distribution
#'                                   (Default: 6)
#' @param violin.ft.size Set size of axis labels on violin plot of 
#'                             ModuleScore distribution (Default: 6)
#' @param step.size Set step size of distribution plots (Default: 0.1)
#' 
#' @import Seurat
#' @import tidyverse
#' @import gridExtra
#' @import quantmod
#' @import grid
#' @import data.table
#' @import utils
#' @importFrom dplyr select
#'   
#' @export
#' @example Do not run: moduleScore(object = seurat,
#'                                  marker.table = immuneCellMarkers,
#'                                  celltypes = c("CD4_T","Treg",Monocytes"),
#'                                  threshold = c(0.1,0.4, 0.3),
#'                                  multi.lvl = FALSE
#'                                  )
#'                                  
#' @example Do not run: moduleScore(object = seurat,
#'                                  marker.table = immuneCellMarkers,
#'                                  celltypes = c("CD4_T","Treg",Monocytes"),
#'                                  threshold = c(0.1,0.4, 0.3),
#'                                  general.class = c("CD_T","Monocytes"),
#'                                  multi.lvl = TRUE,
#'                                  lvl.df = parentChildTable
#'                                  )

#' @return List containing annotated dimension plot with ModuleScore 
#'         distribution of cell marker gene, Seurat Object with cell 
#'         classification metadata

modScore <- function(object, marker.table, ms.threshold, 
                     general.class, lvl.vec = c(), reduction = "tsne", 
                     nbins = 10, gradient.ft.size = 6, 
                     violin.ft.size = 6, step.size = 0.1) 
{
  library(Seurat)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  
  # Function for separating and calling cells by bimodal thresholds
  .modScoreCall <- function(ms.meta, numeric_threshold, reject) {
    thres.ls <- list()
    for (i in 1:ncol(ms.meta)) {
      thres.ls[[i]] <- rep(numeric_threshold[i], nrow(ms.meta))
    }
    thres.df <- data.frame(matrix(unlist(thres.ls), nrow = nrow(ms.meta)))
    thres.filter <- ms.meta > thres.df
    ms.meta.filt <- ms.meta * thres.filter
    max.col.vec <- max.col(ms.meta.filt)
    zero.filt <- as.integer(!apply(ms.meta.filt, 1, function(find_zero_rows) all(find_zero_rows == 0)))
    final.filt <- (max.col.vec * zero.filt) + 1
    append.name <- c(reject, names(ms.meta))
    dupl.data <- ms.meta
    dupl.data[, "MS_Celltype"] <- append.name[final.filt]
    return(dupl.data)
  }
  
  # Upstream processing
  # String split celltype_thresholds - numeric portion
  numeric_threshold <- sapply(stringr::str_split(ms.threshold, " "), function(x) as.numeric(x[2]))
  if (!"Barcode" %in% colnames(object@meta.data)) {
    object@meta.data$Barcode <- rownames(object@meta.data)
  }
  colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", 
                                     colnames(object@meta.data))
  
  # Marker table processing
  marker.tab <- unlist(marker.table)
  celltypes <- sapply(str_split(ms.threshold, " "), function(x) as.character(x[1]))
  marker = select(marker.table, celltypes)
  marker.list = as.list(marker)
  if (sum(unlist(marker.list) %in% rownames(object@assays$SCT@data)) == 
      0) {
    stop("No genes from list was found in data")
  }
  if (length(numeric_threshold) != length(celltypes)) {
    if (sum(numeric_threshold) == 0) {
      numeric_threshold <- rep(0, length(celltypes))
      print("Module Score threshold set to zero - outputing preliminary data")
    } else {
      stop("Threshold length does not match # celltypes to analyze")
    }
  }
  
  # For each celltype, print out present / nonpresent genes, calculate MS and generate plots
  names(numeric_threshold) <- celltypes
  figures <- list()
  exclude_cells <- c()
  h = 0
  j = 1
  for (h in seq_along(marker.list)) {
    print(names(marker.list[h]))
    present = lapply(marker.list[[h]], function(x) x %in% 
                       rownames(object@assays$SCT@data))
    absentgenes = unlist(marker.list[[h]])[present == FALSE]
    absentgenes = absentgenes[is.na(absentgenes) == F]
    presentgenes = unlist(marker.list[[h]])[present == TRUE]
    presentgenes = presentgenes[is.na(presentgenes) == F]
    print(paste0("Genes not present: ", paste0(absentgenes, 
                                               collapse = ",")))
    print(paste0("Genes present: ", paste0(presentgenes, 
                                           collapse = ",")))
    if (length(presentgenes) == 0) {
      print(paste0(names(marker.list[h]), " genes were not found in object and will not be analyzed"))
      exclude_cells[j] <- h
      j = j + 1
    }
  }
  # End of check present / absent genes
  
  if (length(exclude_cells) > 0) {
    marker.list <- marker.list[-exclude_cells]
  } else {
    marker.list <- marker.list
  }
  
  # clean up list, remove NAs for faster run
  marker.list <- lapply(marker.list, na.omit)
  
  # Calculate MS, make density plots
  for (celltype_name in names(marker.list)) {
    object = AddModuleScore(object, marker.list[celltype_name], 
                            name = celltype_name, nbin = nbins, assay = "SCT")
    m = paste0(celltype_name, "1")
    object@meta.data[[m]] <- scales::rescale(object@meta.data[[m]], 
                                             to = c(0, 1))
    
    # Do plots for just general (parent) celltypes
    if (celltype_name %in% general.class){
      clusid = object@meta.data[[m]]
      d <- density(clusid)
      
      # Make dimension reduction tables and plots
      if (reduction == "tsne") {
        p1 <- DimPlot(object, reduction = "tsne", group.by = "ident")
      } else if (reduction == "umap") {
        p1 <- DimPlot(object, reduction = "umap", group.by = "ident")
      } else {
        p1 <- DimPlot(object, reduction = "pca", group.by = "ident")
      }
      
      if (reduction == "tsne") {
        clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$tSNE_1, umap2 = p1$data$tSNE_2, clusid = as.numeric(object@meta.data[[m]]))
      } else if (reduction == "umap") {
        clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$UMAP_1, umap2 = p1$data$UMAP_2, clusid = as.numeric(object@meta.data[[m]]))
      } else {
        clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$PC_1, 
                             umap2 = p1$data$PC_2, clusid = as.numeric(object@meta.data[[m]]))}
      
      clusmat <- mutate(clusmat, sample_clusid = clusmat$clusid)
      umap.pos <- clusmat %>% group_by(clusid) %>% dplyr::summarise(umap1.mean = mean(umap1), umap2.mean = mean(umap2))
      title = as.character(m)
      clusmat <- clusmat %>% dplyr::arrange(clusid)
      clusid.df <- data.frame(id = object@meta.data$orig.ident, 
                              ModuleScore = object@meta.data[[m]])
      
      g <- ggplot(clusmat, aes(x = umap1, y = umap2)) + theme_bw() + 
        theme(legend.title = element_blank()) + geom_point(aes(colour = sample_clusid), alpha = 0.5, shape = 20, size = 1) + scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(0, 
                                                                                                                                                                                                                        numeric_threshold[celltype_name]/2, numeric_threshold[celltype_name], (numeric_threshold[celltype_name] + 1)/2, 
                                                                                                                                                                                                                        1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                 panel.grid.minor = element_blank(), panel.background = element_blank()) + xlab("tsne-1") + ylab("tsne-2")
      
      g1 <- RidgePlot(object, features = m, group.by = "orig.ident") + 
        theme(legend.position = "none", title = element_blank(), 
              axis.text.x = element_text(size = gradient.ft.size)) + 
        geom_vline(xintercept = numeric_threshold[celltype_name], linetype = "dashed", 
                   color = "red3") + scale_x_continuous(breaks = seq(0, 
                                                                     1, step.size))
      
      g2 <- ggplot(clusid.df, aes(x = id, y = ModuleScore)) + 
        geom_violin(aes(fill = id)) + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = rel(0.8)), legend.position = "top", axis.text.y = element_text(size = violin.ft.size)) + 
        guides(colour = guide_legend(override.aes = list(size = 5, 
                                                         alpha = 1))) + geom_hline(yintercept = numeric_threshold[celltype_name], 
                                                                                   linetype = "dashed", color = "red3") + scale_y_continuous(breaks = seq(0, 1, step.size))
      
      g3 <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + 
        xlab("ModuleScore") + ylab("Density") + geom_line() + 
        geom_segment(aes(xend = d$x, yend = 0, colour = x)) + 
        scale_y_log10() + scale_color_gradientn(colours = c("blue4", 
                                                            "lightgrey", "red"), values = scales::rescale(c(0, 
                                                                                                            numeric_threshold[celltype_name]/2, numeric_threshold[celltype_name], (numeric_threshold[celltype_name] + 1)/2, 
                                                                                                            1), limits = c(0, 1))) + geom_vline(xintercept = numeric_threshold[celltype_name], 
                                                                                                                                                linetype = "dashed", color = "red3") + geom_vline(xintercept = numeric_threshold[celltype_name], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0, 1, step.size)) + theme(legend.title = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                  axis.text.x = element_text(size = 6))
      
      figures[[celltype_name]] = arrangeGrob(g, g1, g2, g3, ncol = 2, top = textGrob(paste0(celltype_name," (General Class)"), gp = gpar(fontsize = 14, fontface = "bold")))
    }
  }
  
  # Rename MS columns - get rid of "1" at the end
  colnames(object@meta.data)[colnames(object@meta.data) %in% 
                               paste0(names(marker.list), 1)] <- names(marker.list)
  
  # First annotate general level of celltypes
  general.class <- general.class[general.class %in% colnames(object@meta.data)]
  trunc.meta.gen <- object@meta.data[general.class]
  gen.thrs.vec <- numeric_threshold[general.class]
  call.res <- .modScoreCall(trunc.meta.gen, gen.thrs.vec, reject = "unknown")
  call.res$Barcode <- rownames(call.res)
  
  ### Hierarchical Annotation ###
  # Specific celltypes
  # convert lvl_vec to lvl.df
  if(length(lvl.vec) > 0){
    cell_pop <- str_split(lvl.vec, '-')
    
    # cell_pop: a list of str_split elements from lvl.vec
    entry_list = vector("list", length = length(cell_pop))
    
    for(i in 1:length(cell_pop)){
      # pair cells together
      cnt = 2
      vec = c()
      
      # as long as there are non NAs left, keep updating vec
      while(!is.na(cell_pop[[i]][cnt])){
        
        entry = paste0(cell_pop[[i]][cnt-1],'-',cell_pop[[i]][cnt])
        vec = c(vec, entry)
        cnt = cnt + 1
        
      }
      
      entry_list[[i]] = vec
    }
    
    ### Make levels data.frame ###
    # Pad the shorter entries with NA
    max_len <- max(sapply(entry_list, length))
    padded_list <- lapply(entry_list, function(x) {
      length(x) <- max_len  # implicitly pads with NAs
      return(x)
    })
    
    lvl.df <- as.data.frame(do.call(rbind, padded_list), stringsAsFactors = FALSE)
    
    # clean up lvl.df by removing duplicated entries across columns
    lvl.df <- as.data.frame(lapply(lvl.df, function(col) {
      col[duplicated(col)] <- NA
      col
    }))
    
    for (k in 1:ncol(lvl.df)) {
      
      sub.class.call <- list()
      store.sub.class <- lvl.df[[k]][!is.na(lvl.df[[k]])]
      parent.class <- unique(gsub("(.*)-(.*)", "\\1", store.sub.class))
      
      for (parent in parent.class) {
        sub.class <- store.sub.class[grepl(parent, store.sub.class)]
        children_class <- gsub("(.*)-(.*)", "\\2", sub.class)
        parents <- call.res$Barcode[call.res$MS_Celltype == 
                                      parent]
        trunc.meta.parent <- object@meta.data[parents, 
        ] %>% select(children_class)
        
        gap_ind <- which(names(figures) == parent)
        
        # Stop hierarchical classification in case no parent cell can be called
        if (nrow(trunc.meta.parent) == 0){
          stop(paste0("No ",parent," can be called in ","level ",k-1," classification, try setting more lenient thresholds"))}
        
        for (child in children_class) {
          plot.title <- paste("Density plot for", child, 
                              "Module Scores within", parent, "population", 
                              sep = " ")
          
          figures <- append(figures, list(NA), after = gap_ind)
          
          figures[[gap_ind+1]] <- ggplot(trunc.meta.parent, aes_string(x = child)) + geom_density() + ggtitle(plot.title) + geom_vline(xintercept = numeric_threshold[child], linetype = "dashed", color = "red3") + theme_classic()
          names(figures)[gap_ind+1] <- child
        }
        
        trunc.meta.no.parent <- call.res[!call.res$MS_Celltype == 
                                           parent, ]
        non.parent <- rownames(trunc.meta.no.parent)
        child.thres.vec <- numeric_threshold[children_class]
        
        sub.class.call[[match(parent, parent.class)]] <- .modScoreCall(trunc.meta.parent, child.thres.vec, reject = parent) %>% select(MS_Celltype)}
      
      sub.class.call <- do.call(rbind, sub.class.call)
      sub.class.call$Barcode <- rownames(sub.class.call)
      call.res$temp.call <- sub.class.call$MS_Celltype[match(call.res$Barcode, sub.class.call$Barcode)]
      call.res <- call.res %>% mutate(MS_Celltype = case_when(is.na(temp.call) ~ MS_Celltype, TRUE ~ temp.call))
      call.res$temp.call <- NULL
    }
  }
  
  object@meta.data$MS_Celltype <- call.res$MS_Celltype[match(object@meta.data$Barcode, call.res$Barcode)]
  
  lapply(figures, plot)
  
  return(object)
}
 