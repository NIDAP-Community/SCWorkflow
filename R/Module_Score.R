#' @title Compute ModScore
#' @description Returns Seurat-class object with metadata containing 
#'              ModuleScores and Likely_CellType calls
#' @details Analyzed features are binned based on averaged expression; 
#'          control features are randomly selected from each bin.
#' 
#' @param object Seurat-class object
#' @param samples.subset List of samples to subset the data by
#' @param sample.to.display List of samples to depict on dimension plot, 
#'                          samples not in the list would be colored gray in 
#'                          the background
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
#' @importFrom Seurat AddModuleScore DimPlot subset
#' @importFrom gridExtra arrangeGrob
#' @importFrom stringr str_detect
#' @importFrom quantmod paste0
#' @importFrom purrr map compact
#' @importFrom grid gpar textGrob
#' @importFrom data.table as.factor
#' @importFrom dplyr select group_by summarise mutate arrange
#'   
#' @export
#' @example Do not run: moduleScore(object = seurat,
#'                                  samples.subset = c("mouse1","mouse2"),
#'                                  sample.to.display = c("mouse1"),
#'                                  marker.table = immuneCellMarkers,
#'                                  celltypes = c("CD4_T","Treg",Monocytes"),
#'                                  threshold = c(0.1,0.4, 0.3),
#'                                  multi.lvl = FALSE
#'                                  )
#'                                  
#' @example Do not run: moduleScore(object = seurat,
#'                                  samples.subset = c("mouse1","mouse2"),
#'                                  sample.to.display = c("mouse1"),
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

modScore <- function(object, 
                     samples.subset, 
                     sample.to.display, 
                     marker.table, 
                     cite.seq = FALSE, 
                     celltypes, 
                     threshold = c(0), 
                     general.class, 
                     multi.lvl = FALSE, 
                     lvl.df, 
                     reduction = "tsne", 
                     nbins = 10, 
                     gradient.ft.size = 6, 
                     violin.ft.size = 6, 
                     step.size = 0.1){
  
  # Helper Functions
  # Give each cell an identity based on modulescores and bimodal threshold
  .modScoreCall <- function(ms.meta,threshold,reject){
    
    thres.ls <- lapply(1:ncol(ms.meta), 
                       function(marker_name) rep(threshold[marker_name],nrow(ms.meta)))
    thres.df <- data.frame(matrix(unlist(thres.ls),nrow = nrow(ms.meta)))
    thres.filter <- ms.meta > thres.df
    ms.meta.filt <- ms.meta * thres.filter
    max.col.vec <- max.col(ms.meta.filt)
    zero.filt <- as.integer(!apply(ms.meta.filt, 1, 
                                   function(find_zero_rows) 
                                     all(find_zero_rows == 0)))
    final.filt <- (max.col.vec * zero.filt) + 1
    append.name <- c(reject, names(ms.meta))
    ms.meta[,"Likely_CellType"] <- append.name[final.filt]
    return(ms.meta)
    
  }
  
  # Dimension reduction
  dim_reduction <- function(reduction, object.sub) {
    reduction <- match.arg(reduction, c("tsne", "umap", "pca"))
    DimPlot(object.sub, reduction = reduction, group.by = "ident")
  }
  
  # Prepare data for ggplots
  clus_mat_func <- function(reduction, object.sub, m) {
    clusid = as.numeric(object.sub@meta.data[[m]])
    cols <- c("UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2", "PC_1", "PC_2")
    cols <- cols[startsWith(cols, toupper(reduction))]
    clusmat <- data.frame(ident = object.sub@meta.data$orig.ident, 
                          reduct1 = object.sub[[reduction]]@cell.embeddings[,1],
                          reduct2 = object.sub[[reduction]]@cell.embeddings[,2], 
                          clusid = clusid)
    return(clusmat)
  }
  
  # Main Code Block
  marker.tab <- unlist(marker.table)
  if (!"Barcode" %in% colnames(object@meta.data)) {
    object@meta.data$Barcode <- rownames(object@meta.data)
  }
  
  if (length(samples.subset) == 0) {
    samples.subset = unique(object@meta.data$sample.name)
  }
  
  colnames(object@meta.data) <- gsub("orig_ident","orig.ident", 
                                     colnames(object@meta.data))
  
  Idents(object) <- object$orig.ident
  
  object.sub = subset(object, ident = samples.subset)
  rm(object)
  
  # Retrieve markers from the table and convert to a list
  marker.list <- dplyr::select(marker.table, celltypes) %>% 
    as.list()
  
  # Error checking and messages
  stopifnot(sum(unlist(marker.list) %in% rownames(object.sub@assays$SCT@data)) > 0,
            length(threshold) == length(celltypes) | sum(threshold) == 0)
  
  # Set threshold and names
  threshold <- if(sum(threshold) == 0) rep(0, length(celltypes)) else threshold
  names(threshold) <- celltypes
  
  # Determine present genes for each cell type
  gene.presence <- purrr::map(marker.list, ~ .x %in% rownames(object.sub))
  
  # Iterates over marker.list and gene.presence simultaneously
  marker.list <- purrr::map2(marker.list, gene.presence, function(markers, presence) {
    
    print(paste0("Present: ", paste(markers[presence], collapse=" ")))
    print(paste0("Not Present: ", paste(markers[!presence & markers != ""], collapse=" ")))
    
    if(sum(presence) == 0){
      print("genes were not found in object and will not be analyzed")
    }
    
    return(markers[presence])
    
  })
  
  figures <- list()
  
  # Calculate modulescores for each celltype in marker table
  invisible(
    for (i in seq_along(marker.list)) { 
      object.sub=AddModuleScore(object.sub,
                                marker.list[i],
                                name = names(marker.list[i]),
                                nbin = nbins)
      
      m = paste0(names(marker.list[i]),"1")
      object.sub@meta.data[[m]] <- scales::rescale(object.sub@meta.data[[m]], 
                                                   to=c(0,1))
      
      p1 <- dim_reduction(reduction, object.sub)
      
      clusmat <- clus_mat_func(reduction, object.sub, m)
      
      title=as.character(m)
      clusmat <- clusmat %>% dplyr::arrange(clusid) %>% mutate(
        sample_clusid = clusmat$clusid * grepl(
          paste(sample.to.display, collapse = "|"), 
          clusmat$ident))
      
      # Calculate density data for ModScore vs Number of cells
      d <- density(clusmat$clusid)
      
      # Plot data onto colored dimension reduction plot
      clusid.df <- data.frame(id=object.sub@meta.data$orig.ident,
                              ModuleScore=object.sub@meta.data[[m]])
      
      g <- ggplot(clusmat, aes(x = reduct1, y = reduct2)) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        geom_point(aes(colour = sample_clusid),alpha=0.5,shape = 20,size=1) +
        scale_color_gradientn(colours = c("blue4","lightgrey", "red"), 
                              values = scales::rescale(
                                c(0,threshold[i]/2,threshold[i],
                                  (threshold[i]+1)/2,1), limits = c(0, 1))) + 
        guides(colour = guide_legend(
          override.aes = list(size=5, alpha = 1))) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank()) +
        xlab("tsne-1") + ylab("tsne-2")
      
      g1 <- RidgePlot(object.sub,features=m,group.by="orig.ident") + 
        theme(legend.position = "none", title = element_blank(), 
              axis.text.x = element_text(size = gradient.ft.size)) + 
        geom_vline(xintercept = threshold[i], linetype = "dashed", 
                   color = "red3") +
        scale_x_continuous(breaks = seq(0,1,step.size))
      
      # Violin Plot
      g2 <- ggplot(clusid.df,aes(x=id,y=ModuleScore)) + 
        geom_violin(aes(fill=id)) + 
        theme_classic() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.title = element_blank(), 
              panel.background = element_blank(), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=rel(0.8)),
              legend.position="top", axis.text.y = 
                element_text(size = violin.ft.size)) + 
        guides(colour = guide_legend(
          override.aes = list(size=5, alpha = 1))) +
        geom_hline(yintercept = threshold[i], 
                   linetype = "dashed", color = "red3") +
        scale_y_continuous(breaks = seq(0,1,step.size))
      
      # Color gradient density plot
      g3 <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + 
        xlab("ModuleScore") + 
        ylab("Density") + 
        geom_line() + 
        geom_segment(aes(xend = d$x, yend = 0, colour = x)) + 
        scale_y_log10() +
        scale_color_gradientn(colours = c("blue4","lightgrey", "red"), 
                              values = scales::rescale(
                                c(0,threshold[i]/2,threshold[i],
                                  (threshold[i]+1)/2,1), limits = c(0, 1))) + 
        geom_vline(xintercept = threshold[i], 
                   linetype = "dashed", color = "red3") + 
        geom_vline(xintercept = threshold[i], 
                   linetype = "dashed", color = "red3") + 
        scale_x_continuous(breaks = seq(0,1,step.size)) + 
        theme(legend.title = element_blank(), 
              axis.text.x = element_text(size = 6))
      
      # Set title for arrangeGrob final figure
      figures[[i]] = arrangeGrob(g,g1,g2,g3, 
                                 ncol=2, 
                                 top=textGrob(names(marker.list[i]), 
                                              gp = gpar(fontsize = 14, 
                                                        fontface = "bold")))
    }
  )
  
  # Get rid of "1" at the end of MS columns
  colnames(object.sub@meta.data)[colnames(object.sub@meta.data) %in% paste0(
    names(marker.list),1)] <- names(marker.list)
  
  # Analyze and plot only the class of cells found in the metadata
  general.class <- general.class[general.class %in% 
                                   colnames(object.sub@meta.data)]
  
  # Subset the columns of the metadata containing module scores only
  trunc.meta.gen <- object.sub@meta.data[general.class]
  
  gen.thrs.vec <- threshold[general.class]
  
  # Set elements below threshold to zero. Keep elements above threshold
  call.res <- .modScoreCall(trunc.meta.gen,gen.thrs.vec,reject = "unknown")
  call.res$Barcode <- rownames(call.res)
  
  # Hierarchical Classification
  if (multi.lvl){   
    
    for (k in 1:ncol(lvl.df)){ 
      
      # Temporarily keep results from subpopulation calls
      sub.class.call <- list()
      
      ## Subclass Identification
      # Remove any NAs from comparisons
      store.sub.class <- lvl.df[[k]][!is.na(lvl.df[[k]])]
      
      parent.class <- unique(gsub("(.*)-(.*)","\\1",store.sub.class))
      
      for (parent in parent.class){
        sub.class <- store.sub.class[grepl(parent,store.sub.class)]
        children_class <- gsub("(.*)-(.*)","\\2",sub.class)
        
        # Subset out cells predicted to be parents
        parents <- call.res$Barcode[call.res$Likely_CellType == parent]
        trunc.meta.parent <- object.sub@meta.data[parents,] %>% 
          select(children_class)
        
        # Stores a new density plot containing MS information of child cells 
        for (child in children_class){
          
          plot.title <- paste("Density plot for",
                              child,"Module Scores within", 
                              parent,"population", sep = " ")
          adj.density <- ggplot(trunc.meta.parent, aes_string(
            x = child)) + 
            geom_density() + 
            ggtitle(plot.title) + 
            geom_vline(xintercept = threshold[child], 
                       linetype = "dashed", color = "red3")
          
          figures[[length(figures) + 1]] <- adj.density
        }
        
        # Create output table without information of parent population
        trunc.meta.no.parent <- call.res[!call.res$Likely_CellType == parent,]
        non.parent <- rownames(trunc.meta.no.parent)
        
        # Repeat Module Score Comparison and Cell Prediction with Child Subset:
        child.thres.vec <- threshold[children_class]
        
        sub.class.call[[match(parent,parent.class)]] <- .modScoreCall(
          trunc.meta.parent,child.thres.vec,reject = parent) %>% 
          select(Likely_CellType)
      }
      
      # Reappend subclassification results back to general output 
      sub.class.call <- do.call(rbind,sub.class.call)
      sub.class.call$Barcode <- rownames(sub.class.call)
      
      # Update Likely_CellType column in Calls Output 
      call.res$temp.call <- sub.class.call$Likely_CellType[
        match(call.res$Barcode,sub.class.call$Barcode)]
      
      call.res <- call.res %>% mutate(Likely_CellType = case_when(
        is.na(temp.call) ~ Likely_CellType,
        TRUE ~ temp.call
      ))
      
      # Remove temp.call
      call.res$temp.call <- NULL
    }}
  
  # Updating CellType(s) in metadata with subclass calls
  object.sub$Likely_CellType <- call.res$Likely_CellType[
    match(object.sub$Barcode,call.res$Barcode)]
  
  ms.res <- list(ms.figures = figures,
                 ms.object = object.sub)
  
  return(ms.res)
}
