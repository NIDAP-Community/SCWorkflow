#' @title Violin Plot by Metadata
#' @description Create violin plot of gene expression data across groups
#' @details Takes in a list of genes inputted by the user, displays violin plots
#'          of genes across groups from a layer-assay with (optional) outliers 
#'          removed. Can also choose to scale or transform expression data.
#' 
#' @param object Seurat-class object
#' @param assay Assay to extract gene expression data from (Default: SCT)
#' @param layer Slot to extract gene expression data from (Default: scale.data)
#' @param genes Genes to visualize on the violin plot
#' @param group Split violin plot based on metadata group
#' @param facet.by Split violin plot based on a second metadata group
#' @param filter.outliers Filter outliers from the data (TRUE/FALSE)
#' @param outlier.low Filter lower bound outliers (Default = 0.05)
#' @param outlier.high Filter upper bound outliers (Default = 0.95)
#' @param jitter.points Scatter points on the plot (TRUE/FALSE)
#' @param jitter.dot.size Set size of individual points
#' 
#' @import Seurat 
#' @import reshape2
#' @import tidyverse
#' @import tidyr pivot_longer
#' @import cowplot
#' @import rlang
#' @import ggplot2
#'   
#' @export
#' @examples
#' \dontrun{
#' violinPlot(
#'   object = seurat,
#'   assay = "SCT",
#'   layer = "data",
#'   genes = c("Cd4", "Cd8a"),
#'   group = "celltype",
#'   facet.by = "orig.ident",
#'   filter.outliers = TRUE,
#'   jitter.points = TRUE,
#'   jitter.dot.size = 0.5
#' )
#' }

#' @return violin ggplot2 object

violinPlot <- function (object, 
                            assay, 
                            layer, 
                            genes, 
                            group, 
                            facet.by = "", 
                            filter.outliers = F,
                            outlier.low = 0.05,
                            outlier.high = 0.95,
                            jitter.points, 
                            jitter.dot.size) 
{

  facet_data = facet.by != ""
  
  # for handling orig ident
  if (group == "orig.ident" | group == "orig_ident"){
    
    # grab orig ident, however it is presented in seurat metadata
    group <- colnames(object@meta.data)[grepl("origident",gsub('\\W|_',"",colnames(object@meta.data)))]
  }
  
  available_layers <- if(packageVersion("Seurat") >= "5.0.0"){
    Layers(object[[assay]])
  } else (
    slotNames(object[[assay]])
  )

  if (!assay %in% Assays(object)) {
    stop("expression data type was not found in Seurat object")
  } else if (!layer %in% available_layers) {
    stop("layer not found in Seurat[[assay]] Layers")
  } else if (all(!genes %in% rownames(object[[assay]]))) {
    stop("no genes were found in Seurat object")
  } else if (!group %in% colnames(object@meta.data)) {
    stop("grouping parameter was not found in Seurat object")
  }
  
  # Scale to non-negative for visualization
  gene_mtx <- as.matrix(GetAssayData(object, assay = assay, layer = layer))
  gene_mtx <- scales::rescale(gene_mtx, to = c(0,1))
  
  if(length(genes[!genes %in% rownames(gene_mtx)]) > 0){
    print(paste0(genes[!genes %in% rownames(gene_mtx)], 
                 " not found and will not be displayed"))
  }
  
  genes.present <- genes[genes %in% rownames(gene_mtx)]
  if(facet_data){
    meta_sub <- object@meta.data[,c(group,facet.by)]
  } else {
    meta_sub <- object@meta.data[c(group)]
  }
  
  for (col in genes.present) {
    meta_sub[[col]] <- gene_mtx[col,][match(rownames(meta_sub),names(gene_mtx[col,]))]
  }
  
  data_df <- meta_sub %>% pivot_longer(genes.present, names_to = "Gene", values_to = "Expression")
  
  # set Gene as factor in data_df, so faceted plots will not be alphabetical 
  data_df$Gene <- factor(data_df$Gene, levels = genes.present)
  
  if(facet_data){
    unique_facets <- unique(object@meta.data[,facet.by])
  } else{
    unique_facets <- NULL
  }
  available_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  
  available_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                        "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
                        "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", 
                        "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                        "#76B7B2", "#FF9D9A", "#B07AA1", "#D4A518", 
                        "#DE77AE", "#77AADD", "#EE8866", "#E4CDA7")
  
  # Define the outlier removal function
  .removeOutliers <- function(x, na.rm = TRUE){
    qnt <- quantile(x, probs = c(outlier.low, outlier.high), na.rm = na.rm)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    x[x < (qnt[1] - H) | x > (qnt[2] + H)] <- NA
    x
  }
  
  # Apply only if filtering is enabled
  if (filter.outliers) {
    group.vars <- colnames(data_df)[colnames(data_df) != 'Expression']
    
    data_df <- data_df %>%
      group_by(across(all_of(group.vars))) %>%
      mutate(Expression = .removeOutliers(Expression)) %>%
      ungroup()
  }
  
  if(facet_data){
    # Map the colors to the unique values
    # If there are more unique sets than available colors, this will recycle the colors
    color_mapping <- setNames(rep(available_colors, length.out = length(unique_facets)), unique_facets)
    
    # Set up the common elements of the plot
    g <- ggplot(data_df, aes(x = .data[[group]], y = Expression, fill = .data[[facet.by]])) + 
      geom_violin(scale = "width", position = position_dodge(width = 0.9), trim = TRUE) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
      scale_fill_manual(values = color_mapping) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 1, strip.position = "left") + 
      theme_classic() + 
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text.x = element_text(size = 14, color = "black", face = "bold"),
            strip.text.y = element_text(size = 14, color = "black", face = "bold"),
            strip.placement = "outside")
  } else {
    g <- ggplot(data_df, aes(x = .data[[group]], y = Expression, fill = .data[[group]])) + 
      geom_violin(scale = "width", position = position_dodge(width = 0.9), trim = TRUE) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 1, strip.position = "left") +
      scale_fill_brewer(palette = "Set1") +
      theme_classic() + 
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text.x = element_text(size = 14, color = "black", face = "bold"),
            strip.text.y = element_text(size = 14, color = "black", face = "bold"),
            strip.placement = "outside")
  }
  
  # Add jitter points conditionally
  if (jitter.points) {
    g <- g + geom_jitter(size = jitter.dot.size, shape = 1, position = position_dodge(width = 0.9), alpha = 0.5)
  }
  
  return(list("plots"=g))
}
