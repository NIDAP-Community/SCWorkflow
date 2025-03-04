#' @title Violin Plot by Metadata
#' @description Create violin plot of gene expression data across groups
#' @details Takes in a list of genes inputted by the user, displays violin plots
#'          of genes across groups from a slot-assay with (optional) outliers 
#'          removed. Can also choose to scale or transform expression data.
#' 
#' @param object Seurat-class object
#' @param assay Assay to extract gene expression data from (Default: SCT)
#' @param slot Slot to extract gene expression data from (Default: scale.data)
#' @param group.by Split violin plot based on metadata group
#' @param group.subset Include only a specific subset from group.by
#' @param genes.of.interest Genes to visualize on the violin plot
#' @param filter.outliers Filter outliers from the data (TRUE/FALSE)
#' @param scale.data Scale data from 0 to 1 (TRUE/FALSE)
#' @param log.scale.data Transform data onto a log10 scale (TRUE/FALSE)
#' @param reorder.ident Numeric data will be ordered naturally by default. 
#'                      Toggling this option will order the groups to match the
#'                      group list if non-numeric, and will have no effect if 
#'                      otherwise.
#' @param rename.ident Give alternative names to group.by displayed on 
#'                     the graph
#' @param ylimit Y-axis limit
#' @param plot.style Choose between grid, labeled, and row
#' @param outlier.low.lim Filter lower bound outliers (Default = 0.1)
#' @param outlier.up.lim Filter upper bound outliers (Default = 0.9)
#' @param jitter.points Scatter points on the plot (TRUE/FALSE)
#' @param jitter.width Set spread of jittered points 
#' @param jitter.dot.size Set size of individual points
#' @param print.outliers Print outliers as points in your graph that may be 
#'                       redundant to jitter 

#' @import Seurat 
#' @import reshape2
#' @import tidyverse
#' @import cowplot
#' @import rlang
#' @import ggplot2
#'   
#' @export
#' @example Do not run: violinPlot(object = seurat,
#'                                 group.by = "celltype",
#'                                 group.subset = c("CD4_Tcell","CD8_Tcell")
#'                                 genes.of.interest = c("Cd4","Cd8a"))

#' @return violin ggplot2 object

violinPlot <- function (object, assay, slot, genes, group, facet_by = "", jitter_points, jitter_dot_size) 
{
  library(Seurat)
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(broom)
  
  facet_data = facet_by != ""
  
  if (!assay %in% Assays(object)) {
    stop("expression data type was not found in Seurat object")
  } else if (!slot %in% slotNames(object[[assay]])) {
    stop("slot not found in Seurat[[assay]]")
  } else if (all(!genes %in% rownames(object[[assay]]))) {
    stop("no genes were found in Seurat object")
  } else if (!group %in% colnames(object@meta.data)) {
    stop("grouping parameter was not found in Seurat object")
  }
  
  # Scale to non-negative for visualization
  gene_mtx <- as.matrix(GetAssayData(object, assay = assay, slot = slot))
  gene_mtx <- scales::rescale(gene_mtx, to = c(0,1))
  
  print(paste0(genes[!genes %in% rownames(gene_mtx)], 
               " not found and will not be displayed"))
  
  genes.present <- genes[genes %in% rownames(gene_mtx)]
  if(facet_data){
    meta_sub <- object@meta.data[,c(group,facet_by)]
  } else {
    meta_sub <- object@meta.data[c(group)]
  }
  
  for (col in genes.present) {
    meta_sub[[col]] <- gene_mtx[col,]
  }
  
  data_df <- meta_sub %>% pivot_longer(genes.present, names_to = "Gene", values_to = "Expression")
  
  # set Gene as factor in data_df, so faceted plots will not be alphabetical 
  data_df$Gene <- factor(data_df$Gene, levels = genes.present)
  
  if(facet_data){
    unique_facets <- unique(object@meta.data[,facet_by])
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
  
  if(facet_data){
    # Map the colors to the unique values
    # If there are more unique sets than available colors, this will recycle the colors
    color_mapping <- setNames(rep(available_colors, length.out = length(unique_facets)), unique_facets)
    
    # Set up the common elements of the plot
    g <- ggplot(data_df, aes(x = .data[[group]], y = Expression, fill = .data[[facet_by]])) + 
      geom_violin(scale = "width", position = position_dodge(width = 0.9), trim = TRUE) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
      scale_fill_manual(values = color_mapping) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 3, strip.position = "left") + 
      theme_classic() + 
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 14, color = "black", face = "bold"),
            strip.text.y = element_text(size = 14, color = "black", face = "bold"),
            strip.placement = "outside")
  } else {
    g <- ggplot(data_df, aes(x = .data[[group]], y = Expression, fill = .data[[group]])) + 
      geom_violin(scale = "width", position = position_dodge(width = 0.9), trim = TRUE) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
      scale_fill_brewer(palette = "Set1") +
      theme_classic() + 
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 14, color = "black", face = "bold"),
            strip.text.y = element_text(size = 14, color = "black", face = "bold"),
            strip.placement = "outside")
  }
  
  # Add jitter points conditionally
  if (jitter_points) {
    g <- g + geom_jitter(size = jitter_dot_size, shape = 1, position = position_dodge(width = 0.9), alpha = 0.5)
  }
  
  # anova for subgroups if facet_data is on
  p_values_df <- NULL
  if(facet_data){
    # Function to calculate p-values for a single gene within a cell type
    calculate_p_values <- function(data, data_group, data_gene) {
      # Subset data for the specific cell type and gene
      data_sub <- data[data[,group] == data_group & data[,"Gene"] == data_gene,]
      
      # Perform ANOVA and Tukey HSD
      fit <- aov(as.formula(paste("Expression ~", facet_by)), data = data_sub)
      tukey_result <- TukeyHSD(fit)
      
      # Tidy up the results and add metadata
      tidy_tukey_result <- tidy(tukey_result)
      tidy_tukey_result$gene <- data_gene
      tidy_tukey_result$group <- data_group
      
      return(tidy_tukey_result)
    }
    
    # List unique cell types
    unique_groups <- unique(data_df[[group]])
    
    # Filter out 
    facet_df <- table(data_df[[group]], data_df[[facet_by]])
    
    # Find rows with more than one non-zero column
    count_non_zero <- function(row) {
      sum(row != 0)
    }
    non_zero_counts <- apply(facet_df, 1, count_non_zero)
    
    # Use rownames whose values are in more than 1 column 
    unique_groups <- names(non_zero_counts)[non_zero_counts > 1]
    
    # Calculate p-values for each cell type and gene
    p_values_list <- list()
    for (indv_group in unique_groups) {
      p_values_list[[indv_group]] <- do.call(rbind, lapply(genes.present, function(gene) calculate_p_values(data_df, data_group = indv_group, data_gene = gene)))
    }
    
    # Combine the results into a single data frame
    p_values_df <- do.call(rbind, p_values_list)
  }
  
  final_res <- list(fig = g, stat = p_values_df)
  
  return(final_res)
  
}
