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
#' @importFrom Seurat GetAssayData Idents subset FetchData DefaultAssay
#' @importFrom reshape2 melt
#' @importFrom dplyr filter group_by mutate
#' @importFrom cowplot plot_grid
#' @importFrom rlang "%||%"
#' @importFrom ggplot2 ggplot aes_string geom_violin theme_classic labs theme 
#'  element_text ylim geom_jitter scale_y_log10 geom_boxplot 
#'  facet_grid facet_wrap
#'   
#' @export
#' @example Do not run: violinPlot(object = seurat,
#'                                 group.by = "celltype",
#'                                 group.subset = c("CD4_Tcell","CD8_Tcell")
#'                                 genes.of.interest = c("Cd4","Cd8a"))

#' @return violin ggplot2 object

violinPlot <- function(object, 
                       assay = "SCT", 
                       slot = "scale.data", 
                       group.by, 
                       group.subset = c(), 
                       genes.of.interest,
                       filter.outliers = FALSE, 
                       scale.data = TRUE, 
                       log.scale.data = FALSE, 
                       reorder.ident = TRUE,
                       rename.ident = "", 
                       ylimit = 0, 
                       plot.style = "grid", 
                       outlier.low.lim = 0.1, 
                       outlier.up.lim = 0.9, 
                       jitter.points = FALSE,
                       jitter.width = 0.05, 
                       jitter.dot.size = 0.5, 
                       print.outliers = TRUE){
  
  # Error Messages
  gene.filter <- genes.of.interest %in% rownames(GetAssayData(
    object = object, slot = slot, assay = assay))
  missing.genes <- genes.of.interest[!gene.filter]
  genes.of.interest <- genes.of.interest[gene.filter]
  
  if(length(missing.genes) > 0){
    print(paste("The following genes are missing from the dataset:", 
                missing.genes, sep = " "))
  }
  
  if(length(genes.of.interest) == 0){
    stop("No query genes were found in the dataset.")
  }
  
  if(!group.by %in% colnames(object@meta.data)){
    colnames(object@meta.data) <- gsub("\\.","_",colnames(object@meta.data))
    if(!group.by %in% colnames(object@meta.data)){
      stop("Unable to find ident of interest in metadata.")
    }
  }
  
  group.filter <- group.subset %in% object@meta.data[[group.by]]
  group.subset <- group.subset[group.filter]
  missing.groups <- group.subset[!group.filter]
  if(length(missing.groups) > 0){
    cat("The following groups are missing from the selected ident:\n")
    print(missing.groups)
  }
  
  if(rename.ident %in% c("Gene","Expression","scaled")){
    stop("New ident name cannot be one of Gene, Expression, or scaled.")
  }
  
  # Helper Functions
  # splitFacet helper function comes from https://stackoverflow.com/questions/30510898/split-facet-plot-into-list-of-plots/52225543
  .splitFacet <- function(x){
    facet_vars <- names(x$facet$params$facets)         # 1
    x$facet    <- ggplot2::ggplot()$facet              # 2
    datasets   <- split(x$data, x$data[facet_vars])    # 3
    new_plots  <- lapply(datasets,function(new_data) { # 4
      x$data <- new_data
      x})
  }
  
  # from rapportools
  .isEmpty <- function(x, trim = TRUE, ...) {
    if (length(x) <= 1) {
      if (is.null(x))
        return (TRUE)
      if (length(x) == 0)
        return (TRUE)
      if (is.na(x) || is.nan(x))
        return (TRUE)
      if (is.character(x) && nchar(ifelse(trim, .trimSpace(x), x)) == 0)
        return (TRUE)
      if (is.logical(x) && !isTRUE(x))
        return (TRUE)
      if (is.numeric(x) && x == 0)
        return (TRUE)
      return (FALSE)
    } else
      sapply(x, .isEmpty, trim = trim, ...)
  }
  
  # from rapportools
  .trimSpace <- function(x, what = c('both', 'leading', 'trailing', 'none'), 
                         space.regex = '[:space:]', ...){
    if (missing(x))
      stop('nothing to trim spaces to =(')
    re <- switch(match.arg(what),
                 both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
                 leading  = sprintf('^[%s]+', space.regex),
                 trailing = sprintf('[%s]+$', space.regex),
                 none     = {
                   return (x)
                 })
    .vgsub(re, '', x, ...)
  }
  
  # from rapportools
  .vgsub <- function(pattern, replacement, x, ...){
    for(i in 1:length(pattern))
      x <- gsub(pattern[i], replacement[i], x, ...)
    x
  }
  
  .removeOutliers <- function(x, na.rm = TRUE){
    qnt <- quantile(x, probs=c(outlier.low.lim,outlier.up.lim), na.rm = na.rm)
    H <- 1.5*IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  # Main Code Block
  # deal with limits
  if(ylimit == 0){
    ylimit <- NULL
  }
  
  Idents(object) <- object@meta.data[[group.by]]
  if(!is.null(group.subset)){
  object.sub <- subset(object, idents = group.subset)
  } else {
    object.sub <- object
  }
  
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object.sub, vars = genes.of.interest, slot = slot)
  
  append <- object.sub@meta.data[[group.by]]
  data[[group.by]] <- append[match(rownames(data),colnames(object.sub))]
  
  df.melt <- reshape2::melt(data)
  
  if(!.isEmpty(rename.ident)){
    group.by <- rename.ident
  }
  colnames(df.melt) <- c(group.by,"Gene","Expression")
  
  #check to see if ident of interest looks numeric
  if(suppressWarnings(all(!is.na(as.numeric(as.character(df.melt[[
    group.by]])))))){
    ident.values <- strtoi(df.melt[[group.by]])
    ident.levels <- unique(ident.values)[order(unique(ident.values))]
    df.melt[[group.by]] <- factor(ident.values, levels = ident.levels)
  }else if(reorder.ident){
    # if non-numeric, place in order of groups of interests
    if(!is.null(group.subset)){
    df.melt[[group.by]] <- factor(df.melt[[group.by]], 
                                  levels = group.subset)
    } else {
      df.melt[[group.by]] <- factor(df.melt[[group.by]])
    }        
  }
  
  # Filter outliers
  if(filter.outliers){
    for(gene in genes.of.interest){
      for(group in group.subset){
        current.ind <- which(df.melt[["Gene"]] == gene & df.melt[[
          group.by]] == group)
        df.melt[current.ind,"Expression"] <- .removeOutliers(df.melt[
          current.ind,"Expression", drop = TRUE])
      }
    }
  }
  
  # Scale data to y limit
  if(scale.data){
    expression.data = "scaled"
    axis.title.y = "Expression (scaled)"
    ylimit <- ylimit %||% 1
    df.melt <- df.melt %>% group_by(Gene) %>% mutate(
      scaled = scales::rescale(Expression, to=c(0,ylimit)))
  }else{
    expression.data <- axis.title.y <- "Expression"
  }
  
  g <- ggplot(df.melt, aes_string(x=group.by, y=expression.data)) +
    geom_violin(aes_string(fill = group.by), scale="width", 
                trim = FALSE, show.legend = FALSE) + 
    theme_classic() + 
    labs(y=axis.title.y) +
    theme(strip.text.y = element_text( 
      color="blue", face="bold.italic", angle = -90))
  
  if(!is.null(ylimit)){
    g <- g + ylim(0,ylimit)
  }
  
  if(jitter.points){
    g <- g + geom_jitter(height = 0, width = jitter.width, size=jitter.dot.size)
  }
  
  if(log.scale.data){
    g <- g + scale_y_log10()
  }
  
  # Plot after jitter if wanted
  g <- g + geom_boxplot(width=0.1, fill="white", outlier.shape = 
                          ifelse(print.outliers,19,NA))
  
  # Plot styles
  ncol = ceiling(length(unique(df.melt$Gene))^0.5)
  nrow = ceiling(length(unique(df.melt$Gene)) / ncol)
  if(plot.style == "rows"){
    g <- g + facet_grid(rows = vars(Gene))
  }else{
    g <- g + facet_wrap(~Gene, nrow = nrow, ncol = ncol)
    if(plot.style == "labeled"){
      plots <- .splitFacet(g)
      plots <- lapply(seq_along(plots), function(i) plots[[i]] + 
                        ggtitle(genes.of.interest[i]) + 
                        theme(plot.title = element_text(hjust = 0.5)) )
      g <- plot_grid(plotlist = plots, nrow = nrow, ncol = ncol, 
                     labels = LETTERS[seq( from = 1, to = length(plots) )])
    }
  }
  
  return(g)
}
