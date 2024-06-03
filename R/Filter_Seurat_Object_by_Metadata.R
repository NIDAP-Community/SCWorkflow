#' @title Filter Seurat Object by Metadata
#' @description Filter and subset your Seurat object based on metadata column
#' @details This is a downstream template that should be loaded after
#' Step 5 of the pipeline (SingleR Annotations on Seurat Object)
#'
#' @param object A dataset containing your SingleR annotated/merged seurat object
#' @param samples.to.include Select which samples to include
#' @param sample.name Sample Name Column
#' @param category.to.filter What kind of metadata you want to subset by.
#' This should be one column in your Metadata table
#' @param values.to.filter One or more values where you want to filter
#' @param keep.or.remove TRUE to filter down to selected values,
#' FALSE to filter out selected values. Default is TRUE
#' @param greater.less.than Decide if you want to keep cells above or below
#' the threshold. Default is "greater than"
#' @param seed Set seed for colors
#' @param cut.off The cut-off you want to use for your
#' greater than/less than filter. Default os 0.5
#' @param legend.position Select "none" if legend takes up too much
#' space in plot. Default is "top"
#' @param reduction What kind of clustering visualization you would like to use
#' for summary plot (umap, tsne, pca, protein_tsne, protein_umap, protein_pca).
#' Default is "umap"
#' @param plot.as.interactive.plot TRUE for interactive, FALSE for static
#' @param legend.symbol.size Your legend symbol size. Default is 2
#' @param colors User-selected colors
#' from palette of 62 unique colors from ColorBrewer.
#' @param dot.size Size of dots on TSNE/UMAP projection plot. Default is 0.1
#' @param number.of.legend.columns Default is 1.
#' If legend is too long, provide more legend columns
#' @param dot.size.highlighted.cells Dot size for cells in the after-filter plot
#' which have been highlighted. Default is 0.5
#' @param use.cite.seq.data TRUE if you would like to plot Antibody clusters
#' from CITEseq instead of scRNA.

#'
#' @import Seurat
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @import gridBase
#' @import cowplot
#' @import RColorBrewer
#' @import colorspace
#' @import tidyverse
#'
#'
#'
#' @export
#'
#' @return a subset Seurat object


filterSeuratObjectByMetadata <- function(object,
                                         samples.to.include,
                                         sample.name,
                                         category.to.filter,
                                         values.to.filter,
                                         keep.or.remove = TRUE,
                                         greater.less.than = "greater than",
                                         seed = 10,
                                         cut.off = 0.5,
                                         legend.position = "top",
                                         reduction = "umap",
                                         plot.as.interactive.plot = FALSE,
                                         legend.symbol.size = 2,
                                         colors = c(
                                           "aquamarine3",
                                           "salmon1",
                                           "lightskyblue3",
                                           "plum3",
                                           "darkolivegreen3",
                                           "goldenrod1",
                                           "burlywood2",
                                           "gray70",
                                           "firebrick2",
                                           "steelblue",
                                           "palegreen4",
                                           "orchid4",
                                           "darkorange1",
                                           "yellow",
                                           "sienna",
                                           "palevioletred1",
                                           "gray60",
                                           "cyan4",
                                           "darkorange3",
                                           "mediumpurple3",
                                           "violetred2",
                                           "olivedrab",
                                           "darkgoldenrod2",
                                           "darkgoldenrod",
                                           "gray40",
                                           "palegreen3",
                                           "thistle3",
                                           "khaki1",
                                           "deeppink2",
                                           "chocolate3",
                                           "paleturquoise3",
                                           "wheat1",
                                           "lightsteelblue",
                                           "salmon",
                                           "sandybrown",
                                           "darkolivegreen2",
                                           "thistle2",
                                           "gray85",
                                           "orchid3",
                                           "darkseagreen1",
                                           "lightgoldenrod1",
                                           "lightskyblue2",
                                           "dodgerblue3",
                                           "darkseagreen3",
                                           "forestgreen",
                                           "lightpink2",
                                           "mediumpurple4",
                                           "lightpink1",
                                           "thistle",
                                           "navajowhite",
                                           "lemonchiffon",
                                           "bisque2",
                                           "mistyrose",
                                           "gray95",
                                           "lightcyan3",
                                           "peachpuff2",
                                           "lightsteelblue2",
                                           "lightyellow2",
                                           "moccasin",
                                           "gray80",
                                           "antiquewhite2",
                                           "lightgrey"
                                         ),
                                         dot.size = 0.1,
                                         number.of.legend.columns = 1,
                                         dot.size.highlighted.cells = 0.5,
                                         use.cite.seq.data = FALSE) {
  ## --------------- ##
  ## Parameters      ##
  ## --------------- ##
  
  
  ## --------------- ##
  ## Functions       ##
  ## --------------- ##
  
  # Drawing TSNE/UMAP/PCA plot
  .drawtsne <- function(SO, reduction, scale.col, col.grad) {
    SO.clus <- SO@meta.data[[category.to.filter]]
    
    plot1 <- DimPlot(SO, reduction = reduction, group.by = "ident")
    class(plot1$data$ident) <- "numeric"
    
    if (reduction == "tsne") {
      clusmat = data.frame(
        umap1 = plot1$data$tSNE_1,
        umap2 = plot1$data$tSNE_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    } else if (reduction == "umap") {
      clusmat = data.frame(
        umap1 = plot1$data$UMAP_1,
        umap2 = plot1$data$UMAP_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    } else if (reduction == "pca") {
      clusmat = data.frame(
        umap1 = plot1$data$PC_1,
        umap2 = plot1$data$PC_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    } else if (reduction == "protein_tsne") {
      clusmat = data.frame(
        umap1 = plot1$data$protein_tsne_1,
        umap2 = plot1$data$protein_tsne_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    } else if (reduction == "protein_umap") {
      clusmat = data.frame(
        umap1 = plot1$data$protein_umap_1,
        umap2 = plot1$data$protein_umap_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    } else {
      clusmat = data.frame(
        umap1 = plot1$data$protein_pca_1,
        umap2 = plot1$data$protein_pca_2,
        clusid = as.numeric(SO@meta.data[[category.to.filter]])
      )
    }
    
    # Preparing 
    clusmat %>% group_by(clusid) %>% summarise(umap1.mean = mean(umap1),
                                               umap2.mean = mean(umap2)) -> umap.pos
    title = as.character(category.to.filter)
    clusmat %>% dplyr::arrange(clusid) -> clusmat
    
    plot2 <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
      theme_bw() +
      theme(legend.title = element_blank()) +
      geom_point(
        aes(colour = clusid),
        alpha = 0.5,
        shape = 20,
        size = dot.size
      ) +
      scale_color_gradientn(colors = brewer.pal(n = 5, name = col.grad),
                            values = scale.col) +
      guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      ggtitle(title) +
      xlab("umap-1") + ylab("umap-2")
    return(plot2)
  }
  
  .distinctColorPalette <- function(k = 1, seed) {
    current.color.space <- our.color.space@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(current.color.space, k, iter.max = 20)
    colors <- unname(hex(LAB(km$centers)))
    return(colors)
  }
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # Checking if samples are selected
  samples = eval(parse(text = gsub('\\[\\]', 'c()', samples.to.include)))
  
  if (length(samples) == 0) {
    samples = unique(object@meta.data[[sample.name[1]]])
  }
  
  ## Replace dots in metadata column names with underscores.
  colnames(object@meta.data) = gsub("\\.", "_", colnames(object@meta.data))
  new.sample.name <- gsub("\\.", "_", sample.name[1])
  
  ## If you have protien data, then ...
  if (use.cite.seq.data) {
    reduction = paste("protein_", reduction, sep = '')
  }
  
  
  ## Original color-picking code.
  rand.col <- 2e3
  set.seed(seed)
  our.color.space <- colorspace::RGB(runif(rand.col), runif(rand.col), runif(rand.col))
  our.color.space <- as(our.color.space, "LAB")
  
  
  ## User-selected metadata column is used to set idents.
  Filter.orig = object@meta.data[[category.to.filter[1]]]
  colname <- category.to.filter[1]
  
  ident.of.interest = as.factor(object@meta.data[[colname]])
  names(ident.of.interest) = names(object@active.ident)
  object@active.ident <- as.factor(vector())
  object@active.ident <- ident.of.interest
  
  ## Get colors from user parameter and add more if the default list is too short.
  if (class(object@meta.data[[category.to.filter[1]]]) != "numeric") {
    col.length = length(levels(as.factor(Filter.orig)))
    if (length(colors) < col.length) {
      more.cols = .distinctColorPalette(col.length - length(colors), 10)
      colors <- c(colors, more.cols)
    }
    names(colors) <- levels(as.factor(Filter.orig))
    
    ## Keep or remove cells based on user input values.
    if (keep.or.remove) {
      subset.value <- values.to.filter
      meta.col <- unique(object@meta.data[[category.to.filter[1]]])
      print("Missing values:")
      print(setdiff(subset.value, meta.col))
      subset.value <- intersect(meta.col, subset.value)
    } else {
      meta.col <- unique(object@meta.data[[colname]])
      vals.to.remove <- values.to.filter
      subset.value <- setdiff(meta.col, vals.to.remove)
    }
    
    ## Subset Seurat object.
    SO.sub <- subset(object, idents = subset.value)
    
    ## Log output of tables of cell types by samples before and after filtes.
    print("Breakdown of filtered data:")
    print(table(object@meta.data[[category.to.filter[1]]],
                object@meta.data[[new.sample.name]]))
    
    cat("\n")
    print("After Filtering:")
    print(table(SO.sub@meta.data[[category.to.filter[1]]],
                SO.sub@meta.data[[new.sample.name]]))
    
    ## Set filter for the subsetted SO.
    SO.sub@meta.data[[colname]] <-
      as.factor(as.character(SO.sub@meta.data[[colname]])) #Relevel Factors
    
    filter.sub = SO.sub@meta.data[[colname]]
    
    #Set colors for unfiltered and filtered data by sample name:
    filt.length = length(levels(as.factor(filter.sub)))
    idx = vector("list", filt.length)
    names(idx) <- levels(as.factor(filter.sub))
    for (i in 1:filt.length) {
      id = Filter.orig %in% levels(as.factor(filter.sub))[i]
      idx[[i]] <- rownames(object@meta.data)[id]
    }
    cols2 <- colors[levels(as.factor(filter.sub))]
    
    ## Make before and after plots.
    title <-
      paste0("filtered by ",
             category.to.filter[1]#,
             ##             " and split by ",
             ##             category.to.filter[2]
      )
    plot1 = DimPlot(
      object,
      reduction = reduction,
      group.by = colname,
      pt.size = dot.size,
      raster=FALSE
    ) +
      theme_classic() +
      scale_color_manual(values = colors) +
      theme(legend.position = legend.position) +
      guides(colour = guide_legend(
        ncol = number.of.legend.columns,
        override.aes = list(size = legend.symbol.size)
      )) +
      ggtitle(colname)
    plot2 = DimPlot(
      object,
      reduction = reduction,
      cells.highlight = idx,
      cols.highlight = rev(cols2[1:filt.length]),
      sizes.highlight = dot.size.highlighted.cells,
      raster=FALSE
    ) +
      theme_classic() +
      theme(legend.position = legend.position) +
      guides(colour = guide_legend(
        ncol = number.of.legend.columns,
        reverse = TRUE,
        override.aes = list(size = legend.symbol.size)
      )) +
      ggtitle(title)
    
    ## Else, filter on numeric data with a user defined threshold and direction.
  } else {
    filter.direction <- greater.less.than
    meta.col <- unique(object@meta.data[[category.to.filter]])
    value <- cut.off
    if (filter.direction == "greater than") {
      SO.sub <- subset(object, subset = category.to.filter > cut.off)
    } else {
      SO.sub <- subset(object, subset = category.to.filter < cut.off)
    }
    
    
    clusid = object@meta.data[[category.to.filter]]
    maxclus = max(clusid)
    clusmid = 0.01 / maxclus
    min = min(clusid)
    midpt.1 = 0.99 * value
    midpt.0 = value
    midpt.2 = 1.01 * value
    max = max(clusid)
    col.points <- c(min, midpt.1, midpt.0, midpt.2, max)
    col.points <- scales::rescale(col.points, c(0, 1))
    
    
    
    plot1 <- .drawtsne(object, reduction, col.points, "RdBu")
    
    clusid = scales::rescale(SO.sub@meta.data[[category.to.filter]], to = c(0, 1))
    clus.quant = quantile(clusid[clusid > 0], probs = c(0, .25, .5, .75, 1))
    min = clus.quant[1]
    midpt.1 = clus.quant[3]
    midpt.3 = clus.quant[2]
    midpt.4 = clus.quant[4]
    max = clus.quant[5]
    col.points.2 <- c(min, midpt.3, midpt.1, midpt.4, max)
    
    plot2 <- .drawtsne(SO.sub, reduction, col.points.2, "Blues")
    
  }
  
  
  
  result.list <- list("object" = SO.sub,
                      "plot1" = plot1,
                      "plot2" = plot2)
  return(result.list)
}
