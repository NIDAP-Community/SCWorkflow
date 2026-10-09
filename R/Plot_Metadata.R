#' @title Plotting (i.e. coloring with) different columns of your Metadata Table
#' @description For each column selected, this template will produce a plot
#' (UMAP/TSNE/PCA; your choice) using the data in that column to color the cells
#' @details This is a downstream template for the Single-cell RNA-seq workflow
#' (requires dataset where Filter/QC/SingleR annotations have been run first)
#'
#' @param object A combined Seurat Object with metadata to plot
#' @param samples.to.include Which samples you would like to include
#' @param metadata.to.plot The metadata columns from your Metadata table
#' you would like to plot
#' @param columns.to.summarize The columns you would like to summarize
#' @param summarization.cut.off Select the number of categories you want
#' to display, while marking all other cells as "other." Default is 5
#' @param reduction.type What kind of visualization you would like to use
#' to plot your cells and metadata (tsne, umap, pca). Default is tsne
#' @param use.cite.seq TRUE if you would like to plot Antibody clusters
#' from CITEseq instead of scRNA. Default is FALSE
#' @param show.labels Whether to add labels or not to your reduction map.
#' Default is FALSE
#' @param legend.text.size Customize the size of the legend text in your charts.
#' Default is 1
#' @param legend.position Select how you want to align your legend.
#' Default is "right"
#' @param dot.size The size of the dots displayed on the plot. Default os 0.01
#'
#'
#' @import Seurat
#' @import ggplot2
#' @import RColorBrewer
#' @import scales
#' @import tidyverse
#' @import ggrepel
#' @import gdata
#' @import reshape2
#' @import tools
#' @import grid
#' @import gridBase
#' @import gridExtra
#' @import dplyr
#'
#' @export
#'
#' @return a data.frame extracted from the Seurat object and plot
#'
#' @examples
#' \dontrun{
#' out <- plotMetadata(
#'   object = anno_so,
#'   samples.to.include = c("sample1", "sample2"),
#'   metadata.to.plot = c("celltype", "orig.ident"),
#'   columns.to.summarize = c(),
#'   reduction.type = "umap"
#' )
#' }

plotMetadata <- function(
                        #Basic Parameters:
                        object,
                        samples.to.include,
                        metadata.to.plot,
                        columns.to.summarize,
                        summarization.cut.off = 5,
                        reduction.type = "tsne",
                        use.cite.seq = FALSE,
                        show.labels = FALSE,
                        legend.text.size = 1,
                        legend.position = "right",
                        dot.size = 0.01
  ) {
  
  ###################
  ##   Functions   ##
  ###################
  .collapseForMessage <- function(values, max.items = 30) {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0) {
      return("none")
    }
    suffix <- if (length(values) > max.items) ", ..." else ""
    paste0(paste(head(values, max.items), collapse = ", "), suffix)
  }
  
  .drawMetadata <- function(m) {
    #check if there are NaNs in metadata, if there are, catch
    if (any(is.na(meta.df[[m]]))) {
      stop(sprintf(
        "Metadata column '%s' contains NA values and cannot be plotted. Non-NA values in this column include: %s. Available metadata columns without NA values: %s",
        m,
        .collapseForMessage(meta.df[[m]]),
        .collapseForMessage(valid.columns)
      ))
    }
    #Making a plot based on tsne/umap/pca and CiteSeq settings:
    reduction = reduction.type
    
    if (!(use.cite.seq)) {
      if (reduction == "tsne") {
        p1 <- DimPlot(object.sub,
                      reduction = "tsne",
                      group.by = "ident")
      } else if (reduction == "umap") {
        p1 <- DimPlot(object.sub,
                      reduction = "umap",
                      group.by = "ident")
      } else {
        p1 <- DimPlot(object.sub,
                      reduction = "pca",
                      group.by = "ident")
      }
    } else {
      if (reduction == "tsne") {
        p1 <-
          DimPlot(object.sub,
                  reduction = "protein_tsne",
                  group.by = "ident")
      } else if (reduction == "umap") {
        p1 <-
          DimPlot(object.sub,
                  reduction = "protein_umap",
                  group.by = "ident")
      } else {
        p1 <-
          DimPlot(object.sub,
                  reduction = "protein_pca",
                  group.by = "ident")
      }
    }
    
    #Categorical/Qualitative Variables
    if (!is.element(class(meta.df[[m]][1]), c("numeric", "integer"))) {
      if (!(use.cite.seq)) {
        #plot RNA clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        }
        
      } else {
        #else plot Antibody clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            clusid = as.character(object.sub@meta.data[[m]])
          )
        }
      }
      
      clusmat %>% group_by(clusid) %>% dplyr::summarise(umap1.mean = mean(umap1),
                                                 umap2.mean = mean(umap2)) -> umap.pos
      title = as.character(m)
      cols = list()
      
      lab.palette <- colorRampPalette(brewer.pal(12, "Paired"))
      n = length(unique((object.sub@meta.data[[m]])))
      cols[[1]] = brewer.pal(8, "Set3")[-2]  #Alternative
      cols[[2]] = brewer.pal(8, "Set1")
      cols[[3]] = c(cols[[1]], brewer.pal(8, "Set2")[3:6])
      cols[[4]] = c("#F8766D",
                    "#FF9912",
                    "#a100ff",
                    "#00BA38",
                    "#619CFF",
                    "#FF1493",
                    "#010407")
      cols[[5]] = c("blue", "red", "grey")
      cols[[6]] = lab.palette(n)
      cols[[7]] = c("red", "green", "blue", "orange", "cyan", "purple")
      cols[[8]] = c(
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
        "#a9a9a9",
        "#808080",
        "#A9A9A9",
        "#8B7355"
      )
      colnum = 8
      
      n = length(unique(clusmat$clusid))
      qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      qual.col.pals = qual.col.pals[c(7, 6, 2, 1, 8, 3, 4, 5),]
      col = unlist(mapply(
        brewer.pal,
        qual.col.pals$maxcolors,
        rownames(qual.col.pals)
      ))
      
      #Select to add labels to plot or not:
      if (show.labels) {
        g <- ggplot(clusmat) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          geom_point(aes(
            x = umap1,
            y = umap2,
            colour = clusid
          ), size = dot.size) +
          scale_color_manual(values = col) +
          xlab(paste(reduction, "-1")) + ylab(paste(reduction, "-2")) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            legend.text = element_text(size = rel(legend.text.size))
          ) +
          guides(colour = guide_legend(override.aes = list(
            size = 5, alpha = 1
          ))) +
          ggtitle(title) +
          geom_label_repel(
            data = umap.pos,
            aes(
              x = umap1.mean,
              y = umap2.mean,
              label = umap.pos$clusid
            ),
            size = 4
          )
      }
      else{
        g <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          geom_point(aes(colour = clusid), size = dot.size) +
          scale_color_manual(values = col) +
          xlab(paste(reduction, "-1")) + ylab(paste(reduction, "-2")) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            legend.text = element_text(size = rel(legend.text.size))
          ) +
          guides(colour = guide_legend(override.aes = list(
            size = 5, alpha = 1
          ))) +
          ggtitle(title)
      }
      
      
    } else {
      ##THIS IS THE PART WE PLOT QUANTITATIVE DATA
      m = as.character(m)
      clusid = object.sub@meta.data[[m]]
      clusid = scales::rescale(object.sub@meta.data[[m]], to = c(0, 1))
      clus.quant = quantile(clusid[clusid > 0], probs = c(.1, .5, .9))
      midpt.1 = clus.quant[2]
      midpt.2 = clus.quant[1]
      midpt.3 = clus.quant[3]
      #hist(clusid[!is.na(clusid)], breaks=100, main=m)
      #abline(v=midpt,col="red",lwd=2)
      
      if (!(use.cite.seq)) {
        #plot RNA clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        }
        
      } else {
        #else plot Antibody clusters
        if (reduction == "tsne") {
          clusmat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        } else if (reduction == "umap") {
          clusmat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        } else {
          clusmat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            clusid = as.numeric(object.sub@meta.data[[m]])
          )
        }
      }
      
      clusmat %>% group_by(clusid) %>% dplyr::summarise(umap1.mean = mean(umap1),
                                                 umap2.mean = mean(umap2)) -> umap.pos
      title = as.character(m)
      clusmat %>% dplyr::arrange(clusid) -> clusmat
      g <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        geom_point(aes(colour = clusid), size = 1) +
        #scale_color_gradient2(low = "blue4", mid = "white", high = "red",
        #          midpoint = midpt[[p]], na.value="grey",limits = c(0, 1)) +
        scale_color_gradientn(
          colours = c("blue4", "lightgrey", "red"),
          values = scales::rescale(c(0, midpt.2, midpt.1, midpt.3, 1),
                                   limits = c(0, 1))
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        ) +
        ggtitle(title) +
        xlab(paste0(reduction.type, "-1")) + ylab(paste0(reduction.type, "-2"))
    }
    
    return(g)
  }
  
    
  ###################
  ##   MAIN CODE   ##
  ###################    
    
    summarize.cut.off <- min(summarization.cut.off, 20)

    metadata.columns <- colnames(object@meta.data)
    sample.metadata.column <- if ("orig.ident" %in% metadata.columns) {
      "orig.ident"
    } else if ("orig_ident" %in% metadata.columns) {
      "orig_ident"
    } else {
      stop(
        "Found neither orig.ident nor orig_ident in object metadata. Please provide an object with one of these metadata column names."
      )
    }
    
    # checking for samples included:
    samples <- samples.to.include
    if (is.character(samples) && any(grepl('c\\(|\\[\\]', samples))) {
      samples <- eval(parse(text = gsub('\\[\\]', 'c()', samples)))
    }
    
    if (length(samples) == 0) {
      samples = unique(object@meta.data[[sample.metadata.column]])
    }
    
    if ("active.ident" %in% slotNames(object)) {
      sample_name = as.factor(object@meta.data[[sample.metadata.column]])
      names(sample_name) = names(object@active.ident)
      object@active.ident <- as.factor(vector())
      object@active.ident <- sample_name
      object.sub = subset(object, ident = samples)
    } else {
      sample_name = as.factor(object@meta.data[[sample.metadata.column]])
      names(sample_name) = names(object@active.ident)
      object@active.ident <- as.factor(vector())
      object@active.ident <- sample_name
      object.sub = subset(object, ident = samples)
    }

    meta.df <- object.sub@meta.data
    
    available.metadata.columns <- colnames(object.sub@meta.data)
    available.metadata.columns <-
      available.metadata.columns[!grepl("Barcode", available.metadata.columns)]
    available.metadata.columns.message <-
      if (length(available.metadata.columns) > 0) {
        paste(available.metadata.columns, collapse = ", ")
      } else {
        "none"
      }
    
    
    # checking metadata for sanity
    if (is.character(metadata.to.plot) && any(grepl('c\\(|\\[\\]', metadata.to.plot))) {
      m = eval(parse(text = gsub('\\[\\]', 'c()', metadata.to.plot)))
    }else{
      m=metadata.to.plot
    }
    m <- trimws(m)
    m <- m[!is.na(m) & nzchar(m)]
    
    m = m[!grepl("Barcode", m)]
    if (length(m) == 0) {
      stop(sprintf(
        "metadata.to.plot must include at least one metadata column. Available metadata columns: %s",
        available.metadata.columns.message
      ))
    }
    missing.metadata.columns <- setdiff(m, colnames(object.sub@meta.data))
    if (length(missing.metadata.columns) > 0) {
      stop(sprintf(
        "metadata.to.plot contains metadata columns not found in the object: %s. Available metadata columns: %s",
        paste(missing.metadata.columns, collapse = ", "),
        available.metadata.columns.message
      ))
    }
    
    #ERROR CATCHING
    #collect valid names of valid columns
    valid.columns <- character()
    for (i in colnames(meta.df)) {
      if (!any(is.na(meta.df[[i]]))) {
        valid.columns <- c(valid.columns, i)
      }
    }
    
    # Checking for content of "Columns to Summarize"
    cols.to.summarize <-
      eval(parse(text = gsub('\\[\\]', 'c()', columns.to.summarize)))
    cols.to.summarize <- trimws(cols.to.summarize)
    cols.to.summarize <- cols.to.summarize[!is.na(cols.to.summarize) & nzchar(cols.to.summarize)]
    m = unique(c(m, cols.to.summarize))
    
    if (length(cols.to.summarize) > 0) {
      #Optional Summarization of Metadata
      for (i in cols.to.summarize) {
        col <- meta.df[[i]]
        val.count <- length(unique(col))
        
        if ((val.count >= summarize.cut.off) &
            (i != 'Barcode') &
            (!is.element(class(meta.df[[i]][1]), c("numeric", "integer")))) {
          freq.vals <- as.data.frame(-sort(-table(col)))$col[1:summarize.cut.off]
          summarized.col = list()
          count <- 0
          for (j in col) {
            if (is.na(j) || is.null(j) || (j == "None")) {
              count <- count + 1
              summarized.col[count] <- "NULLorNA"
            } else if (j %in% freq.vals) {
              count <- count + 1
              summarized.col[count] <- j
            } else {
              count <- count + 1
              summarized.col[count] <- "Other"
            }
          }
          meta.df[[i]] <- summarized.col
        }
      }
      #assign new metadata
      object.sub@meta.data <- meta.df
    }
    
    
  
  grobs <- lapply(m, function(x) .drawMetadata(x))
  cat("Available metadata columns: ", available.metadata.columns.message, "\n", sep = "")
  

  result.list <- list(
    "object" = object, 
    "plots" = grobs
    )
  return(result.list)
  
}
