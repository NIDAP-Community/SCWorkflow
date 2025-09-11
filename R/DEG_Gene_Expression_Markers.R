#' @title DEG (Gene Expression Markers)
#' @description This function performs a DEG (differential expression of genes)
#' analysis on a merged Seurat object to identify expression markers
#' between different groups of cells (contrasts).
#' @details The recommended input is a merged Seurat object
#' with SingleR annotations, along with its associated sample names and metadata
#'
#' @param object Seurat-class object
#' @param samples Samples to be included in the analysis
#' @param contrasts Contrasts in the "A-B" format
#' @param parameter.to.test Select the metadata column that you would like
#' to use to perform your DEG analysis and construct your contrasts from.
#' Default is "orig_ident"
#' @param test.to.use The kind of algorithm you would like to use
#' to perform your DEG analysis. Default is the MAST algorithm
#' (wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2).
#' @param log.fc.threshold The minimum log fold-change between contrasts
#' that you would like to analyze. Default is 0.25
#' @param use.spark Opt to use Spark to parallelize computations.
#' Default is FALSE
#' @param assay.to.use The assay to use for your DEG analysis.
#' Default is SCT, but can use linearly scaled data by selecting RNA instead


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
#' @import parallel
#' @import MAST
#'
#'
#' @export
#'
#' @return a dataframe with DEG.



degGeneExpressionMarkers <- function(object,
                                     samples,
                                     contrasts,
                                     parameter.to.test = "orig_ident",
                                     test.to.use = "MAST",
                                     log.fc.threshold = 0.25,
                                     use.spark = FALSE,
                                     assay.to.use = "SCT"
                                     ) {
  ## --------------- ##
  ## Functions       ##
  ## --------------- ##
  
  #define and call function for running DEG
  .getDegTable <- function(n) {
    first.cluster <- unlist(n)[1]
    second.cluster <- unlist(n)[2]
    second.cluster.name <- second.cluster
    
    if (second.cluster == "all") {
      second.cluster <- NULL
    }
    
    Idents(object.sub) <- param.to.test
    
    markers = FindMarkers(
      object.sub,
      ident.1 = first.cluster,
      ident.2 = second.cluster,
      test.use = test.to.use,
      logfc.threshold = log.fc.threshold,
      verbose = FALSE,
      assay = assay.to.use,
      slot = "counts"
    )
    
    
    colnames(markers) <-
      chartr(
        old = " ",
        new = "_",
        paste(
          colnames(markers),
          first.cluster,
          "vs",
          second.cluster.name,
          sep = "_"
        )
      )
    
    return(markers)
  }
  
  
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # Getting metadata and checking sample names:
  metadata.table <- object@meta.data
  samples = eval(parse(text = gsub('\\[\\]', 'c()', samples)))
  
  if (length(samples) == 0) {
    samples = unique(object@meta.data$sample_name)
  }
  
  colnames(object@meta.data) <-
    gsub("orig_ident", "orig.ident", colnames(object@meta.data))
  if ("active.ident" %in% slotNames(object)) {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples)
  } else {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples)
  }
  
  print("selected samples:")
  print(object.sub)
  
  colnames(object.sub@meta.data) = gsub("\\.", "_",
                                        colnames(object.sub@meta.data))
  
  #define contrasts
  new.cont <- list()
  for (i in 1:length(contrasts)) {
    new.cont[[i]] <- c(paste(unlist(strsplit(contrasts[i], "-"))))
  }
  contrasts <- new.cont
  
  #ERROR CATCHING
  #collect valid names of valid columns
  valid.columns <- character()
  for (i in colnames(metadata.table)) {
    if (!any(is.na(metadata.table[[i]]))) {
      valid.columns <- c(valid.columns, i)
    }
  }
  
  param.to.test <- parameter.to.test
  
  if (param.to.test == "") {
    mcols = colnames(object.sub@meta.data)
    param.to.test <- mcols[grepl("RNA_snn", mcols)][[1]]
    print(paste("No parameter selected, defaulting to", param.to.test))
  }
  
  contrast.target <- object.sub@meta.data[[param.to.test]]
  contrast.type <- param.to.test
  contrast.counts = as.data.frame(table(contrast.target))
  valid.contrasts = subset(contrast.counts, Freq > 2)[[1]]
  
  #catch malformed contrasts
  for (i in contrasts) {
    if (!(i[[1]] %in% contrast.target)) {
      print(paste(
        i[[1]],
        "is not a valid contrast for contrast type:",
        contrast.type
      ))
      print("Please see below for an example of valid contrasts
             for your selected contrast type.")
      print(valid.contrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[2]] %in% contrast.target) & (i[[2]] != "all")) {
      print(paste(
        i[[2]],
        "is not a valid contrast for contrast type:",
        contrast.type
      ))
      print("Please see below for an example of valid contrasts
             for your selected contrast type.")
      print(valid.contrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (length(i) > 2) {
      print("Contrasts are as follows..")
      print(i)
      stop(
        "The console says there are too many inputs in your contrasts.
         A contrast should only contain Group1-Group2,
         but the console thinks you have inputed Group1-Group2-Group3"
      )
    } else if (!(i[[2]] %in% valid.contrasts) & (i[[2]] != "all")) {
      print(
        paste(
          i[[2]],
          "has two few values (less than 3 cells) to contrast against.
           Please see below for contrasts with enough cells:",
          valid.contrasts
        )
      )
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[1]] %in% valid.contrasts)) {
      print(
        paste(
          i[[1]],
          "has two few values (less than 3 cells) to contrast against.
           Please see below for contrasts with enough cells:",
          valid.contrasts
        )
      )
      stop("You have entered an invalid group to contrast against.")
    }
  }
  
  #print out contrast cell contrast.counts
  for (i in seq_along(contrasts)) {
    first.group <- contrasts[[i]][[1]]
    first.group.count <-
      subset(contrast.counts, contrast.target == first.group)$Freq
    if (contrasts[[i]][[2]] != "all") {
      second.group <- contrasts[[i]][[2]]
      second.group.count <-
        subset(contrast.counts, contrast.target == second.group)$Freq
      print(
        paste(
          "Contrast No.",
          i,
          "contrasts cluster",
          first.group,
          "with",
          first.group.count,
          "cells vs. cluster",
          second.group,
          "with",
          second.group.count,
          "cells."
        )
      )
    } else {
      second.group.count <- ncol(object.sub) - first.group.count
      print(
        paste(
          "Contrast No.",
          i,
          "contrasts cluster",
          first.group,
          "with",
          first.group.count,
          "cells vs. all other clusters, totalling",
          second.group.count,
          "cells."
        )
      )
    }
  }
  
  
  
  if (use.spark) {
    deg.tables <- spark.lapply(contrasts, .getDegTable)
  } else {
    deg.tables <- lapply(contrasts, .getDegTable)
  }
  
  for (i in seq_along(deg.tables)) {
    degtab <- deg.tables[[i]]
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
    print(paste0(
      "The number of upregulated genes at p<0.05 in contrast number ",
      i,
      " is:"
    ))
    print(pos[1])
    print(paste0(
      "The number of downregulated genes at p<0.05 in contrast number ",
      i,
      " is:"
    ))
    print(neg[1])
  }
  
  #Merge the deg tables together
  out.df <- NULL
  for (i in deg.tables) {
    if (is.null(out.df)) {
      out.df <- deg.tables[1]
      out.df <- as.data.frame(out.df)
    } else {
      out.df <- merge(out.df, i, by = "row.names", all = TRUE)
      rownames(out.df) <- out.df$Row.names #set the rownames
      out.df$Row.names <-
        NULL #drop the row.names columns which we no longer need
    }
  }
  
  out.df$Gene <- rownames(out.df)
  out.df$Row.names <- NULL
  out.df <- out.df %>% dplyr::select(Gene, everything())
  #  return(out.df)
  
  result.list <- list("df" = out.df)
  return(result.list)
  
}
