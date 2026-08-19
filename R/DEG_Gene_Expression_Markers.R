#' @title DE with Find Markers [CCBR] [scRNA-seq]
#' @description This function performs DE (differential expression) analysis on
#' a merged Seurat object to identify expression markers between different
#' groups of cells (contrasts).
#' @details The recommended input is a merged Seurat object
#' with SingleR annotations, along with its associated sample names and metadata
#'
#' @param object Seurat-class object
#' @param samples Samples to be included in the analysis. Leave blank to
#' include all samples.
#' @param contrasts Contrasts in the "A-B" format
#' @param parameter.to.test Select the metadata column that you would like
#' to use to perform your DEG analysis and construct your contrasts from.
#' Leave blank to stop and print categorical/discrete metadata columns.
#' @param test.to.use The kind of algorithm you would like to use
#' to perform your DEG analysis. Default is the MAST algorithm
#' (wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2).
#' @param log.fc.threshold The minimum log fold-change between contrasts
#' that you would like to analyze. Default is 0.25
#' @param use.spark Opt to use Spark to parallelize computations.
#' Default is FALSE
#' @param assay.to.use The assay to use for your DEG analysis.
#' Default is SCT, but can use linearly scaled data by selecting RNA instead
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
#' @import parallel
#' @import MAST
#'
#'
#' @export
#'
#' @return a dataframe with DEG.
#'
#' @examples
#' \dontrun{
#' deg <- degGeneExpressionMarkers(
#'   object = anno_so,
#'   samples = c("sample1", "sample2"),
#'   contrasts = c("A-B"),
#'   parameter.to.test = "cluster"
#' )
#' }
#'
#'
#'
degGeneExpressionMarkers <- function (object, samples = c(""), contrasts, parameter.to.test = "", 
    test.to.use = "MAST", log.fc.threshold = 0.25, use.spark = FALSE, 
    assay.to.use = "SCT") 
{
  
  ## --------------- ##
  ## Functions       ##
  ## --------------- ##
  
    .getDegTable <- function(n) {
        first.cluster <- unlist(n)[1]
        second.cluster <- unlist(n)[2]
        second.cluster.name <- second.cluster
        if (second.cluster == "all") {
            second.cluster <- NULL
        }
        Idents(object.sub) <- param.to.test
        markers = FindMarkers(object.sub, ident.1 = first.cluster, 
            ident.2 = second.cluster, test.use = test.to.use, 
            logfc.threshold = log.fc.threshold, verbose = FALSE, 
            assay = assay.to.use, slot = "counts")
        colnames(markers) <- chartr(old = " ", new = "_", paste(colnames(markers), 
            first.cluster, "vs", second.cluster.name, sep = "_"))
        return(markers)
    }

    .isDiscreteMetadataColumn <- function(values) {
        values <- values[!is.na(values)]
        if (length(values) == 0) {
            return(FALSE)
        }
        if (is.factor(values) || is.character(values) || is.logical(values)) {
            return(TRUE)
        }
        if (is.numeric(values) || is.integer(values)) {
            unique.values <- unique(values)
            return(length(unique.values) <= 100 && all(abs(unique.values - round(unique.values)) < sqrt(.Machine$double.eps)))
        }
        FALSE
    }
    
    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
    
    # Getting metadata and checking sample names:
    metadata.table <- object@meta.data
    
    if(any(grepl('c\\(|\\[\\]',samples))) {
      samples = eval(parse(text = gsub('\\[\\]', 'c()', samples)))
    }else{
      samples=samples
    }
    samples <- as.character(samples)
    samples <- samples[nzchar(trimws(samples))]

    colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", 
        colnames(object@meta.data))
    sample.column <- if ("orig.ident" %in% colnames(object@meta.data)) {
        "orig.ident"
    } else if ("sample_name" %in% colnames(object@meta.data)) {
        "sample_name"
    } else {
        stop("No sample column found. Expected 'orig.ident', 'orig_ident', or 'sample_name' in object metadata.")
    }
    possible.samples <- sort(unique(as.character(object@meta.data[[sample.column]])))
    possible.samples <- possible.samples[nzchar(trimws(possible.samples))]
    print("Possible sample names:")
    print(possible.samples)
    if (length(samples) == 0) {
        samples = possible.samples
    }
    if ("active.ident" %in% slotNames(object)) {
        sample.name = as.factor(object@meta.data[[sample.column]])
        names(sample.name) = names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        object.sub = subset(object, ident = samples)
    }
    else {
        sample.name = as.factor(object@meta.data[[sample.column]])
        names(sample.name) = names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        object.sub = subset(object, ident = samples)
    }
    print("selected samples:")
    print(object.sub)
    colnames(object.sub@meta.data) = gsub("\\.", "_", colnames(object.sub@meta.data))
    
    #define contrasts
    new.cont <- list()
    for (i in 1:length(contrasts)) {
        new.cont[[i]] <- c(paste(unlist(strsplit(contrasts[i], 
            "-"))))
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
    metadata.columns <- colnames(object.sub@meta.data)
    discrete.metadata.columns <- metadata.columns[vapply(object.sub@meta.data, .isDiscreteMetadataColumn, logical(1))]
    discrete.metadata.column.options <- paste0("  - ", discrete.metadata.columns, collapse = "\n")
    param.to.test <- gsub("\\.", "_", parameter.to.test)
    param.to.test <- as.character(param.to.test)
    param.to.test <- param.to.test[nzchar(trimws(param.to.test))]
    
    if (length(param.to.test) == 0) {
        stop(paste0(
            "parameter.to.test is required.\n",
            "Possible parameter.to.test columns with categorical/discrete values:\n",
            if (length(discrete.metadata.columns) == 0) "  - <none>" else discrete.metadata.column.options
        ), call. = FALSE)
    }
    param.to.test <- param.to.test[[1]]
    if (!(param.to.test %in% metadata.columns)) {
        print(paste(param.to.test, "is not a valid parameter.to.test metadata column."))
        print("Possible parameter.to.test columns with categorical/discrete values:")
        cat(discrete.metadata.column.options, "\n")
        stop("You have entered an invalid metadata column for parameter.to.test.")
    }
    contrast.target <- object.sub@meta.data[[param.to.test]]
    possible.parameter.values <- sort(unique(as.character(contrast.target)))
    possible.parameter.values <- possible.parameter.values[nzchar(trimws(possible.parameter.values))]
    print(paste0("Possible values for ", param.to.test, ":"))
    cat(paste0("  - ", possible.parameter.values, collapse = "\n"), "\n")
    contrast.type <- param.to.test
    contrast.counts = as.data.frame(table(contrast.target))
    valid.contrasts = subset(contrast.counts, Freq > 2)[[1]]
    
    #catch malformed contrasts
    for (i in contrasts) {
        if (!(i[[1]] %in% contrast.target)) {
            print(paste(i[[1]], "is not a valid contrast for contrast type:", 
                contrast.type))
            print("Please see below for an example of valid contrasts\n             for your selected contrast type.")
            print(valid.contrasts)
            stop("You have entered an invalid group to contrast against.")
        }
        else if (!(i[[2]] %in% contrast.target) & (i[[2]] != 
            "all")) {
            print(paste(i[[2]], "is not a valid contrast for contrast type:", 
                contrast.type))
            print("Please see below for an example of valid contrasts\n             for your selected contrast type.")
            print(valid.contrasts)
            stop("You have entered an invalid group to contrast against.")
        }
        else if (length(i) > 2) {
            print("Contrasts are as follows..")
            print(i)
            stop("The console says there are too many inputs in your contrasts.\n         A contrast should only contain Group1-Group2,\n         but the console thinks you have inputed Group1-Group2-Group3")
        }
        else if (!(i[[2]] %in% valid.contrasts) & (i[[2]] != 
            "all")) {
            print(paste(i[[2]], "has two few values (less than 3 cells) to contrast against.\n           Please see below for contrasts with enough cells:", 
                valid.contrasts))
            stop("You have entered an invalid group to contrast against.")
        }
        else if (!(i[[1]] %in% valid.contrasts)) {
            print(paste(i[[1]], "has two few values (less than 3 cells) to contrast against.\n           Please see below for contrasts with enough cells:", 
                valid.contrasts))
            stop("You have entered an invalid group to contrast against.")
        }
    }
    
    #print out contrast cell contrast.counts
    for (i in seq_along(contrasts)) {
        first.group <- contrasts[[i]][[1]]
        first.group.count <- subset(contrast.counts, contrast.target == 
            first.group)$Freq
        if (contrasts[[i]][[2]] != "all") {
            second.group <- contrasts[[i]][[2]]
            second.group.count <- subset(contrast.counts, contrast.target == 
                second.group)$Freq
            print(paste("Contrast No.", i, "contrasts cluster", 
                first.group, "with", first.group.count, "cells vs. cluster", 
                second.group, "with", second.group.count, "cells."))
        }
        else {
            second.group.count <- ncol(object.sub) - first.group.count
            print(paste("Contrast No.", i, "contrasts cluster", 
                first.group, "with", first.group.count, "cells vs. all other clusters, totalling", 
                second.group.count, "cells."))
        }
    }
    if (use.spark) {
        deg.tables <- spark.lapply(contrasts, .getDegTable)
    }
    else {
        deg.tables <- lapply(contrasts, .getDegTable)
    }
    for (i in seq_along(deg.tables)) {
        degtab <- deg.tables[[i]]
        pos <- degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 
            0) %>% dim()
        neg <- degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 
            0) %>% dim()
        print(paste0("The number of upregulated genes at p<0.05 in contrast number ", 
            i, " is:"))
        print(pos[1])
        print(paste0("The number of downregulated genes at p<0.05 in contrast number ", 
            i, " is:"))
        print(neg[1])
    }
    
    #Merge the deg tables together
    out.df <- NULL
    for (i in deg.tables) {
        if (is.null(out.df)) {
            out.df <- deg.tables[1]
            out.df <- as.data.frame(out.df)
        }
        else {
            out.df <- merge(out.df, i, by = "row.names", all = TRUE)
            rownames(out.df) <- out.df$Row.names
            out.df$Row.names <- NULL
        }
    }
    out.df$Gene <- rownames(out.df)
    out.df$Row.names <- NULL
    out.df <- out.df %>% dplyr::select(Gene, everything())

### Purpose of this edit is to bring output table in-line with 
### column name expectations from GSEA downstream. Column names 
### should be more similar to those from Bulk DEG Analysis output
### when we are done.

# Original column names
original_colnames <- colnames(out.df)

# Function to rename columns based on pattern matching
rename_columns <- function(colnames, pattern, suffix) {
  # Identify columns that match the pattern
  matched <- grepl(pattern, colnames)
  
  # Apply the transformation only to the matching columns
  colnames[matched] <- gsub(
    pattern = paste0("^", pattern),             # Match the specific pattern at the start
    replacement = "C_",                         # Replace pattern with "C_"
    x = colnames[matched]
  )
  
  # Append the suffix to the columns that were matched
  colnames[matched] <- gsub(
    pattern = "C_(.*)_vs_(.*)",                 # Capture the parts after "C_"
    replacement = paste0("C_\\1_vs_\\2", suffix), # Reconstruct the name and append the suffix
    x = colnames[matched]
  )
  
  return(colnames)
}

# Apply the renaming function for various patterns with their corresponding suffixes
new_colnames <- original_colnames
new_colnames <- rename_columns(new_colnames, "p_val_adj_", "_adjpval")
new_colnames <- rename_columns(new_colnames, "avg_log2FC_", "_logFC")
new_colnames <- rename_columns(new_colnames, "pct.1_", "_pct1")
new_colnames <- rename_columns(new_colnames, "pct.2_", "_pct2")
new_colnames <- rename_columns(new_colnames, "p_val_", "_pval")

# Update the column names in the dataframe
colnames(out.df) <- new_colnames

### End of "shift" edits.



    result.list <- list("data" = list("DEG_Table"=out.df))
    return(result.list)
}
