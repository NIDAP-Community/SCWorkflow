#' @title Heatmap of transcript and/or protein expression values in single cells
#' @description This method provides a heatmap of single cell data from a Seurat
#'  object given a set of genes and optionally orders by various metadata and/or
#'  gene or protein expression levels. Method is based on ComplexHeatmap::pheatmap
#'
#' @param object Seurat-class object
#' @param sample.names Sample names
#' @param metadata Metadata column to plot
#' @param transcripts Transcripts to plot
#' @param proteins Proteins to plot (default is NULL)
#' @param heatmap.color Color for heatmap. Choices are "Cyan to Mustard",
#'   "Blue to Red", "Red to Vanilla", "Violet to Pink", "Bu Yl Rd", 
#'   "Bu Wt Rd" (default is "Bu Yl Rd")
#' @param plot.title Title of plot (default is "Heatmap")
#' @param add.gene.or.protein Add Gene or protein annotations (default is FALSE)
#' @param protein.annotations Protein annotations to add (defulat is NULL)
#' @param rna.annotations Gene annotations to add (default is NULL)
#' @param arrange.by.metadata Arrange by metadata (default is TRUE)
#' @param add.row.names Add row names (default is TRUE)
#' @param add.column.names Add column names (default is FALSE)
#' @param row.font Font size for rows (default is 5)
#' @param col.font Font size for columns (default is 5)
#' @param legend.font Font size for legend (default is 5)
#' @param row.height Height of row. If NA, adjust to plot size (default is 15)
#' @param set.seed Seed for colors (default is 6)
#' @param scale.data Perform z-scaling on rows (default is TRUE)
#' @param trim.outliers Remove outlier data (default is TRUE)
#' @param trim.outliers.percentage Set outlier percentage (default is 0.01)
#' @param order.heatmap.rows Order heatmap rows (default is FALSE)
#' @param row.order Gene vector to set row order. If NULL, use cluster order
#'  (default is NULL)
#'
#' @import Seurat
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom dendsort dendsort
#' @importFrom dplyr filter arrange across all_of mutate_if select
#' @importFrom dendextend rotate
#' @importFrom stats hclust kmeans as.hclust
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace_all str_wrap
#' @importFrom colorspace RGB diverge_hcl heat_hcl hex
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#' 
#' @return This function returns a heatmap plot and the data underlying the 
#'  heatmap.
#'
heatmapSC <- function(object,
                      sample.names,
                      metadata,
                      transcripts,
                      use_assay = 'SCT',
                      proteins = NULL,
                      heatmap.color = "Bu Yl Rd",
                      plot.title = "Heatmap",
                      add.gene.or.protein = FALSE,
                      protein.annotations = NULL,
                      rna.annotations = NULL,
                      arrange.by.metadata = TRUE,
                      add.row.names = TRUE,
                      add.column.names = FALSE,
                      row.font = 5,
                      col.font = 5,
                      legend.font = 5,
                      row.height = 15,
                      set.seed = 6,
                      scale.data = TRUE,
                      trim.outliers = TRUE,
                      trim.outliers.percentage = 0.01,
                      order.heatmap.rows = FALSE,
                      row.order = c()) {
  
  #### Functions ####
  n <- 2e3
  set.seed(set.seed)
  color.space <- colorspace::RGB(runif(n), runif(n), runif(n))
  color.space <- as(color.space, "LAB")
  
  #function to create large palette of colors for annotation tracks 
  .distinctColorPalette <- function(k = 1, seed) {
    current.color.space <- color.space@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(current.color.space, k, iter.max = 20)
    colors <- unname(colorspace::hex(colorspace::LAB(km$centers)))
    return(colors)
  }
  
  ## Function to create cyan to mustard palette
  .pal <- function (n,
                    h = c(237, 43),
                    c = 100,
                    l = c(70, 90),
                    power = 1,
                    fixup = TRUE,
                    gamma = NULL,
                    alpha = 1,
                    ...) {
    if (n < 1L)
      return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- colorspace::hex(colorspace::polarLUV(
      L = l[2L] - diff(l) * abs(rval) ^ power[2L],
      C = c * abs(rval) ^ power[1L],
      H = ifelse(rval > 0, h[1L], h[2L])
    ),
    fixup = fixup,
    ...)
    if (!missing(alpha)) {
      alpha <- pmax(pmin(alpha, 1), 0)
      alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
                      width = 2L,
                      upper.case = TRUE)
      rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  }
  
  #Color selections for heatmap:
  np0 = .pal(100) #Cyan to Mustard
  np1 = colorspace::diverge_hcl(100,
                    c = 100,
                    l = c(30, 80),
                    power = 1)  #Blue to Red
  np2 = colorspace::heat_hcl(
    100,
    c = c(80, 30),
    l = c(30, 90),
    power = c(1 / 5, 2)
  ) #Red to Vanilla
  np3 = rev(colorspace::heat_hcl(
    100,
    h = c(0, -100),
    c = c(40, 80),
    l = c(75, 40),
    power = 1
  )) #Violet to Pink
  np4 = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)) #Bu Yl Rd
  np5 = grDevices::colorRampPalette(c("steelblue", "white", "red"))(100)
  #Steelblue to White to Red
  
  np = list(np0, np1, np2, np3, np4, np5)
  names(np) = c("Cyan to Mustard",
                "Blue to Red",
                "Red to Vanilla",
                "Violet to Pink",
                "Bu Yl Rd",
                "Bu Wt Rd")
  
  # Function to set various heatmap parameters and run pheatmap
  .doHeatmap <- function(dat, clus.col, clus.row, rn, cn, col) {
    col.pal <- np[[col]]
    minx = min(dat)
    maxx = max(dat)
    tree.height <- 25
    breaks = seq(minx, maxx, length = 100)
    breaks = sapply(breaks, signif, 4)
    leg.breaks = seq(minx, maxx, length = 5)
    leg.breaks = sapply(leg.breaks, signif, 4)
    
    #Run cluster method:
    hc = hclust(dist(t(dat)), method = "complete")
    hcrow = hclust(dist(dat), method = "complete")
    
    if (clus.col == TRUE) {
      .sortHClust <-
        function(...)
          as.hclust(rev(dendsort::dendsort(as.dendrogram(...))))
    } else {
      .sortHClust <- function(...)
        as.hclust(dendsort::dendsort(as.dendrogram(...)))
    }
    
    if (clus.row == TRUE) {
      row.clus <- .sortHClust(hcrow)
    } else {
      row.clus = FALSE
    }
    
    plot.title <- stringr::str_replace_all(plot.title, "_", " ")
    plot.title <- stringr::str_wrap(plot.title, 50)
    labels.col <- colnames(dat)
    
    #Set up heatmap parameters for running pheatmap
    hm.parameters <- list(
      dat,
      color = col.pal,
      legend_breaks = leg.breaks,
      scale = "none",
      treeheight_col = tree.height,
      treeheight_row = tree.height,
      kmeans_k = NA,
      breaks = breaks,
      fontsize = legend.font,
      fontsize_row = row.font,
      fontsize_col = col.font,
      cellheight = row.height,
      show_rownames = rn,
      show_colnames = cn,
      main = plot.title,
      clustering_method = "complete",
      cluster_rows = row.clus,
      cluster_cols = clus.col,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      annotation_col = annotation.col,
      annotation_colors = annot.col,
      labels_col = labels.col
    )
    
    mat = t(dat)
    
    callback = function(hc, mat) {
      dend = rev(dendsort::dendsort(as.dendrogram(hc)))
      #dend %>% dendextend::rotate(c(1:length(dend))) -> dend
      as.hclust(dend)
    }
    do.call(ComplexHeatmap::pheatmap, c(hm.parameters, list(clustering_callback = callback)))
  }
  
  ##### Main Code Block ####
  
  # Select samples to include
  samples.to.include = sample.names
  samples.to.include <- samples.to.include[samples.to.include != ""]
  samples.to.include <- gsub("-", "_", samples.to.include)
  
  #Error messaging for metadata 
  
  if(is.null(metadata)){
    stop("Error: You should choose at least one annotation track under metadata_to_plot")
  }
  
  if(sum(grepl("Barcode",metadata,ignore.case=TRUE)) > 0){
    sprintf("Annotation Track cannot include Barcode")
    metadata <- metadata[!grepl('Barcode', metadata, ignore.case=TRUE)]
  }
  
  #Clean up transcript names and print missing genes:
  transcripts = gsub(" ", "", transcripts)
  
  l1 <- length(transcripts)
  p1 <- length(proteins)
  
  if(l1 + p1 == 0){
    stop(sprintf("At least 1 transcript and/or protein is needed for plotting"))
  }
  
  dups <- transcripts[duplicated(transcripts)]
  transcripts <- transcripts[!duplicated(transcripts)]
  
  
  l2 <- length(transcripts)
  sprintf("There are %s total unique genes/proteins in the dataset", l2)
  if (l1 > l2) {
    warning(sprintf("\n\nThe following duplicate genes were removed: %s",
                    dups))
  }
  
  assay.for.transcripts <- if (use_assay == "Harmony") "Harmony" else "SCT"
  assay.features <- rownames(GetAssayData(object, assay = assay.for.transcripts, layer = "scale.data"))
  missing.genes <- character(0)
  if (length(transcripts) > 0 && transcripts[1] != "") {
    missing.genes = setdiff(transcripts, assay.features)
    if (length(missing.genes) > 0) {
      print(paste("missing genes:", missing.genes))
    }
  }
  
  #missing.genes <- transcripts[!transcripts %in% rownames(object)]
  l3 <- length(missing.genes)
  missing.genes <- paste(shQuote(missing.genes), collapse = ", ")
  if (l3 == l2) {
    stop("No genes listed are found in dataset.")
  }
  if (l3 > 0) {
    warning(
      sprintf(
        "There are %s gene(s) absent from dataset:
            %s. Possible reasons are that gene is not official
            gene symbol or gene is not highly expressed and
            has been filtered.\n ",
        l3,
        missing.genes
      )
    )
  }
  
  transcripts <- transcripts[transcripts %in% assay.features]
  
  .assayNames <- function(so) {
    a <- Assays(so)
    n <- names(a)
    if (!is.null(n) && length(n) > 0) {
      return(n)
    }
    tryCatch(as.character(a), error = function(e) character(0))
  }

  #Clean up protein names and print missing proteins:
  if ("Protein" %in% .assayNames(object)) {
    proteins = gsub(" ", "", proteins)
    if (proteins[1] != "") {
      protmiss = setdiff(proteins, rownames(GetAssayData(object, assay="Protein", layer="scale.data")))
      if (length(protmiss) > 0) {
        sprintf("missing proteins: %s", protmiss)
      }
    }
    proteins = proteins[proteins %in% rownames(GetAssayData(object, assay="Protein", layer="scale.data"))]
  }
  
  #Error messaging for protein annotation tracks:
  
  if(add.gene.or.protein == FALSE & (!is.null(protein.annotations) | !is.null(rna.annotations))) {
    stop("Error: You should choose to add gene or protein annotation tracks if you add protein or rna annotations")
  }
  
  
  #collect transcript expression data from SCT / Harmony slot
  df.mat1 = NULL
  if (length(transcripts) > 0) {
    if (length(transcripts) == 1) {
      if(use_assay == 'SCT'){
        df.mat1 <-
          vector(mode = "numeric",
                 length = length(GetAssayData(object, assay="SCT", layer="scale.data")[transcripts,]))
        df.mat1 <- GetAssayData(object, assay="SCT", layer="scale.data")[transcripts,]
      } else if (use_assay == 'Harmony'){
        df.mat1 <-
          vector(mode = "numeric",
                 length = length(GetAssayData(object, assay="Harmony", layer="scale.data")[transcripts,]))
        df.mat1 <- GetAssayData(object, assay="Harmony", layer="scale.data")[transcripts,]
        }
    } else {
      if(use_assay == 'SCT'){
        df.mat1 <- as.matrix(GetAssayData(object, assay="SCT", layer="scale.data")[transcripts,])
      } else if (use_assay == 'Harmony'){
        df.mat1 <- as.matrix(GetAssayData(object, assay="Harmony", layer="scale.data")[transcripts,])
        }
    }
  }
  
  #collect protein expression values
  df.mat2 = NULL
  if ("Protein" %in% .assayNames(object)) {
    if (length(proteins) > 0) {
      if (length(proteins) == 1) {
        df.mat2 <-
          vector(mode = "numeric",
                 length = length(GetAssayData(object, assay="Protein", layer="scale.data")[proteins,]))
        df.mat2 <- GetAssayData(object, assay="Protein", layer="scale.data")[proteins,]
        protname <- paste0(proteins, "_Prot")
      } else {
        df.mat2 <- as.matrix(GetAssayData(object, assay="Protein", layer="scale.data")[proteins,])
        protname <- paste0(proteins, "_Prot")
        rownames(df.mat2) <- protname
      }
    }
  }
  
  #Put together transcript and protein data
  df.mat <- rbind(df.mat1, df.mat2)
  if (!is.null(df.mat1)) {
    rownames(df.mat)[rownames(df.mat) == "df.mat1"] <- transcripts
  }
  
  if (!is.null(df.mat2)) {
    rownames(df.mat)[rownames(df.mat) == "df.mat2"] <- protname
  }
  
  df.mat <- df.mat[sort(rownames(df.mat)),]
  
  #Set row order to specific genes/proteins
  if (order.heatmap.rows == TRUE) {
    row.order = row.order[row.order %in% rownames(df.mat)]
    row.order = c(row.order, setdiff(rownames(df.mat), row.order))
    df.mat <- df.mat[row.order,]
    cluster.rows <- FALSE
  } else {
    cluster.rows <- TRUE
  }
  
  metadata <- sub("orig_ident", "orig.ident", metadata)
  
  #Set up annotation track:
  metadata.table <- object@meta.data %>%
    dplyr::filter(orig.ident %in% samples.to.include)
  
  annot <- metadata.table %>% select(all_of(metadata))
  a = dim(annot)[2] - 1
  
  #Addition of gene and/or protein expression data to annotation tracks:
  if (add.gene.or.protein == TRUE) {
    if (length(protein.annotations) > 0) {
      annot1 <-
        as.matrix(GetAssayData(object, assay="Protein", layer="scale.data")[protein.annotations,])
      if (length(protein.annotations) == 1) {
        annot1 <- annot1[match(rownames(annot), rownames(annot1))]
        protname <- paste0(protein.annotations, "_Prot")
      } else {
        annot1 <- annot1[, match(rownames(annot), colnames(annot1))]
        annot1 <- t(annot1)
        colnames(annot1) = paste0(colnames(annot1), "_Prot")
      }
    }
    
    if (length(rna.annotations) > 0) {
      annot2 <- as.matrix(GetAssayData(object, assay="SCT", layer="scale.data")[rna.annotations,])
      if (length(rna.annotations) == 1) {
        annot2 <- annot2[match(rownames(annot), rownames(annot2))]
      } else {
        annot2 <- annot2[, match(rownames(annot), colnames(annot2))]
        annot2 <- t(annot2)
      }
    }
  }
  
  if (exists("annot1")) {
    annot <- cbind(annot, annot1)
    colnames(annot)[colnames(annot) == "annot1"] <- protname
  }
  if (exists("annot2")) {
    annot <- cbind(annot, annot2)
    colnames(annot)[colnames(annot) == "annot2"] <- rna.annotations
  }
  
  #Arrange columns by metadata tracks:
  if (arrange.by.metadata == TRUE) {
    annot <- annot %>% arrange(across(all_of(colnames(annot)))) 
    df.mat <- df.mat[, match(rownames(annot), colnames(df.mat))]
    df.mat <-
      df.mat[, apply(df.mat, 2, function(x)
        ! any(is.na(x)))]
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }
  
  
  #Set up annotation track and colors:
  annotation.col = as.data.frame(unclass(annot[, !(names(annot) %in%
                                                     "Barcode")]))
  annotation.col <- annotation.col %>%
    mutate_if(is.logical, as.factor)
  rownames(annotation.col) <- rownames(annot)
  if (dim(annot)[2] == 1) {
    annottitle = colnames(annot)[1]
    colnames(annotation.col) = annottitle
  }
  annot.col = list()
  groups = colnames(annotation.col)
  colnames(annot) <- groups
  
  q = 0
  for (i in 1:dim(annotation.col)[2]) {
    if (class(annot[[groups[i]]]) != "numeric") {
      annotation.col[, i] = factor(annotation.col[, i])
      q = q + length(levels(annotation.col[, i]))
    }
  }
  
  colors = .distinctColorPalette(q, 5)
  
  b = 1
  i = 1
  nam = NULL
  col <- NULL
  annot.col <- NULL
  for (i in 1:length(groups)) {
    nam <- groups[i]
    if (class(annotation.col[, i]) != "numeric") {
      grp <- as.factor(annotation.col[, i])
      c <- b + length(levels(grp)) - 1
      col = colors[b:c]
      names(col) <- levels(grp)
      assign(nam, col)
      annot.col = append(annot.col, mget(nam))
      b = c + 1
      i = i + 1
    }
    else{
      grp <- annotation.col[, i]
      np5 = colorRampPalette(c("steelblue", "white", "red"))(length(grp))
      col = np5
      assign(nam, col)
      annot.col = append(annot.col, mget(nam))
    }
  }
  
  sprintf("The total number of genes in heatmap: %s", nrow(df.mat))
  
  # Apply scaling to data (optional)
  if (scale.data == TRUE) {
    tmean.scale = t(scale(t(df.mat)))
    tmean.scale = tmean.scale[is.finite(rowSums(tmean.scale)),]
  } else {
    tmean.scale = df.mat
  }
  
  # Apply outlier trimming to data (optional)
  if (trim.outliers == TRUE) {
    quantperc <- trim.outliers.percentage
    upperquant <- 1 - quantperc
    for (i in 1:nrow(tmean.scale)) {
      data <- tmean.scale[i,]
      dat.quant <- quantile(data, probs = c(quantperc, upperquant))
      data[data > dat.quant[2]] <- dat.quant[2]
      data[data < dat.quant[1]] <- dat.quant[1]
      tmean.scale[i,] <- data
    }
  }
  
  #Draw heatmap
  p = .doHeatmap(
    dat = tmean.scale,
    clus.col = cluster.cols,
    clus.row = cluster.rows,
    rn = add.row.names,
    cn = add.column.names,
    col = heatmap.color
  )
  
  #Add legends and title formatting:
  p@matrix_color_mapping@name <- " "
  p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
  p@column_title_param$gp$fontsize <- 10
  
  #Return expression matrix used in heatmap
  heatmap.df <- as.data.frame(tmean.scale) %>%
    rownames_to_column("gene")
    
  heatmap.res <- list("data" = heatmap.df, 
                      "plots" = p )
  return(heatmap.res)
}
