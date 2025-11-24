#' @title Visualize gene expression for provided Genes
#' across your cells as a heatmap
#' @description You should see one plot (TSNE or UMAP, your choice)
#' per gene name provided. The intensity of the red color will be relative
#' to the expression of the gene in that cell
#' @details This function must be run downstream of the Sample Names function,
#' as well as be provided a combined Seurat Object
#' such as the one produced by the SingleR Cell Annotation function
#'
#' @param object Object of class Seurat
#' @param samples.to.include Samples to be included in the analysis
#' @param gene Genes which you would like to visualize
#' @param reduction.type  Select the kind of clustering visualization
#' you would like to use to visualize the cell type results
#' ("umap", "tsne", "pca"). Default is "umap"
#' @param number.of.rows The number of rows you want to arrange your plots into
#' @param return.seurat.object Set to FALSE if you want only a geneset
#' (and not a Seurat object) to be returned
#' @param color The color you want to use in your heatmap (default "red")
#' @param point.size The size of the points representing
#' each cell in your visualization. Default is 1
#' @param point.shape The code for your point shape (R "pch" argument).
#' Default is 16
#' @param point.transparency Set the transparency. Default is 0.5
#' @param use.cite.seq.data TRUE if you would like to plot Antibody clusters
#' from CITEseq instead of scRNA.
#' @param assay Select Assay to Plot (default is SCT)
#'
#' @import Seurat
#' @import gridExtra
#' @import ggplot2
#' @import tidyverse
#'
#'
#' @export
#'
#' @return a Seurat object with additional metadata or gene table and plot



colorByGene <- function(object,
                        samples.to.include,
                        gene,
                        reduction.type = "umap",
                        number.of.rows = 0,
                        return.seurat.object = FALSE,
                        color = "red",
                        point.size = 1,
                        point.shape = 16,
                        point.transparency = 0.5,
                        use.cite.seq.data = FALSE,
                        assay = "SCT") {
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  
  print(object)
  # checking for samples
  if(any(grepl('c\\(|\\[\\]',samples))) {
    samples = eval(parse(text = gsub('\\[\\]', 'c()', samples)))
  }else{
    samples=samples
  }
  # if none specified, using ALL
  if (length(samples) == 0) {
    samples = unique(object@meta.data$orig.ident)
  }
  
  # Fix for underscore
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
  
  #Check input for missing genes
  no.gene = gene[!gene %in% rownames(object.sub[[assay]]@scale.data)]
  
  if (!is.null(no.gene)) {
    print("Gene(s) missing from dataset:")
    print(no.gene)
  }
  
  gene = gene[gene %in% rownames(object.sub[[assay]]@scale.data)]
  
  if (length(gene) > 0) {
    .plotGene <- function(gene) {
      gene.mat = object.sub[[assay]]@scale.data[gene,]
      gene.quant = quantile(gene.mat[gene.mat > 1], probs = c(.1, .5, .90))
      gene.mat[gene.mat > gene.quant[3]] = gene.quant[3]
      gene.mat[gene.mat < gene.quant[1]] = 0
      
      if (!(use.cite.seq.data)) {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object.sub,
                        reduction = "tsne",
                        group.by = "ident")
          colnames(p1$data) <- gsub("tsne_","tSNE_",colnames(p1$data))
          clus.mat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            gene = gene.mat
          )
        }
        else if (reduction.type == "umap") {
          p1 <- DimPlot(object.sub,
                        reduction = "umap",
                        group.by = "ident")
          colnames(p1$data) <- gsub("umap_","UMAP_",colnames(p1$data))
          clus.mat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            gene = gene.mat
          )
        } else {
          p1 <- DimPlot(object.sub,
                        reduction = "pca",
                        group.by = "ident")
          colnames(p1$data) <- gsub("pc_","PC_",colnames(p1$data))
          clus.mat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            gene = gene.mat
          )
        } #if CITEseq is chosen then:
      } else {
        if (reduction.type == "tsne") {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_tsne",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            gene = gene.mat
          )
        }
        else if (reduction.type == "umap") {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_umap",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            gene = gene.mat
          )
        } else {
          p1 <-
            DimPlot(object.sub,
                    reduction = "protein_pca",
                    group.by = "ident")
          clus.mat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            gene = gene.mat
          )
        }
      }
      
      # Set reduction type x and y coordinates
      reduction.type.x <- paste0(reduction.type, "-1")
      reduction.type.y <- paste0(reduction.type, "-2")
      
      clus.mat %>% dplyr::arrange(gene) -> clus.mat
      g <- ggplot(clus.mat, aes(x = umap1, y = umap2)) +
        theme_bw() +
        theme(legend.title = element_blank()) +
        ggtitle(gene) +
        geom_point(
          aes(colour = gene),
          alpha = point.transparency,
          shape = point.shape,
          size = point.size
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = rel(0.5))
        ) +
        scale_color_gradient(
          limits = c(0, gene.quant[3]),
          low = "lightgrey",
          high = color
        ) +
        xlab(reduction.type.x) + ylab(reduction.type.y)
      return(g)
    }
    
    if (number.of.rows == 0) {
      n = ceiling(length(gene) ^ 0.5)
    } else {
      n = number.of.rows
    }
    
    
    grob <- lapply(seq_along(gene), function(x)
      .plotGene(gene[x]))
    ##    gridExtra::grid.arrange(grobs=grob,nrow=n,newpage=F)
    ###    plots <- gridExtra::grid.arrange(grobs=grob,nrow=n,newpage=F)
    
    if (return.seurat.object) {
      result.list <- list("object" = object, "plot" = grob)
      return(result.list)
    } else {
      gene = as.data.frame(gene)
      result.list <- list("object" = gene, "plot" = grob)
      return(result.list)
    }
    
  } else {
    print("No genes found in dataset")
    return(NULL)
  }
  
}
