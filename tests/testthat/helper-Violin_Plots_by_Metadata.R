selectViolin <- function(dataset) {
  
  if (dataset == "TEC"){
    
    object = selectCRObject("TEC")
    
    object = selectCRObject("TEC")
    group = "orig_ident"
    assay = 'SCT'
    layer = 'scale.data'
    jitter.points = T
    jitter.dot.size = 4
    filter.outliers = F
    outlier.low = 0.05
    outlier.high = 0.95
    facet.by = ""
    set.seed(81)
    genes = sample(rownames(Seurat::GetAssayData(object, assay = assay, layer = layer)), 5, 
             replace = FALSE)
    } else if (dataset == "Chariou"){

    object = selectCRObject("Chariou")
    group = "orig_ident"
    assay = 'SCT'
    layer = 'scale.data'
    jitter.points = T
    jitter.dot.size = 4
    filter.outliers = F
    outlier.low = 0.05
    outlier.high = 0.95
    facet.by = ""
    set.seed(82)
    genes = sample(rownames(Seurat::GetAssayData(object, assay = assay, layer = layer)), 5,
             replace = FALSE)

  # } else if (dataset == "Chariou.allgroups"){
  # 
  #   object = selectCRObject("Chariou")
  #   group = "orig_ident"
  #   assay = 'SCT'
  #   slot = 'scale.data'
  #   jitter.points = T
  #   jitter.dot.size = 4
  #   filter.outliers = F
  #   outlier.low = 0.05
  #   outlier.high = 0.95
  #   facet.by = ""
  #   set.seed(821)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5,
  #                              replace = FALSE)
  # 
  # 
  # }else if (dataset == "Chariou.subgroup"){
  # 
  #   object = selectCRObject("Chariou")
  #   group = "orig_ident"
  #   assay = 'SCT'
  #   slot = 'scale.data'
  #   jitter.points = T
  #   jitter.dot.size = 4
  #   filter.outliers = F
  #   outlier.low = 0.05
  #   outlier.high = 0.95
  #   facet.by = ""
  #   set.seed(822)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5,
  #                              replace = FALSE)

  } else if (dataset == "pbmc.single"){

    object = selectCRObject("pbmc-single")
    group = "orig_ident"
    assay = 'SCT'
    layer = 'scale.data'
    jitter.points = T
    jitter.dot.size = 4
    filter.outliers = F
    outlier.low = 0.05
    outlier.high = 0.95
    facet.by = ""
    set.seed(83)
    genes = sample(rownames(Seurat::GetAssayData(object, assay = assay, layer = layer)), 5,
             replace = FALSE)


  } else if (dataset == "nsclc.multi"){

    object = selectCRObject("nsclc-multi")
    group = "orig_ident"
    assay = 'SCT'
    layer = 'scale.data'
    jitter.points = T
    jitter.dot.size = 4
    filter.outliers = F
    outlier.low = 0.05
    outlier.high = 0.95
    facet.by = ""
    set.seed(84)
    genes = sample(rownames(Seurat::GetAssayData(object, assay = assay, layer = layer)), 5,
             replace = FALSE)

  } else if (dataset == "brca"){

    object = selectCRObject("BRCA")
    group = "orig_ident"
    assay = 'SCT'
    layer = 'scale.data'
    jitter.points = T
    jitter.dot.size = 4
    filter.outliers = F
    outlier.low = 0.05
    outlier.high = 0.95
    facet.by = ""
    set.seed(85)
    genes = sample(rownames(Seurat::GetAssayData(object, assay = assay, layer = layer)), 5,
             replace = FALSE)}
  
  return(list("object" = object, 
              "group" = group,
              "assay" = assay,
              "layer" = layer,
              "jitter.points" = jitter.points,
              "jitter.dot.size" = jitter.dot.size,
              "genes" = genes))
  }

.drawViolin <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
