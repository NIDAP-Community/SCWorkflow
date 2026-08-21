getParamHM <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    sample.names <-
      c("1_Embryo_13_5", "2_Embryo_15", "3_Newborn", "4_Adult")
    metadata <- c("orig.ident", "Phase")
    set.seed(15)
    transcripts <- sample(rownames(object), 10, replace = FALSE)
    proteins <- NULL
    add.gene.or.protein <- TRUE
    rna.annotations <- transcripts[c(1, 2)]
    protein.annotations <- NULL
    plot.title <- "Heatmap_TEC_test"
    
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    sample.names <- c("PBS", "CD8dep", "ENT", "NHSIL12", "Combo")
    metadata <- c("orig.ident")
    set.seed(15)
    add.gene.or.protein <- TRUE
    transcripts <- sample(rownames(object), 10, replace = FALSE)
    proteins <- NULL
    rna.annotations <- transcripts[c(1, 2)]
    protein.annotations <- NULL
    plot.title <- "Heatmap_Chariou_test"
    
  } else if (data == "Chariou2") {
    object <- selectCRObject("Chariou")
    sample.names <- c("PBS", "CD8dep", "ENT", "NHSIL12", "Combo")
    metadata <- c("orig.ident","Phase")
    set.seed(15)
    add.gene.or.protein <- FALSE
    transcripts <- sample(rownames(object), 10, replace = FALSE)
    proteins <- NULL
    rna.annotations <- NULL
    protein.annotations <- NULL
    plot.title <- "Heatmap_Chariou_test"
  
  } else if (data == "Chariou3") {
    object <- selectCRObject("Chariou")
    sample.names <- c("PBS", "CD8dep", "ENT", "NHSIL12", "Combo")
    metadata <- c("orig.ident","Phase")
    set.seed(15)
    add.gene.or.protein <- FALSE
    transcripts <- NULL
    proteins <- NULL
    rna.annotations <- NULL
    protein.annotations <- NULL
    plot.title <- "Heatmap_Chariou_test"
    
  } else if (data == "pbmc-single") {
    object <- selectSRObject("pbmc-single")
    sample.names <- c("PBMC_Single")
    metadata <- c("HPCA_main", "BP_encode_main")
    set.seed(15)
    transcripts <- sample(rownames(object), 10, replace = FALSE)
    proteins <- NULL
    add.gene.or.protein <- TRUE
    rna.annotations <- transcripts[c(1, 2)]
    protein.annotations <- NULL
    
    plot.title <- "Heatmap_Single_NSCLC"
  } else if (data == "nsclc-multi") {
    object <- selectCRObject("nsclc-multi")
    sample.names <-
      c("Donor_1",
        "Donor_2",
        "Donor_3",
        "Donor_4",
        "Donor_5",
        "Donor_6",
        "Donor_7")
    metadata <- c("orig.ident", "Phase")
    set.seed(15)
    transcripts <- sample(rownames(object), 10, replace = FALSE)
    proteins <- c("CD3", "CD8")
    add.gene.or.protein <- TRUE
    rna.annotations <- transcripts[c(1, 2)]
    protein.annotations <-  c("CD3", "CD8")
    plot.title <- "Heatmap_Multi_NSCLC"
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
    sample.names <-
      c("CID4471",
        "CID4290A",
        "CID44971",
        "CID4040",
        "CID4513",
        "CID4535")
    metadata <- c("orig.ident", "SCT_snn_res.0.2")
    set.seed(15)
    transcripts <-
     sample(rownames(object@assays$SCT@scale.data), 10,
            replace = FALSE)
    proteins <- NULL
    add.gene.or.protein <- TRUE
    rna.annotations <- transcripts[c(1, 2)]
    protein.annotations <- NULL
    plot.title <- "Heatmap_BRCA"
  }
  
  
  
  return(
    list(
      "object" = object,
      "sample.names" = sample.names,
      "metadata" = metadata,
      "transcripts" = transcripts,
      "proteins" = proteins,
      "plot.title" = plot.title,
      "add.gene.or.protein" = add.gene.or.protein,
      "rna.annotations" = rna.annotations,
      "protein.annotations" = protein.annotations
    )
  )
}

.drawHeatPng <- function(x,width=1200, height = 500){
  path <- tempfile(fileext = ".png")
  png(path,
      width=width,
      height=height,
      )
  on.exit(dev.off())
  print(x)
  path
}

