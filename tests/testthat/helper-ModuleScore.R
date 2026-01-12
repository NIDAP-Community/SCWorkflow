getModuleScoreParam <- function(data){

  if(data == "tec"){
    
      object = selectCRObject("TEC")
      marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      ms_threshold = paste(colnames(marker.table)[1:3], rep(0, ncol(marker.table)))
      use_columns = colnames(marker.table)[1:3]
      general.class = colnames(marker.table)[1:3]
      # lvl.vec = c('Pan_Tcells-CD4_T-Tregs','Pan_Tcells-CD4_T-New','Pan_Tcells-CD8_T')

  } else if (data == "chariou") {
    
      object = selectCRObject("Chariou")
      marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      ms_threshold = paste(colnames(marker.table)[1:3], rep(0, ncol(marker.table)))
      use_columns = colnames(marker.table)[1:3]
      general.class = colnames(marker.table)[1:3]
      # lvl.vec = c('Pan_Tcells-CD4_T-Tregs','Pan_Tcells-CD4_T-New','Pan_Tcells-CD8_T')
      
  } else if (data == "pbmc.single") {
    
    object = selectCRObject("pbmc-single")
    set.seed(114)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    ms_threshold = paste(colnames(marker.table), rep(0, ncol(marker.table)))
    use_columns = colnames(marker.table)
    general.class = colnames(marker.table)
    # lvl.vec = c('Pan_Tcells-CD4_T-Tregs','Pan_Tcells-CD4_T-New','Pan_Tcells-CD8_T')
    
  } else if (data == "nsclc.multi") {

    object = selectCRObject("nsclc-multi")
    set.seed(214)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5,
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    ms_threshold = paste(colnames(marker.table), rep(0, ncol(marker.table)))
    use_columns = colnames(marker.table)
    general.class = colnames(marker.table)
    # lvl.vec = c('Pan_Tcells-CD4_T-Tregs','Pan_Tcells-CD4_T-New','Pan_Tcells-CD8_T')

  } else if (data == "brca") {

    object = selectCRObject("BRCA")
    set.seed(314)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    ms_threshold = paste(colnames(marker.table), rep(0, ncol(marker.table)))
    use_columns = colnames(marker.table)
    general.class = colnames(marker.table)
    # lvl.vec = c('Pan_Tcells-CD4_T-Tregs','Pan_Tcells-CD4_T-New','Pan_Tcells-CD8_T')
    
  }
  
  return(list("object" = object, 
              "ms_threshold"= ms_threshold, 
               # "lvl.vec" = lvl.vec, 
              "marker.table" = marker.table, 
              "general.class" = general.class,
              "use_columns"=use_columns
              ))  
}

.drawMSfig <- function(x, width = 10, height = 10, index = 1){
  target <- x
  if (is.list(x) && all(c("object", "figures") %in% names(x))) {
    if (length(x$figures) < index) {
      stop("Requested index exceeds available figures")
    }
    target <- x$figures[[index]]
  }
  path <- tempfile(fileext = ".png")
  ggsave(path, target, width = 10, height = 10)
  print(path)
}
