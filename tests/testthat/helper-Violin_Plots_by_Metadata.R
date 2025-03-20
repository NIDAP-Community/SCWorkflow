selectViolin <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    object = selectCRObject("TEC")
    group = "orig.ident"
    assay = 'SCT'
    slot = 'scale.data'
    jitter_points = T
    jitter_dot_size = 4
    facet_by = ""
    set.seed(81)
    genes = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    }
    
  # } else if (dataset == "Chariou"){
  #   
  #   print("selected Chariou dataset") 
  #   object = selectCRObject("Chariou")
  #   group.by = "orig_ident"
  #   group.subset =  unique(object$orig.ident)[1:3]
  #   set.seed(82)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  #   
  # } else if (dataset == "Chariou.allgroups"){
  #   
  #   print("selected Chariou dataset") 
  #   object = selectCRObject("Chariou")
  #   group.by = "Gender"
  #   group.subset = c()
  #   set.seed(821)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  #   
  # }else if (dataset == "Chariou.subgroup"){
  #   
  #   print("selected Chariou dataset") 
  #   object = selectCRObject("Chariou")
  #   group.by = "Gender"
  #   group.subset = c("M")
  #   set.seed(822)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  #   
  # } else if (dataset == "pbmc.single"){
  #   
  #   print("selected nsclc_single dataset") 
  #   object = selectCRObject("pbmc-single")
  #   group.by = "orig.ident"
  #   group.subset =  unique(object$orig.ident)
  #   set.seed(83)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  #   
  # } else if (dataset == "nsclc.multi"){
  #   
  #   print("selected nsclc_multi dataset") 
  #   object = selectCRObject("nsclc-multi")
  #   group.by = "orig.ident"
  #   group.subset =  unique(object$orig.ident)[1:3]
  #   set.seed(84)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  # } else if (dataset == "brca"){
  #   
  #   print("selected BRCA dataset") 
  #   object = selectCRObject("BRCA")
  #   group.by = "orig.ident"
  #   group.subset =  unique(object$orig.ident)[1:3]
  #   set.seed(85)
  #   genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
  #                              replace = FALSE)
  #   
  # }
  object = selectCRObject("TEC")
  group = "orig_ident"
  assay = 'SCT'
  slot = 'scale.data'
  jitter_points = T
  jitter_dot_size = 4
  
  return(list("object" = object, 
              "group" = group,
              "assay" = assay,
              "slot" = slot,
              "jitter_points" = jitter_points,
              "jitter_dot_size" = jitter_dot_size,
              "genes" = genes))}

.drawViolin <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
