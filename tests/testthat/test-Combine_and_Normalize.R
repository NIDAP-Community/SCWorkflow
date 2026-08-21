
for (data in c(
  'TEC','Chariou',
  'NSCLC_Single')) {

  test_that(
    paste0("Test Combine & Renormalize - Standard (",data," dataset)"), {


      data.run <- getParamCN(data)
      if(length(data.run$object)>2){
      data.run$object=data.run$object[c(1,2)]}
      combine.renormalize.out <- do.call(combineNormalize, data.run)


      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
      # plot slot contains data
      expect_true(object.size(combine.renormalize.out$plots) > 0)

      skip_on_ci()
      expect_snapshot_file(
        .drawFig(combine.renormalize.out$plots$TSNE),
        paste0(data,"_Standard_TSNE.png")
      )
      ## Identical test runs returns different rds files
      # expect_snapshot_file(
      #   .saveSO(combine.renormalize.out$object),
      #   paste0(data,"_Standard.rds")
      # )


    })
}


for (data in c('NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - Standard (",data," dataset)"), {
      
      
      data.run <- getParamCN(data)
      data.run$object=data.run$object[1:2]
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # Data contains Protein Slot
      expect_setequal(names(combine.renormalize.out$object@assays), 
                      c("RNA","Protein","SCT"))
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
      # plot slot contains data
      expect_true(object.size(combine.renormalize.out$plots) > 0)
      expect_true(object.size(combine.renormalize.out$plots$`CITEseq UMAP`) > 0)
      
      skip_on_ci()
      expect_snapshot_file(
        .drawFig(combine.renormalize.out$plots$TSNE),
        paste0(data,"_Standard_TSNE.png")
      )
      ## Identical test runs returns different rds files
      # expect_snapshot_file(
      #   .saveSO(combine.renormalize.out$object),
      #   paste0(data,"_Standard.rds")
      # )
      
      
    })
}



for (data in c('TEC')) {

  test_that(
    paste0("Test Combine & Renormalize - only.var.genes (",data," dataset)"), {


      data.run <- getParamCN(data)
      data.run$object=data.run$object[c(2,3)]
      data.run$only.var.genes = TRUE
      combine.renormalize.out <- do.call(combineNormalize, data.run)


      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
      # plot slot contains data
      expect_true(object.size(combine.renormalize.out$plots) > 0)

      skip_on_ci()
      expect_snapshot_file(
        .drawFig(combine.renormalize.out$plots$TSNE),
        paste0(data,"_only.var.genes_TSNE.png")
      )
      # expect_snapshot_file(
      #   .saveSO(combine.renormalize.out$object),
      #   paste0(data,"_only.var.genes.rds")
      # )

    })
}



for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - SCT level - sample (",data," dataset)"), {
      
      
      data.run <- getParamCN(data)
      data.run$object=data.run$object[c(2,3)]
      data.run$SCT.level = 'Sample'
      data.run$object = data.run$object
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
      # plot slot contains data
      expect_true(object.size(combine.renormalize.out$plots) > 0)

      skip_on_ci()
      expect_snapshot_file(
        .drawFig(combine.renormalize.out$plots$TSNE),
        paste0(data,"_SCTlevel_TSNE.png")
      )

      # expect_snapshot_file(
      #   .saveSO(combine.renormalize.out$object),
      #   paste0(data,"_SCTlevel.rds")
      # )
      
    })
} 


for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - exclude.sample (",data," dataset)"), {
      
      data.run <- getParamCN(data)
      data.run$object=data.run$object[c(2,3)]
      data.run$exclude.sample=data.run$object[1]%>%names
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
      # plot slot contains data
      expect_true(object.size(combine.renormalize.out$plots) > 0)
      # sample is removed
      expect_false(
        data.run$exclude.sample%in%
          unique(combine.renormalize.out$object@meta.data$orig.ident)
                  )
      
      skip_on_ci()
      expect_snapshot_file(
        .drawFig(combine.renormalize.out$plots$TSNE),
        paste0(data,"_exclude.sample_TSNE.png")
      )
      # expect_snapshot_file(
      #   .saveSO(combine.renormalize.out$object),
      #   paste0(data,"_exclude.sample.rds")
      # )
      # 
    })
}


for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - selection.method = mean.var.plot (",
           data," dataset)"), {
             
             
             data.run <- getParamCN(data)
             data.run$object=data.run$object[c(2,3)]
             data.run$selection.method = "mean.var.plot"
             combine.renormalize.out <- do.call(combineNormalize, data.run)
             
             
             # create output
             expected.elements = c("object","plots")
             expect_setequal(names(combine.renormalize.out), expected.elements)
             # figure slot is a grob
             expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
             # SO slot contains data
             expect_true(nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
             # plot slot contains data
             expect_true(object.size(combine.renormalize.out$plots) > 0)
             
             skip_on_ci()
             expect_snapshot_file(
               .drawFig(combine.renormalize.out$plots$TSNE),
               paste0(data,"_selection.method_TSNE.png")
             )
             # expect_snapshot_file(
             #   .saveSO(combine.renormalize.out$object),
             #   paste0(data,"_selection.method.rds")
             # )
             
           })
}


for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - Select PCs plots"), {
             
             
             data.run <- getParamCN(data)
             data.run$object=data.run$object[c(2,3)]
             data.run$methods.pca = c('Elbow','Marchenko-Pastur')
             data.run$pca.reg.plot=T
             data.run$vars.to.regress = c("percent.mt")
             data.run$jackstraw = T
             data.run$jackstraw.dims=2
             
             combine.renormalize.out <- do.call(combineNormalize, data.run)
             
             
             # create output
             expected.elements = c("object","plots")
             expect_setequal(names(combine.renormalize.out), expected.elements)
             # figure slot is a grob
             expect_equal(
               class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
             # SO slot contains data
             expect_true(
               nrow(combine.renormalize.out$object@assays$RNA@counts) > 0)
             # plot slot contains data
             expect_true(
               object.size(combine.renormalize.out$plots$TSNE) > 0)
             expect_true(
               object.size(combine.renormalize.out$plots$UMAP) > 0)
             expect_true(
               object.size(combine.renormalize.out$plots$`Elbow Plot`) > 0)
             expect_true(object.size(
               combine.renormalize.out$plots$`Regression Plots`[[1]]) > 0)
             expect_true(
               object.size(combine.renormalize.out$plots$JackStraw) > 0)
             
             skip_on_ci()
             expect_snapshot_file(
               .drawFig(combine.renormalize.out$plots$TSNE),
               paste0(data,"_SelectPCsPlot_TSNE.png")
             )
             expect_snapshot_file(
               .drawFig(combine.renormalize.out$plots$`Elbow Plot`),
               paste0(data,"_SelectPCsPlot_Elbow.png")
             )
             expect_snapshot_file(
               .drawFig(combine.renormalize.out$plots$JackStraw),
               paste0(data,"_SelectPCsPlot_JackStraw.png")
             )
             expect_snapshot_file(
               .drawFig(combine.renormalize.out$plots$`Regression Plots`[[1]]),
               paste0(data,"_SelectPCsPlot_PCAreg.png")
             )
             # expect_snapshot_file(
             #   .saveSO(combine.renormalize.out$object),
             #   paste0(data,"_selection.method.rds")
             # )
             
           })
}




