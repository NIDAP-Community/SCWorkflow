for (data in c('TEC','Chariou','PBMC_Single','PBMC_Multi')) {
  
  test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","FilteringMeta","plots")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(data.run$object))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plots$PostFilterCombined)[2], 'ggplot')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plots),'> 0' )
    
  })
  
}



################################################################

for (data in c('TEC')) {
  # data='   } else if (data == "pbmc-single'
  
  test_that(paste0("Test Filter and QC - filter.vdj.genes (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    data.run$filter.vdj.genes <- "TRUE"
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","FilteringMeta","plots")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(data.run$object))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plots$PostFilterCombined)[2], 'ggplot')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plots),'> 0' )
    
  })
  
}




