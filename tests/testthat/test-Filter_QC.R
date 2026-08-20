expect_objects_have_counts <- function(objects) {
  count.dims <- vapply(objects, function(so) {
    counts <- so@assays$RNA@counts
    c(
      rows = nrow(counts),
      columns = ncol(counts),
      nonzero = Matrix::nnzero(counts),
      size = as.numeric(object.size(counts))
    )
  }, numeric(4))

  expect_true(all(count.dims["rows", ] > 0))
  expect_true(all(count.dims["columns", ] > 0))
  expect_true(all(count.dims["nonzero", ] > 0))
  expect_true(all(count.dims["size", ] > 0))
}

for (data in c('TEC','Chariou','PBMC_Single','NSCLC_Multi')) {
  
  test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","data","plots")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(data.run$object))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plots$PostFilterCombined)[2], 'ggplot')
    # SO slot contains data
    expect_objects_have_counts(filter.qc.out$object)
    # plot slot contains data
    expect_true(object.size(filter.qc.out$plots) > 0)
    
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(filter.qc.out$plots$PostFilterCombined),
      paste0(data,"_Standard_combFig.png")
    )
    
    
  })
  
}


################################################################

for (data in c('Chariou')) {
  # data='   } else if (data == "pbmc-single'
  
  test_that(paste0("Test Filter and QC - filter.vdj.genes (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    data.run$filter.vdj.genes <- "TRUE"
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","data","plots")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(data.run$object))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plots$PostFilterCombined)[2], 'ggplot')
    # SO slot contains data
    expect_objects_have_counts(filter.qc.out$object)
    # plot slot contains data
    expect_true(object.size(filter.qc.out$plots) > 0)
    # Check if VDJ genes are removed
    expect_equal(
     sum(filter.qc.out$data$FilteringCounts$`VDJ Genes Removed` > 0),
     5)

    
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(filter.qc.out$plots$PostFilterCombined),
      paste0(data,"_filter.vdj.genes_combFig.png")
    )
    
    
  })
  
}




