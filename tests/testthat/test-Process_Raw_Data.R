
for (data in c('TEC','NSCLC_Multi','NSCLC_Single')) {
  
  test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
    
    
    data.run <- getParamRaw(data)
    Raw.out <- do.call(processRawData, data.run)

    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(Raw.out$object),length(data.run$input))
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'= 0' )
    
    
    # Check for Identical files
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(Raw.out$plots$CombinedQC),
      paste0(data,"_Standard_combFig.png")
    )
    # expect_snapshot_file( # Test failed each run with no changes
    #   .saveSO(Raw.out$object),
    #   paste0(data,"_Standard.rds")
    # )
    
    
  })
  
}

for (data in c('BRCA')) {
  
  test_that(paste0("Test Split h5 (",data," dataset)"), {
    
    
    data.run <- getParamRaw(data)
    Raw.out <- do.call(processRawData, data.run)
    
    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO is Split into multiple files
    expect_false(isTRUE(all.equal(length(Raw.out$object),
                                  length(data.run$input))))
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'= 0' )
    
    # Check for Identical files
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(Raw.out$plots$CombinedQC),
      paste0(data,"_Standard_combFig.png")
    )
    # expect_snapshot_file( # Test failed each run with no changes
    #   .saveSO(Raw.out$object),
    #   paste0(data,"_Standard.rds")
    # )
    
  })
  
}


################################################################
################################################################

for (data in c('Chariou')) {
  
  test_that(paste0("Test Filter and QC - TCR data for (",data," dataset)"), {
    
    
    data.run <- getParamRaw(data)
    Raw.out <- do.call(processRawData, data.run)
    
    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(Raw.out$object),length(grep('\\.h5',data.run$input,value = T)))
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'= 0' )
    
    
    TCRmeta=c("ab_pair","cell_beta_seq_list","cell_beta_reads_list",
              "cell_unique_betas","cell_TRBV_list",
              "cell_TRBJ_list","cell_alpha_seq_list","cell_alpha_reads_list",
              "cell_unique_alphas", "isPolyAlphaCell",
              "isPolyBetaCell","cell_top_beta","cell_TRBV",
              "cell_TRBJ","cell_top_alpha",
              "cell_TRAV","cell_TRAJ","clonotype_id",
              "summarized_cell_top_alpha", "summarized_cell_top_beta",
              "summarized_clonotype_id")
    for (n in names(Raw.out$object)) {
      so=Raw.out$object[[n]]
      # check if TCR metadata columns have been added
      expect_true(all(TCRmeta%in%colnames(so@meta.data)),
                  label = paste0("Sample: ",n,
                                 " Metadata table not updated with TCR data"))
      # Check if TCR Metadata Columns contain data
      
      TCRtbl=so@meta.data[!is.na(so@meta.data[,TCRmeta[1]]),TCRmeta,drop=F]
      # Check if TCR data is added 
      expect_gt(nrow(TCRtbl),0,
             label = paste0("Sample: ",n, " TCR data Not Added"))
      # Check if TCR data is added 
      expect_equal(nrow(TCRtbl[rowSums(is.na(TCRtbl))>0,]),0 ,
                   label = paste0("Sample: ",n,
                                  " TCR data contains missing values"))

    }
    
    # Check for Identical files
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(Raw.out$plots$CombinedQC),
      paste0(data,"_Standard_combFig.png")
    )
    # expect_snapshot_file( # Test failed each run with no changes
    #   .saveSO(Raw.out$object),
    #   paste0(data,"_Standard.rds")
    # )
    
  })
  
}



################################################################

for (data in c('TEC')) {
  
  test_that(paste0("Test Filter and QC - Regex (",data," dataset)"), {
    
    data.run <- getParamRaw(data)
    data.run$keep = TRUE
    data.run$file.filter.regex = "1|2"
    Raw.out <- do.call(processRawData, data.run)
    
    
    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(Raw.out$object),2)
    # test that samples 
    
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'>0' )
    

    
    # Check for Identical files
    skip_on_ci()
    expect_snapshot_file(
      .drawFig(Raw.out$plots$CombinedQC),
      paste0(data,"_Regex_combFig.png")
    )
    # expect_snapshot_file( # Test failed each run with no changes
    #   .saveSO(Raw.out$object),
    #   paste0(data,"_Standard.rds")
    # )
    
    
  })
  
}


################################################################

# for (data in c('TEC')) {
#   
#   test_that(paste0("Test Filter and QC - Test organsim ERROR:  (",data," dataset)"), {
#     
#     data.run <- getParamRaw(data)
#     data.run$organism = "Human"
#     Raw.out <- do.call(processRawData, data.run)
#     
#     expect_error(do.call(processRawData, data.run), 
#           c("No Mitochondrial Genes Detetcted: Wrong Organism may be selected.
#            Supported Organisms are Human or Mouse"), ignore.case = TRUE)
#     
#   })
#   
# }



