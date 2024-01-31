getParamRaw <- function(data){
  
  if(data == "TEC"){
    
    data.table(Sample_Name=c('5_ABSC_E13_5_TEC','6_ABSC_Newborn_cTEC',
                              '7_ABSC_Adult_cTEC','8_ABSC_E15_cTEC'),
               Gender=c('M','F','N','F'),
               Rename=c("S1_Embryo_13_5","S3_Newborn",
                        "S4_Adult","S2_Embryo_15"))%>%
      write.table(
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
      sep = '\t')
    
    
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    split.h5=F
    

  
    } else if (data == "Chariou") {
      
      data.table(Sample_Name=c('SCAF1713_1_1','SCAF1714_2_1','SCAF1715_3_1',
                                'SCAF1716_4_1','SCAF1717_5_1'),
                 Rename=c("PBS","ENT","NHSIL12","Combo","CD8dep"))%>%
        write.table(
          test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
          sep = '\t')
      
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),"",full.names = T)
    # input=list.files(
    #   test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    split.h5=F
    
    
    
  } else if (data == "NSCLC_Single") {
    
    data.table(Sample_Name=c('20k_PBMC_3p_HT_nextgem_Chromium_X'),
               Rename=c("PBMC_Single"))%>%
      write.table(
        test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
        sep = '\t')
    
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    split.h5=F
    
    
  } else if (data == "NSCLC_Multi") {
  
      data.table(Sample_Name=
        c('NSCLC_40k_DTC_3p_HT_nextgem_donor_1_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_2_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_3_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_4_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_5_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_6_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_7_count_sample_feature_bc_matrix'),
                 Rename=c("Donor_1","Donor_2",
                          "Donor_3","Donor_4",
                          "Donor_5","Donor_6","Donor_7"))%>%
      write.table(
        test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
        sep = '\t')      
    input=list.files(  
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    split.h5=F
    
    
    
    
  } else if (data == "BRCA") {
    h5= '/rstudio-files/ccbr-data/data/singlecell/BRCA/NG.h5'
    fixture=test_path(paste0("fixtures/",data,"/h5files/NG.h5"))
    if(file.exists(fixture)==F){
      file.copy(h5,fixture)
    }
    
    data.table(My_Sample_Names=
                c('CID3586','CID3921','CID45171','CID3838','CID4066','CID44041',
                  'CID4465','CID4495','CID44971','CID44991','CID4513','CID4515',
                  'CID4523','CID3946','CID3963','CID4461','CID4463','CID4471',
                  'CID4530N','CID4535','CID4040','CID3941','CID3948','CID4067',
                  'CID4290A','CID4398'),
               My_Renamed_Samples=
                c('CID3586','CID3921','CID45171','CID3838','CID4066','CID44041',
                  'CID4465','CID4495','CID44971','CID44991','CID4513','CID4515',
                  'CID4523','CID3946','CID3963','CID4461','CID4463','CID4471',
                  'CID4530N','CID4535','CID4040','CID3941','CID3948','CID4067',
                  'CID4290A','CID4398'),
               My_Variable_1=
                 c('M','F','M','F','F','F','M','M','F','M','M','M','M','F','F',
                   'M','F','F','M','F','F','F','M','M','F',NA),
               My_Variable_2=
                 c(1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,NA)
               )%>%
      write.table(
        test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
        sep = '\t')      
    
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='My_Sample_Names'
    rename.col='My_Renamed_Samples'
    split.h5=T
    
    
  }
  
  return(list("input" = input, 
              "organism"=organism, 
              "sample.metadata.table"=sample.metadata.table,
              'sample.name.column'=sample.name.column,
              "rename.col"=rename.col,
              "split.h5"=split.h5
              
              ))  
}


.drawFig <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
.saveSO <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".rds")
  saveRDS(x, file = path)
  print(path)
}
