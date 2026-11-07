#This helper script will return parameters for each test dependent on the input:

getParamDGEM <- function(data){
  if(data == "TEC"){
      object <- selectCRObject("TEC")
      samples = 'c("1_Embryo_13_5","2_Embryo_15")'
      parameter.to.test = "SCT_snn_res_0_2"
      contrasts =  c("1-2","2-all")
  } else if (data == "Chariou") {
      object <- selectCRObject("Chariou")
      samples = 'c("CD8dep","Combo","ENT","NHSIL12","PBS")'
      parameter.to.test = "SCT_snn_res_2_8"
      contrasts =  c("0-1","0-all")
  } else if (data == "NSCLC_Single") {
      object <- selectCRObject("NSCLC_Single")
      samples = paste0('c("', unique(as.character(object$orig.ident)), '")')
      parameter.to.test = "SCT_snn_res_0_2"
      contrasts =  c("0-1","0-all")
  } else if (data == "nsclc-multi") {
      object <- selectCRObject("nsclc-multi")
      samples = 'c("Donor_1",
                   "Donor_2",
                   "Donor_3",
                   "Donor_4",
                   "Donor_5",
                   "Donor_6",
                   "Donor_7")'
      parameter.to.test = "SCT_snn_res_0_2"
      contrasts =  c("0-1","0-all")
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      samples = 'c("CID3586",
                   "CID3838",
                   "CID3921",
                   "CID3941",
                   "CID3946",
                   "CID3948",
                   "CID3963",
                   "CID4040",
                   "CID4066",
                   "CID4067",
                   "CID4290A",
                   "CID4398",
                   "CID44041",
                   "CID4461",
                   "CID4463",
                   "CID4465",
                   "CID4471",
                   "CID4495",
                   "CID44971",
                   "CID44991",
                   "CID4513",
                   "CID4515",
                   "CID45171",
                   "CID4523",
                   "CID4530N",
                   "CID4535")'
      parameter.to.test = "SCT_snn_res_0_2"
      contrasts =  c("0-1","0-all")
  }

  return(list("object" = object,
              "samples"=samples, 
              "parameter.to.test" = parameter.to.test,
              "contrasts" = contrasts))  
}