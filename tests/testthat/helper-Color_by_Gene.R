#Returning parameters for each test depending on the data input:

getParamCBG <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    samples.to.include <-
      'c("1_Embryo_13_5",
         "2_Embryo_15",
         "3_Newborn",
         "4_Adult")'
    gene <- c("Gapdh",
              "GAPDH",
              "Foxp3",
              "Il2ra",
              "Cd4",
              "Cd8a")
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    samples.to.include = 'c("CD8dep",
                             "Combo",
                             "ENT",
                             "NHSIL12",
                             "PBS")'
    gene = c("Gapdh",
             "GAPDH",
             "Foxp3",
             "Il2ra",
             "Cd4",
             "Cd8a")
  } else if (data == "NSCLC_Single") {
    object <- selectCRObject("NSCLC_Single")
    samples.to.include = paste0('c("', unique(as.character(object$orig.ident)), '")')
    gene = c("Gapdh",
             "GAPDH",
             "FOXP3",
             "IL2RA",
             "CD4",
             "CD8A",
             "MYC")
  } else if (data == "nsclc-multi") {
    object <- selectCRObject("nsclc-multi")
    samples.to.include = 'c("Donor_1",
                            "Donor_2",
                            "Donor_3",
                            "Donor_4",
                            "Donor_5",
                            "Donor_6",
                            "Donor_7")'
    gene = c("Gapdh",
             "GAPDH",
             "FOXP3",
             "IL2RA",
             "CD4",
             "CD8A",
             "MYC",
             "AMER2",
             "CD9",
             "CEBPD")
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
    samples.to.include = 'c("CID3586",
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
    gene = c("Gapdh",
             "GAPDH",
             "FOXP3",
             "IL2RA",
             "CD4",
             "CD8A",
             "MYC",
             "PARP8")
  }
  
  
  
  
  return(list(
    "object" = object,
    "samples.to.include" = samples.to.include,
    "gene" = gene
  ))
}
