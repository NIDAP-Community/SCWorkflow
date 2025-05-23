#This helper script will return parameters for each test dependent on the input:

getParamPM <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    samples.to.include <-
      'c("1_Embryo_13_5",
      "2_Embryo_15",
      "3_Newborn",
      "4_Adult")'
    metadata.to.plot <-
      'c("SCT_snn_res_0_2",
         "SCT_snn_res_0_4",
         "SCT_snn_res_0_6",
         "SCT_snn_res_0_8",
         "SCT_snn_res_1",
         "SCT_snn_res_1_2")'
    columns.to.summarize = "c()"
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    samples.to.include = 'c("CD8dep",
                            "Combo",
                            "ENT",
                            "NHSIL12",
                            "PBS")'
    metadata.to.plot = 'c("SCT_snn_res_2_4",
                          "SCT_snn_res_2_6",
                          "SCT_snn_res_2_8")'
    columns.to.summarize = "c()"
  } else if (data == "pbmc-single") {
    object <- selectCRObject("pbmc-single")
    samples.to.include = 'c("PBMC_Single")'
    metadata.to.plot = 'c("SCT_snn_res_0_2",
                          "SCT_snn_res_0_4",
                          "SCT_snn_res_0_6",
                          "SCT_snn_res_0_8",
                          "SCT_snn_res_1",
                          "SCT_snn_res_1_2")'
    columns.to.summarize = "c()"
  } else if (data == "nsclc-multi") {
    object <- selectCRObject("nsclc-multi")
    samples.to.include = 'c("Donor_1",
                            "Donor_2",
                            "Donor_3",
                            "Donor_4",
                            "Donor_5",
                            "Donor_6",
                            "Donor_7")'
    metadata.to.plot = 'c("SCT_snn_res_0_2",
                          "SCT_snn_res_0_4",
                          "SCT_snn_res_0_6",
                          "SCT_snn_res_0_8",
                          "SCT_snn_res_1",
                          "SCT_snn_res_1_2")'
    columns.to.summarize = "c()"
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
    metadata.to.plot = 'c("SCT_snn_res_0_2",
                          "SCT_snn_res_0_4",
                          "SCT_snn_res_0_6",
                          "SCT_snn_res_0_8",
                          "SCT_snn_res_1",
                          "SCT_snn_res_1_2")'
    columns.to.summarize = "c()"
  }
  
  return(
    list(
      "object" = object,
      "samples.to.include" = samples.to.include,
      "metadata.to.plot" = metadata.to.plot,
      "columns.to.summarize" = columns.to.summarize
    )
  )
}
