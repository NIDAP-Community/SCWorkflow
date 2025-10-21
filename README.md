# SCWorkflow

[![Gitflow Action for R Package Development](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml/badge.svg)](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml) 
[![Version](https://img.shields.io/github/v/release/nidap-community/scworkflow)](https://github.com/NIDAP-Community/SCWorkflow/releases/latest)
[![Docker Image Version](https://img.shields.io/docker/v/nciccbr/scworkflow?label=docker)](https://hub.docker.com/r/nciccbr/scworkflow)



[![](https://raw.githubusercontent.com/NIDAP-Community/SCWorkflow/GalaxyCLI/vignettes/SCWorkflow.png)](https://lucid.app/lucidchart/c7b852ad-72dc-4821-90d5-e45bed0c4199/view)


<center>**Click Figure to Navigate Workflow**.</center>

<br>
<br>

The Single Cell Workflow streamlines the analysis of multimodal Single Cell RNA-Seq data produced from 10x Genomics.  It can be run in a docker container, and for biologists, in user-friendly web-based interactive notebooks (NIDAP, Palantir Foundry). Much of it is based on the Seurat workflow in Bioconductor, and supports CITE-Seq data.  It incorporates a cell identification step (ModScore) that utilizes module scores obtained from Seurat and also includes Harmony for batch correction.


For further documentation see our detailed [Docs Website](https://nidap-community.github.io/SCWorkflow/)





## Install Package

```
# install.packages("remotes")
# remotes::install_github("NIDAP-Community/SCWorkflow", dependencies = TRUE)

library(SCWorkflow)
```
<br>
Future Developments include addition of support for multiomics (TCR-Seq, ATAC-Seq) single cell data and integration with spatial transcriptomics data.

