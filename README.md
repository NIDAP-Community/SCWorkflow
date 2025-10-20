# SCWorkflow

Workflow Package for Analysis of Single Cell Data

```mermaid
flowchart TD
    C["Import Data & Quality Control"]
    D["Annotate Cells"]
     C:::Peach
     D:::Sky
    classDef Peach stroke-width:1px, stroke-dasharray:none, stroke:#FBB35A, fill:#FFEFDB, color:#8F632D
    classDef Sky stroke-width:1px, stroke-dasharray:none, stroke:#374D7C, fill:#E2EBFF, color:#374D7C
    click C "https://github.com/NIDAP-Community/SCWorkflow/blob/GalaxyCLI/vignettes/Getting%20Started%20and%20Quality%20Control.html"
```


[![Gitflow Action for R Package Development](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml/badge.svg)](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml)
[![Version](https://img.shields.io/github/v/release/nidap-community/scworkflow)](https://github.com/NIDAP-Community/SCWorkflow/releases/latest)
[![Docker Image Version](https://img.shields.io/docker/v/nciccbr/scworkflow?label=docker)](https://hub.docker.com/r/nciccbr/scworkflow)

The Single Cell Workflow streamlines the analysis of multimodal Single Cell RNA-Seq data produced from 10x Genomics.  It can be run in a docker container, and for biologists, in user-friendly web-based interactive notebooks (NIDAP, Palantir Foundry). Much of it is based on the Seurat workflow in Bioconductor, and supports CITE-Seq data.  It incorporates a cell identification step (ModScore) that utilizes module scores obtained from Seurat and also includes Harmony for batch correction.

Some of the steps in the workflow:

<img src="scWorkflow_image.png">


Future Developments include addition of support for multiomics (TCR-Seq, ATAC-Seq) single cell data and integration with spatial transcriptomics data.

