# SCWorkflow

Workflow Package for Analysis of Single Cell Data

[![Gitflow Action for R Package Development](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml/badge.svg)](https://github.com/NIDAP-Community/SCWorkflow/actions/workflows/gitflow-R-action.yml)
[![Version](https://img.shields.io/github/v/release/nidap-community/scworkflow)](https://github.com/NIDAP-Community/SCWorkflow/releases/latest)
[![Docker Image Version](https://img.shields.io/docker/v/nciccbr/scworkflow?label=docker)](https://hub.docker.com/r/nciccbr/scworkflow)

The Single Cell Workflow streamlines the analysis of multimodal Single Cell RNA-Seq data produced from 10x Genomics.  It can be run in a docker container, and for biologists, in user-friendly web-based interactive notebooks (NIDAP, Palantir Foundry). Much of it is based on the Seurat workflow in Bioconductor, and supports CITE-Seq data.  It incorporates a cell identification step (ModScore) that utilizes module scores obtained from Seurat and also includes Harmony for batch correction.

Some of the steps in the workflow:

<img src="scWorkflow_image.png">


Future Developments include addition of support for multiomics (TCR-Seq, ATAC-Seq) single cell data and integration with spatial transcriptomics data.

