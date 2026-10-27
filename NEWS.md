# SCWorkflow (development version)

## Dependencies

* Declare `Matrix`, `shiny`, and `patchwork` as required dependencies (#107).
* Declare optional analysis dependencies already used by the package, including trajectory-analysis libraries (#107).
* Require `ggalluvial` for cell-population comparison plots (#107).

## Fixes

* Resolve package namespace conflicts and correct Seurat cell-cycle data references (#107).

## compareCellPopulations

* Add `compareCellPopulations()` to compare cell population frequencies or counts across experimental groups using alluvial and box plots, with customizable group order and colors (#107).

## plotMetadata

* Improved metadata selection, validation, and run logging for `plotMetadata()` (#95).
* Improve sample-selection guidance and log selected sample names and plotted/input cell counts (#106).
* Keep cells with missing metadata in `plotMetadata()` plots, colored grey with an `NA` legend entry, including summarized metadata (#106).
* Validate requested RNA and CITE-seq reductions and log the resolved reduction after sample selection and cell counts (#106).
* Respect `dot.size` for numeric metadata plots as well as categorical plots (#106).
* Warn about missing requested metadata columns before listing available columns and continue plotting valid selections; empty or wholly invalid selections still error (#106).
* Name returned metadata plots by their retained columns to prevent mislabeling when missing selections are skipped (#106).
* Correct the documented list return value and clarify that sample selection and metadata summarization affect plots only (#107).

## reclusterSeuratObject

* Preserve requested metadata columns by their exact names during reclustering (#96).
* Clarify validation and preservation of prior clustering columns during reclustering (#104).
* Improve detection, preservation, and guidance for prior clustering columns during reclustering (#105).

## heatmapSC

* Add `sample.name.column` to select the metadata column used for sample filtering, defaulting to `orig.ident` (#100).
* Log requested sample filtering and displayed filtered/total cell counts (#101).
* Use all samples when sample selection is blank (#102).
* Report found Heatmap features before missing requested features without duplicate warnings (#103).
* Identify unavailable requested samples and list available sample names.

## Documentation

* Complete help for `filterQC()` and module-score plotting, and correct population-comparison documentation (#107).
* Correct vignette configuration in source packages (#107).

# SCWorkflow 1.0.1 (2024-02-01)

## Features

* Add the initial filtering and quality-control workflow ([a843d4d](https://github.com/NIDAP-Community/SCWorkflow/commit/a843d4dd33ea83edabb71f441e2cfa04788b67ba)).
* Add heatmap, dotplot, 3D t-SNE, and dual-labeling functions ([0ed4df7](https://github.com/NIDAP-Community/SCWorkflow/commit/0ed4df7b93dc0ffdf2c7baabcda82be769852e7f)).
* Add cluster naming ([dbe4d15](https://github.com/NIDAP-Community/SCWorkflow/commit/dbe4d157527f9ad86080d32dd035e83759ad3f95)), sample naming ([9f16588](https://github.com/NIDAP-Community/SCWorkflow/commit/9f165886a1cf416db814fcd88d5487ee86f7cd1e)), and metadata tables ([42c3569](https://github.com/NIDAP-Community/SCWorkflow/commit/42c3569db6e71563a56387ac78876dcd4a42a347)).
* Introduce the `cluster` variable ([e26f2aa](https://github.com/NIDAP-Community/SCWorkflow/commit/e26f2aaf7d79a2ff57e202758739bd40b69c8460)).
* Add color selection for 3D t-SNE plots ([2882b76](https://github.com/NIDAP-Community/SCWorkflow/commit/2882b768758596c7233f73235eafcc4e8cc3c41e)).
* Add a cached celldex option ([2cabb8d](https://github.com/NIDAP-Community/SCWorkflow/commit/2cabb8dfcc1b0f2567efc92c6d0d2d77a621d021)).
* Add a contingency-table plot to dual labeling ([0b17248](https://github.com/NIDAP-Community/SCWorkflow/commit/0b17248be847b81555c97600fe8e2caff9db7d6b)).

## Fixes

* Fix heatmap plotting and sorting by annotation ([db3ee4e](https://github.com/NIDAP-Community/SCWorkflow/commit/db3ee4e74d8118df4e8fdd77aa3486200e7f4adf), [8d09781](https://github.com/NIDAP-Community/SCWorkflow/commit/8d097818a0ebe01be9fd9b077a451663c7af9ed9)).
* Fix dual labeling ([9eede5b](https://github.com/NIDAP-Community/SCWorkflow/commit/9eede5bc4f310ea22207cec057c617f16202b6d4)).
* Update violin-plot subsetting ([0609899](https://github.com/NIDAP-Community/SCWorkflow/commit/06098990cda9235bec91260bea697405304cefb9)).
* Use barcodes rather than clusters for merging ([f561e2f](https://github.com/NIDAP-Community/SCWorkflow/commit/f561e2fe5c425efa25636f817943a0156427019d)).
* Fix the `All` cluster label ([18bbdc8](https://github.com/NIDAP-Community/SCWorkflow/commit/18bbdc841e98fc8e28d6a3fcfac5e3b7255c47d4)).
* Rename `select_crobject` to `selectCRObject` ([b7446a6](https://github.com/NIDAP-Community/SCWorkflow/commit/b7446a68cd3fef6e06f1fa1ea2cddc48669298b4)).
* Update Seurat importing in raw-data processing ([846ac5d](https://github.com/NIDAP-Community/SCWorkflow/commit/846ac5d5b71ae0d188bd96d50846d9210c3cd453)).
* Update latent-variable and `second.clust` handling ([67c3aea](https://github.com/NIDAP-Community/SCWorkflow/commit/67c3aea82d8775e61954cf1135973b4ad3b21f35)).
* Use `ComplexHeatmap::pheatmap` for heatmaps ([5c9316b](https://github.com/NIDAP-Community/SCWorkflow/commit/5c9316bf95e1779a6ef56f9e6d41a97d8a56e711)).
* Improve dotplot error messages ([9742c5f](https://github.com/NIDAP-Community/SCWorkflow/commit/9742c5f72f93f3f0a802a7aa7c13988f51e25752)).
* Update cell-cycle gene calls ([797aca1](https://github.com/NIDAP-Community/SCWorkflow/commit/797aca16c7da61a0af131138f90075e88c3f69f8)).
* Adjust PNG rendering for NIDAP display ([0b9e4d2](https://github.com/NIDAP-Community/SCWorkflow/commit/0b9e4d28e2e764c85135caada0075521320549ed)).
* Fix the "after" plot visualization ([d7bb90f](https://github.com/NIDAP-Community/SCWorkflow/commit/d7bb90f5a7cb674ffb4605a9378b55e4e2d9f48e)).
* Update package function loading ([fe7ad65](https://github.com/NIDAP-Community/SCWorkflow/commit/fe7ad65efba87c382e91d1adfa6cbace3b707ace)).
* Update metadata plotting and module scoring ([573b57a](https://github.com/NIDAP-Community/SCWorkflow/commit/573b57a227bed380ccbeb7fe5e48b84a21627c72), [4367b28](https://github.com/NIDAP-Community/SCWorkflow/commit/4367b28db3b71aa627b4f3e33dfda52dff6cf6de)).

## Documentation

* Add the original pseudobulk citation ([4293f96](https://github.com/NIDAP-Community/SCWorkflow/commit/4293f96b1ba14693905e854085270919b67eaa6c)).
