# SCWorkflow (development version)

## Dependencies

* Align Gitflow validation with the frozen R 4.1.3 runtime and isolated CI-tool library while preserving existing deployment configuration (#107).
* Declare optional trajectory dependencies already used by the package (#107).
* Declare `patchwork` in `Imports` to match its existing namespace import (#107).
* Keep `ggalluvial` optional so SCWorkflow can load without it, and require it only when comparing cell populations (#107).
* Validate package changes with the preserved R 4.1.3 environment and a separate checksum-pinned CI-tool library, without updating existing dependencies (#107).

## plotMetadata

* Improved metadata selection, validation, and run logging for `plotMetadata()` (#95).
* Improve sample-selection guidance and log selected sample names and plotted/input cell counts (#106).
* Keep cells with missing metadata in `plotMetadata()` plots, colored grey with an `NA` legend entry, including summarized metadata (#106).
* Validate requested RNA and CITE-seq reductions and log the resolved reduction after sample selection and cell counts (#106).
* Respect `dot.size` for numeric metadata plots as well as categorical plots (#106).
* Warn about missing requested metadata columns before listing available columns and continue plotting valid selections; empty or wholly invalid selections still error (#106).
* Name returned metadata plots by their retained columns to prevent mislabeling when missing selections are skipped (#106).

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
