# SCWorkflow (development version)

## plotMetadata

* Improved metadata selection, validation, and run logging for `plotMetadata()` (#95).
* Keep cells with missing metadata in `plotMetadata()` plots, colored grey with an `NA` legend entry, including summarized metadata.
* Validate requested RNA and CITE-seq reductions and log the resolved reduction after sample selection and cell counts.
* Respect `dot.size` for numeric metadata plots as well as categorical plots.

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
