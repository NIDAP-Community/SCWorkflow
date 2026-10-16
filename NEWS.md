# SCWorkflow (development version)

## plotMetadata

* Improved metadata selection, validation, and run logging for `plotMetadata()` (#95).

## reclusterSeuratObject

* Preserve requested metadata columns by their exact names during reclustering (#96).

## heatmapSC

* Add `sample.name.column` to select the metadata column used for sample filtering, defaulting to `orig.ident` (#100).
* Log requested sample filtering and displayed filtered/total cell counts (#101).
* Use all samples when sample selection is blank (#102).
* Report found Heatmap features before missing requested features without duplicate warnings (#103).
* Identify unavailable requested samples and list available sample names.
