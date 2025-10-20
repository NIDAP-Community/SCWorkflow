#!/usr/bin/env bash
set -euo pipefail

# add scworkflow executable to the path
scworkflow=$(R -s -e "cat(system.file('exec','scworkflow', package='SCWorkflow'))")
export PATH="$PATH:$(dirname $scworkflow)"

## run functions from SCWorkflow
#scworkflow processRawData --json=json_args/processRawData.json
#scworkflow filterQC --json=json_args/filterQC.json
#scworkflow combineNormalize --json=json_args/combineNormalize.json
#scworkflow annotateCellTypes --json=json_args/annotateCellTypes.json
#scworkflow appendMetadataToSeuratObject --json=json_args/appendMetadataToSeuratObject.json
#scworkflow filterSeuratObjectByMetadata --json=json_args/filterSeuratObjectByMetadata.json
#scworkflow reclusterFilteredSeuratObject --json=json_args/reclusterFilteredSeuratObject.json
#scworkflow reclusterSeuratObject --json=json_args/reclusterSeuratObject.json
#scworkflow colorByGene --json=json_args/colorByGene.json
#scworkflow colorByMarkerTable --json=json_args/colorByMarkerTable.json
#scworkflow modScore --json=json_args/modScore.json
#scworkflow nameClusters --json=json_args/nameClusters.json
#scworkflow plotMetadata --json=json_args/plotMetadata.json
#scworkflow dotPlotMet --json=json_args/dotPlotMet.json
#scworkflow violinPlot_mod --json=json_args/violinPlot_mod.json
#scworkflow heatmapSC --json=json_args/heatmapSC.json
#scworkflow tSNE3D --json=json_args/tSNE3D.json
#scworkflow dualLabeling --json=json_args/dualLabeling.json
#scworkflow degGeneExpressionMarkers --json=json_args/degGeneExpressionMarkers.json

scworkflow harmonyBatchCorrect --json=json_args/harmonyBatchCorrect.json
