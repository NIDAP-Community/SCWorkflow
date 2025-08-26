# load all relevant packages
devtools::load_all()

# load harmonized object containing Harmony assay with backcalculated data
harm_object = readRDS('/data/bianjh/SCWorkflow/tests/testthat/fixtures/chariou_harmonized.rds')

# get metadata for downstream steps
meta <- object@meta.data

# annotate celltypes - doesn't require Harmony assay
source('/data/bianjh/SCWorkflow/R/Annotate_Cell_Types.R')

object <- annotateCellTypes(object,
         species = "Mouse",
         reduction.type = "umap",
         legend.dot.size = 2,
         do.finetuning = FALSE,
         local.celldex = NULL,
         use.clusters = NULL)

# dual labeling
source('/data/bianjh/SCWorkflow/R/Dual_Labeling.R')
object <- dualLabeling(object[[1]], 
            samples = unique(object[[1]]$orig.ident), 
            marker.1 = "Pdcd4", 
            marker.2 = "Postn", 
            marker.1.type = "SCT", 
            marker.2.type = "SCT", 
            data.reduction = "both", 
            point.size = 0.5, 
            point.shape = 16, 
            point.transparency = 0.5, 
            add.marker.thresholds = TRUE, 
            marker.1.threshold = 0.5, 
            marker.2.threshold = 0.5, 
            filter.data = TRUE, 
            marker.1.filter.direction = "greater than", 
            marker.2.filter.direction = "greater than", 
            apply.filter.1 = TRUE, 
            apply.filter.2 = TRUE, 
            filter.condition = TRUE, 
            parameter.name = "Harmony_Assay_Dual_Label_Test", 
            trim.marker.1 = FALSE, 
            trim.marker.2 = FALSE, 
            pre.scale.trim = 0.99, 
            display.unscaled.values = FALSE) 

# ModuleScore
source('/data/bianjh/SCWorkflow/R/ModuleScore.R')

marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
ms.threshold = paste(colnames(marker.table), rep(0, ncol(marker.table)))
general.class = colnames(marker.table)[1:3]

object = modScore(harm_object, marker.table, ms.threshold,
         general.class, lvl.vec = c(), reduction = "tsne", 
         nbins = 10, gradient.ft.size = 6, 
         violin.ft.size = 6, step.size = 0.1) 

# color by genes
source('/data/bianjh/SCWorkflow/R/Color_by_Gene.R')

samples.to.include = 'c("CD8dep",
                         "Combo",
                         "ENT",
                         "NHSIL12",
                         "PBS")'
gene = c("Lox",
         "Eomes",
         "Tspan13",
         "Cd33",
         "Dut",
         "Pde4b")

colorByGene(harm_object,
            samples.to.include,
            gene,
            reduction.type = "umap",
            number.of.rows = 0,
            return.seurat.object = FALSE,
            color = "red",
            point.size = 1,
            point.shape = 16,
            point.transparency = 0.5,
            use.cite.seq.data = FALSE,
            assay = "Harmony")

source('/data/bianjh/SCWorkflow/R/Color_by_Gene_Automatic.R')

samples.subset = unique(object$orig.ident)
samples.to.display = unique(object$orig.ident)
marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
cells.of.interest = colnames(marker.table)[4:6]

cbg_auto <- colorByMarkerTable(harm_object,
                   samples.subset,
                   samples.to.display,
                   marker.table,
                   cells.of.interest,
                   protein.presence = FALSE,
                   assay = "Harmony",
                   reduction.type = "umap",
                   point.transparency = 0.5,
                   point.shape = 16,
                   cite.seq = FALSE)

plot(cbg_auto)

source('/data/bianjh/SCWorkflow/R/Dotplot_by_Metadata.R')
metadata <- "orig.ident"
set.seed(15)
markers <- sample(rownames(object),10, replace = FALSE)
cells <- c("PBS","CD8dep","NHSIL12","ENT","Combo")

harm_dot <- dotPlotMet(harm_object,
           metadata,
           cells,
           markers,
           use_assay = "Harmony",
           plot.reverse = FALSE,
           cell.reverse.sort = FALSE,
           dot.color = "darkblue")

source('/data/bianjh/SCWorkflow/R/Violin_Plots_by_Metadata.R')

group = "orig_ident"
assay = 'Harmony'
slot = 'scale.data'
jitter_points = T
jitter_dot_size = 4
filter_outliers = F
outlier_low = 0.05
outlier_high = 0.95
facet_by = ""
set.seed(82)
genes = sample(rownames(object$Harmony@scale.data), 5,
               replace = FALSE)

violinPlot_mod(harm_object, 
              assay, 
              slot, 
              genes, 
              group, 
              facet_by = "", 
              filter_outliers = F,
              outlier_low = 0.05,
              outlier_high = 0.95,
              jitter_points, 
              jitter_dot_size) 

source('/data/bianjh/SCWorkflow/R/Heatmap.R')

sample.names <- c("PBS", "CD8dep", "ENT", "NHSIL12", "Combo")
metadata <- c("orig.ident")
set.seed(15)
add.gene.or.protein <- TRUE
transcripts <- sample(rownames(object), 10, replace = FALSE)
proteins <- NULL
rna.annotations <- transcripts[c(1, 2)]
protein.annotations <- NULL
plot.title <- "Heatmap_Chariou_test"

harm_heat <- heatmapSC(harm_object,
          sample.names,
          metadata,
          transcripts,
          use_assay = 'Harmony',
          proteins = NULL,
          heatmap.color = "Bu Yl Rd",
          plot.title = "Heatmap",
          add.gene.or.protein = FALSE,
          protein.annotations = NULL,
          rna.annotations = NULL,
          arrange.by.metadata = TRUE,
          add.row.names = TRUE,
          add.column.names = FALSE,
          row.font = 5,
          col.font = 5,
          legend.font = 5,
          row.height = 15,
          set.seed = 6,
          scale.data = TRUE,
          trim.outliers = TRUE,
          trim.outliers.percentage = 0.01,
          order.heatmap.rows = FALSE,
          row.order = c())

plot(harm_heat$plot)




