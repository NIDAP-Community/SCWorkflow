# set uploaded file size limit ~ 2GB
options(shiny.maxRequestSize = 2*1024^3)

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  library(DT)
})

# Try to use SCWorkflow helpers; if package not installed, source locally
has_scworkflow <- requireNamespace("SCWorkflow", quietly = TRUE)
if (!has_scworkflow) {
  # From inst/shiny/ModuleScoreApp to package root is ../../../
  repo_root <- normalizePath(file.path(getwd(), "..", "..", ".."))
  helpers <- file.path(repo_root, "R", "ModuleScoreHelpers.R")
  if (file.exists(helpers)) {
    source(helpers)
  } else {
    stop(sprintf("Helpers not found at: %s", helpers))
  }
  if (!exists("compute_modscore_data", mode = "function")) {
    stop("compute_modscore_data not found after sourcing ModuleScoreHelpers.R")
  }
  if (!exists("build_modscore_plots", mode = "function")) {
    stop("build_modscore_plots not found after sourcing ModuleScoreHelpers.R")
  }
}

parse_marker_text <- function(txt) {
  # Expect lines like: CellType: gene1,gene2,gene3
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    ct <- trimws(parts[1])
    genes <- trimws(parts[2])
    gene_vec <- trimws(strsplit(genes, ",", fixed = TRUE)[[1]])
    gene_vec <- gene_vec[gene_vec != ""]
    out[[ct]] <- gene_vec
  }
  out
}

ui <- fluidPage(
  titlePanel("ModuleScore Explorer"),
  sidebarLayout(
    sidebarPanel(
      helpText("Upload a Seurat .rds object and define markers per celltype."),
      fileInput("obj_file", "Seurat object (.rds)", accept = ".rds"),
      textInput("new_celltype", "Celltype name", placeholder = "e.g., CD4_T"),
      selectizeInput(
        "new_genes", "Markers (autocomplete)",
        choices = NULL, multiple = TRUE,
        options = list(placeholder = "Type to search genes...", maxItems = 50)
      ),
      actionButton("add_marker_row", "Add to marker table"),
      actionButton("clear_marker_row", "Clear current entry"),
      hr(),
      fluidRow(
        column(12,
          h4("Marker table"),
          div(
            actionButton("edit_marker", "Edit selected"),
            actionButton("delete_marker", "Delete selected"),
            style = "margin-bottom: 10px;"
          ),
          DTOutput("marker_table")
        )
      ),
      selectInput("reduction", "Reduction", choices = c("tsne","umap","pca"), selected = "tsne"),
      uiOutput("meta_ui"),
      actionButton("compute", "Compute Module Scores"),
      hr(),
      uiOutput("celltype_ui"),
      
      numericInput("nbins", "nbins", value = 10, min = 5, max = 50),
      numericInput("gradient_size", "Gradient axis text size", value = 6, min = 4, max = 14),
      numericInput("violin_size", "Violin axis text size", value = 6, min = 4, max = 14),
      numericInput("step.size", "Axis step size", value = 0.1, min = 0.05, max = 0.5, step = 0.05)
    ),
    mainPanel(
      uiOutput("plots_ui"),
      hr(),
      h4("Computation log"),
      verbatimTextOutput("compute_log")   # new: log output
    )
  )
)

server <- function(input, output, session) {
  object_rv <- reactiveVal(NULL)
  ms_data_rv <- reactiveVal(NULL)
  celltypes_rv <- reactiveVal(character())
  compute_log_rv <- reactiveVal("")       # new: holds captured messages
  marker_rv <- reactiveVal(list())
  gene_choice_rv <- reactiveVal(character())
  thresholds_rv <- reactiveVal(setNames(numeric(0), character(0)))

  sanitize_id <- function(ct) {
    id <- gsub("[^A-Za-z0-9_]", "_", ct)
    if (id == "") id <- "ct"
    id
  }

  observeEvent(input$obj_file, {
    req(input$obj_file)
    obj <- readRDS(input$obj_file$datapath)
    object_rv(obj)
  })

  # new: populate metadata dropdown after object loads
  output$meta_ui <- renderUI({
    req(object_rv())
    meta_cols <- colnames(object_rv()@meta.data)
    if (length(meta_cols) == 0) return(NULL)
    # prefer common fields if present
    default <- if ("orig.ident" %in% meta_cols) "orig.ident" else if ("seurat_clusters" %in% meta_cols) "seurat_clusters" else meta_cols[1]
    selectInput("meta_var", "Metadata variable", choices = meta_cols, selected = default)
  })

  output$celltype_ui <- renderUI({
    cts <- celltypes_rv()
    if (length(cts) == 0) return(NULL)
    selectInput("celltype", "Celltype", choices = cts, selected = cts[[1]])
  })
  
  # dynamic plot sections with per-celltype threshold sliders
  output$plots_ui <- renderUI({
    res <- ms_data_rv(); cts <- celltypes_rv()
    if (is.null(res) || !length(cts)) return(NULL)
    th <- thresholds_rv()
    tagList(lapply(cts, function(ct) {
      id  <- sanitize_id(ct)
      val <- th[ct]
      if (is.na(val)) val <- 0  # safe default
      tagList(
        h4(ct),
        sliderInput(id, paste0("Threshold: ", ct), min = 0, max = 1, value = val, step = 0.01),
        fluidRow(column(6, plotOutput(paste0(id, "_g"))),
                 column(6, plotOutput(paste0(id, "_g1")))),
        fluidRow(column(6, plotOutput(paste0(id, "_g3")))),
        hr()
      )
    }))
  })

  # populate gene choices after object upload
  observeEvent(object_rv(), {
    req(object_rv()@assays$SCT)
    gene_choice_rv(rownames(object_rv()@assays$SCT@data))

    # new entry input
    updateSelectizeInput(session, "new_genes", choices = gene_choice_rv(),
     server = TRUE)
  })

  observeEvent(input$add_marker_row, {
    ct <- trimws(input$new_celltype)
    genes <- unique(trimws(input$new_genes))
    if (nzchar(ct) && length(genes) > 0) {
      current <- marker_rv()
      current[[ct]] <- genes
      marker_rv(current)

      # initialize threshold (if it doesn't yet exist)
      th <- thresholds_rv()
      if (is.null(th[ct])) th[ct] <- 0
      thresholds_rv(th)

      # Reset entry
      updateTextInput(session, "new_celltype", value = "")
      updateSelectizeInput(session, "new_genes", selected = character(0))
    } else {
      showNotification("Enter a celltype and at least one gene.", type = "warning")
    }
  })

  observeEvent(input$clear_marker_row, {
    updateTextInput(session, "new_celltype", value = "")
    updateSelectizeInput(session, "new_genes", selected = character(0), server = TRUE)
  })

  # render DT table from marker_rv
  marker_df <- reactive({
    mt <- marker_rv()
    if (length(mt) == 0) return(data.frame(celltype = character(), markers = character()))
    data.frame(
      celltype = names(mt),
      markers  = vapply(mt, function(x) paste(x, collapse = ", "), character(1)),
      stringsAsFactors = FALSE
    )
  })

  output$marker_table <- renderDT({
    DT::datatable(
      marker_df(),
      selection = "single",
      rownames = FALSE,
      options = list(pageLength = 5, lengthChange = FALSE)
    )
  })

  # Delete selected row
  observeEvent(input$delete_marker, {
    df <- marker_df()
    sel <- input$marker_table_rows_selected
    if (!length(sel)) {
      showNotification("Select a row to delete.", type = "warning")
      return(NULL)
    }
    ct <- df$celltype[sel]
    current <- marker_rv()
    current[[ct]] <- NULL
    marker_rv(current)

    th <- thresholds_rv()
    thresholds_rv(th[names(th) != ct])
  })

  # initialize value for storing celltype(s) to edit
  editing_ct <- reactiveVal(NULL)

  # Edit selected row via modal
  observeEvent(input$edit_marker, {
    df <- marker_df()
    sel <- input$marker_table_rows_selected
    if (!length(sel)) {
      showNotification("Select a row to edit.", type = "warning")
      return(NULL)
    }
    ct <- df$celltype[sel]
    current <- marker_rv()
    old_genes <- current[[ct]]

    editing_ct(ct)

    showModal(modalDialog(
      title = sprintf("Edit markers for %s", ct),
      textInput("edit_celltype", "Celltype name", value = ct),
      selectizeInput(
        "edit_genes", "Markers", choices = NULL,
        selected = old_genes, multiple = TRUE,
        options = list(placeholder = "Type to search genes...", maxItems = 50)
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_marker_edit", "Save")
      ),
      easyClose = TRUE
    ))

    choices <- isolate(gene_choice_rv())
    selected <- old_genes
  
  session$onFlushed(function(){
    # Guard: avoid reactive reads here; use captured values
    # Also tolerate cases where the modal was closed before this runs
    try(
      updateSelectizeInput(
        session, "edit_genes",
        choices = choices, selected = selected, server = TRUE
      ),
      silent = TRUE
    )
  }, once = TRUE)
})

  # Save edit back to marker_rv
  observeEvent(input$save_marker_edit, {
    req(input$edit_celltype)
    new_ct <- trimws(input$edit_celltype)
    new_genes <- unique(trimws(input$edit_genes))
    if (!nzchar(new_ct) || !length(new_genes)) {
      showNotification("Provide a celltype and at least one gene.", type = "error")
      return(NULL)
    }

    old_ct <- editing_ct()
    if(is.null(old_ct)){
      showNotification("No row selected for edit.", type = "error")
      return(NULL)
    }

    current <- marker_rv()

    # if renaming, remove old key
    if (!identical(old_ct, new_ct)) {
      current[[old_ct]] <- NULL
    }
    current[[new_ct]] <- new_genes
    marker_rv(current)

    # update thresholds mapping on rename
    th <- thresholds_rv()

    if (!identical(old_ct, new_ct)) {
      th <- th[names(th) != old_ct]
      th[new_ct] <- 0
    }

    thresholds_rv(th)

    removeModal()
  })

  observeEvent(input$compute, {
    req(object_rv())
    req(input$meta_var)
    marker.list <- marker_rv()
    if (length(marker.list) == 0) {
      showNotification("Marker table is empty. Add celltypes + genes first.", type = "error")
      return(NULL)
    }

    use.columns <- names(marker.list)
    celltypes_rv(use.columns)

    # reconcile thresholds per current celltypes
    old_th <- thresholds_rv()
    new_th <- setNames(numeric(length(use.columns)), use.columns)
    for (ct in use.columns) {
      new_th[ct] <- if (!is.null(old_th[ct])) old_th[ct] else 0
    }
    thresholds_rv(new_th)

    present_counts <- vapply(marker.list, function(genes)
      sum(genes %in% rownames(object_rv()@assays$SCT@data)), numeric(1))
    if (any(present_counts == 0)) {
      missing <- names(present_counts)[present_counts == 0]
      showNotification(sprintf("Skipping celltypes with all genes missing: %s",
                               paste(missing, collapse = ", ")), type = "warning")
      marker.list <- marker.list[present_counts > 0]
      use.columns <- names(marker.list)
      celltypes_rv(use.columns)
    }

    # capture printed output and messages from compute
    log_lines <- character(0)
    con <- textConnection("log_lines", "w")
    sink(con)
    sink(con, type = "message")
    on.exit({
      sink(type = "message")
      sink()
      close(con)
    }, add = TRUE)

    compute_fun <- if (has_scworkflow) getFromNamespace("compute_modscore_data", "SCWorkflow") else compute_modscore_data
    res <- compute_fun(object_rv(), marker.list, use.columns,
               reduction = input$reduction,
               nbins = input$nbins,
               group.var = input$meta_var)
    compute_log_rv(paste(log_lines, collapse = "\n"))

    if (inherits(res, "try-error")) {
      showNotification(paste("Compute failed:", attr(res, "condition")$message), type = "error")
      return(NULL)
    }

    ms_data_rv(res)
  })

  output$compute_log <- renderText({
    compute_log_rv()
  })

  # Render dynamic plots per celltype
  observe({
    res <- ms_data_rv()
    cts <- celltypes_rv()
    req(res, length(cts) > 0, input$meta_var)
    build_fun <- if (has_scworkflow) getFromNamespace("build_modscore_plots", "SCWorkflow") else build_modscore_plots
    for (ct in cts) {
      local({
        ct_local <- ct
        entry <- res[[ct_local]]
        req(entry)
        id <- sanitize_id(ct_local)
        output[[paste0(id, "_g")]] <- renderPlot({
          thr <- input[[id]]
          if (is.null(thr)) {
            th <- thresholds_rv()
            thr <- if (!is.null(th[ct_local])) th[ct_local] else 0
          }
          p <- build_fun(
            object = entry$object,
            m = entry$m,
            coords = entry$coords,
            clusid.df = entry$clusid.df,
            d = entry$density,
            threshold = thr,
            gradient.ft.size = input$gradient_size,
            violin.ft.size = input$violin_size,
            step.size = input$step.size,
            reduction = input$reduction,
            group.var = input$meta_var
          )
          p$g
        })
        output[[paste0(id, "_g1")]] <- renderPlot({
          thr <- input[[id]]
          if (is.null(thr)) {
            th <- thresholds_rv()
            thr <- if (!is.null(th[ct_local])) th[ct_local] else 0
          }
          p <- build_fun(
            object = entry$object,
            m = entry$m,
            coords = entry$coords,
            clusid.df = entry$clusid.df,
            d = entry$density,
            threshold = thr,
            gradient.ft.size = input$gradient_size,
            violin.ft.size = input$violin_size,
            step.size = input$step.size,
            reduction = input$reduction,
            group.var = input$meta_var
          )
          p$g1
        })
        output[[paste0(id, "_g3")]] <- renderPlot({
          thr <- input[[id]]
          if (is.null(thr)) {
            th <- thresholds_rv()
            thr <- if (!is.null(th[ct_local])) th[ct_local] else 0
          }
          p <- build_fun(
            object = entry$object,
            m = entry$m,
            coords = entry$coords,
            clusid.df = entry$clusid.df,
            d = entry$density,
            threshold = thr,
            gradient.ft.size = input$gradient_size,
            violin.ft.size = input$violin_size,
            step.size = input$step.size,
            reduction = input$reduction,
            group.var = input$meta_var
          )
          p$g3
        })
      })
    }
  })
}

shinyApp(ui, server)
