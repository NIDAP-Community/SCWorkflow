#' Create Volcano Plots for Differential Expression Analysis
#'
#' Creates enhanced volcano plots for differential expression analysis using the
#' EnhancedVolcano package. Generates both static and interactive (plotly) versions
#' for each comparison specified in the input data.
#'
#' @param DEGAnalysis Data frame containing differential expression analysis results.
#'   Must include columns for gene/feature IDs, log fold changes, and p-values.
#' @param label.col Character string specifying the column name containing feature
#'   IDs (e.g., gene names). Default is "Gene".
#' @param sig.col Character vector of column names containing significance values
#'   (p-values or adjusted p-values). Can specify multiple comparisons.
#' @param lfc.col Character vector of column names containing log2 fold change values.
#'   Must have the same length as \code{sig.col}.
#' @param pCutoff Numeric value for the significance threshold. Default is 0.001.
#' @param FCcutoff Numeric value for the log2 fold change threshold. Default is 1.0.
#' @param value_to_sort_the_output_dataset Character string specifying how to sort
#'   features for labeling. Options are "p-value" or "fold-change". Default is "p-value".
#' @param no_genes_to_label Integer specifying the number of top features to label
#'   on the plot. Default is 30.
#' @param use_only_addition_labels Logical indicating whether to label only the
#'   features specified in \code{additional_labels}. Default is FALSE.
#' @param additional_labels Character string of additional feature names to label,
#'   separated by commas or spaces. Default is "".
#' @param is_red Logical indicating whether to restrict labeled features to only
#'   those passing both significance and fold change thresholds. Default is TRUE.
#' @param labSize Numeric value for the size of feature labels. Default is 4.
#' @param change_sig_name Character string for custom y-axis label. Default is "p-value".
#' @param change_lfc_name Character string for custom x-axis label. Default is "log2FC".
#' @param title Character string for the plot title. Default is "Volcano Plots".
#' @param use_custom_lab Logical indicating whether to use custom axis labels.
#'   Default is FALSE.
#' @param ylim Numeric value for the maximum y-axis limit. Set to 0 for automatic
#'   scaling. Default is 0.
#' @param custom_xlim Character string for custom x-axis limits. Can be empty for
#'   automatic scaling, a single number for symmetrical limits, or two comma-separated
#'   numbers for asymmetrical limits. Default is "".
#' @param xlim_additional Numeric value to add padding to x-axis limits. Default is 0.
#' @param ylim_additional Numeric value to add padding to y-axis limits. Default is 0.
#' @param axisLabSize Numeric value for the size of axis labels. Default is 24.
#' @param pointSize Numeric value for the size of points on the plot. Default is 2.
#'
#' @return A named list containing:
#'   \item{data}{Original input data frame with added rank columns for each comparison}
#'   \item{comparison_static}{Static ggplot2 volcano plot for each comparison}
#'   \item{comparison_interactive}{Interactive plotly volcano plot for each comparison}
#'
#' @examples
#' \dontrun{
#' # Create sample DEG data
#' set.seed(123)
#' deg_data <- data.frame(
#'   Gene = paste0("Gene", 1:100),
#'   `WT-KO_logFC` = rnorm(100, 0, 2),
#'   `WT-KO_pval` = runif(100, 0, 0.1),
#'   check.names = FALSE
#' )
#'
#' # Generate volcano plots
#' results <- Volcano_Plot(
#'   DEGAnalysis = deg_data,
#'   label.col = "Gene",
#'   sig.col = "WT-KO_pval",
#'   lfc.col = "WT-KO_logFC",
#'   pCutoff = 0.05,
#'   FCcutoff = 1.0
#' )
#'
#' # Access the plots
#' print(results$`WT-KO_static`)
#' print(results$`WT-KO_interactive`)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_split
#' @importFrom ggrepel geom_text_repel
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom plotly ggplotly
#' @importFrom grid grid.newpage
#' @export
Volcano_Plot <- function(DEGAnalysis,
                         label.col = "Gene",
                         sig.col,
                         lfc.col,
                         pCutoff = 0.001,
                         FCcutoff = 1.0,
                         value_to_sort_the_output_dataset = "p-value",
                         no_genes_to_label = 30,
                         use_only_addition_labels = FALSE,
                         additional_labels = "",
                         is_red = TRUE,
                         labSize = 4,
                         change_sig_name = "p-value",
                         change_lfc_name = "log2FC",
                         title = "Volcano Plots",
                         use_custom_lab = FALSE,
                         ylim = 0,
                         custom_xlim = "",
                         xlim_additional = 0,
                         ylim_additional = 0,
                         axisLabSize = 24,
                         pointSize = 2) {

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  df.orig <- DEGAnalysis
  rank <- list()
  plot_list <- list()

  for(i in 1:length(lfc.col)){
    lfccol <- lfc.col[i]
    sigcol <- sig.col[i]
    columns_of_interest <- c(label.col, lfc.col[i], sig.col[i])
    df <- df.orig %>% 
      dplyr::select(one_of(columns_of_interest)) %>% 
      mutate(!!sym(lfccol) := replace_na(!!sym(lfccol), 0)) %>%
      mutate(!!sym(sigcol) := replace_na(!!sym(sigcol), 1))

    if (use_custom_lab == TRUE){
      if (nchar(change_lfc_name) == 0){lfc_name = lfc.col[i]}
      if (nchar(change_sig_name) == 0){sig_name = sig.col[i]}
      colnames(df) <- c(label.col, change_lfc_name, sig_name)
    } else {
      lfc_name = lfc.col[i]
      sig_name = sig.col[i]
    }

    group <- gsub("_pval|p_val_", "", sig_name)
    rank[[i]] <- -log10(df[[sig_name]]) * sign(df[[lfc_name]]) 
    names(rank)[i] <- paste0("C_", group, "_rank")

    cat(paste0("Genes in initial dataset: ", nrow(df), "\n"))

    # Select top genes by logFC or Significance
    if (value_to_sort_the_output_dataset == "fold-change") {
      df <- df %>% dplyr::arrange(desc(.data[[lfc_name]]))
    } else if (value_to_sort_the_output_dataset == "p-value") {
      df <- df %>% dplyr::arrange(.data[[sig_name]]) 
    }

    if (is_red) {
      df_sub <- df[df[[sigcol]] <= pCutoff & abs(df[[lfccol]]) >= FCcutoff, ]
    } else {
      df_sub <- df
    }

    genes_to_label <- as.character(df_sub[1:no_genes_to_label, label.col])

    # Modifying Additional Labels List:
    # Replace commas with spaces and split the string
    split_values <- unlist(strsplit(gsub(",", " ", additional_labels), " "))
    additional_labels_vec <- split_values[split_values != ""]

    filter <- additional_labels_vec %in% df[, label.col]
    missing_labels <- additional_labels_vec[!filter]
    additional_labels_vec <- additional_labels_vec[filter]

    if(length(missing_labels) > 0){
      cat("Could not find:\n")
      print(missing_labels)
    }

    if(use_only_addition_labels){
      genes_to_label <- additional_labels_vec
    } else {
      genes_to_label <- unique(append(genes_to_label, additional_labels_vec))
    }

    significant = vector(length = nrow(df))
    significant[] = "Not significant"
    significant[which(abs(df[, 2]) > FCcutoff)] = "Fold change only"
    significant[which(df[, 3] < pCutoff)] = "Significant only"
    significant[which(abs(df[, 2]) > FCcutoff & df[, 3] < pCutoff)] = "Significant and fold change"
    print(table(significant))

    # Fix pvalue == 0
    shapeCustom <- rep(19, nrow(df))
    maxy <- max(-log10(df[[sig_name]]), na.rm = TRUE)
    if(ylim > 0){
      maxy <- ylim
    }

    cat(paste0("Maxy: ", maxy, "\n"))
    if(maxy == Inf){
      # Sometimes, pvalues == 0
      keep <- df[[sig_name]] > 0
      df[[sig_name]][!keep] <- min(df[[sig_name]][keep])
      shapeCustom[!keep] <- 17

      maxy <- -log10(min(df[[sig_name]][keep]))
      cat("Some p-values equal zero. Adjusting y-limits.\n")
      cat(paste0("Maxy adjusted: ", maxy, "\n"))
    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[sig_name]]) <= maxy
    df[[sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom) <- rep("Exact", length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"

    # Remove if nothin' doin'
    if(all(shapeCustom == 19)){
      shapeCustom <- NULL
    }

    maxy <- ceiling(maxy)

    if (grepl("log", lfc.col[i])){
      xlab <- bquote(~Log[2]~ "fold change")
    } else {
      xlab <- "Fold change"
    }
    if (grepl("adj", sig.col[i])){
      ylab <- bquote(~-Log[10]~ "FDR")
    } else {
      ylab <- bquote(~-Log[10]~ "p-value")
    }
    if(use_custom_lab){
      if(lfc_name != lfc.col[i]){
        xlab <- gsub("_", " ", lfc_name)
      }
      if (sig_name != sig.col[i]){ 
        ylab <- gsub("_", " ", sig_name)
      }
    }

    # X-axis custom range change:
    if (custom_xlim == "") {
      xlim = c(floor(min(df[, lfc_name])) - xlim_additional, 
               ceiling(max(df[, lfc_name])) + xlim_additional)
    } else if (grepl(",", custom_xlim) == FALSE) {
      xlim = c(-1 * as.numeric(trimws(custom_xlim)), 
               as.numeric(trimws(custom_xlim)))
    } else {
      split_values <- strsplit(custom_xlim, ",")[[1]]
      x_min <- as.numeric(trimws(split_values[1]))
      x_max <- as.numeric(trimws(split_values[2]))
      xlim <- c(x_min, x_max)
    }

    # Create static plot
    p <- EnhancedVolcano(df, x = lfc_name, y = sig_name,
                         lab = df[, label.col],
                         selectLab = genes_to_label,
                         title = title,
                         subtitle = group,
                         xlab = xlab,
                         ylab = ylab,
                         xlim = xlim,
                         ylim = c(0, maxy + ylim_additional),
                         pCutoff = pCutoff,
                         FCcutoff = FCcutoff,
                         axisLabSize = axisLabSize,
                         labSize = labSize,
                         pointSize = pointSize,
                         shapeCustom = shapeCustom)

    # Create interactive plot with no labels
    p_empty <- EnhancedVolcano(df, x = lfc_name, y = sig_name,
                               lab = rep("", nrow(df)),
                               selectLab = NULL,
                               title = title,
                               subtitle = group,
                               xlab = xlab,
                               ylab = ylab,
                               xlim = xlim,
                               ylim = c(0, maxy + ylim_additional),
                               pCutoff = pCutoff,
                               FCcutoff = FCcutoff,
                               axisLabSize = axisLabSize,
                               labSize = labSize,
                               pointSize = pointSize,
                               shapeCustom = shapeCustom)

    # Extract the data used for plotting
    plot_data <- ggplot_build(p_empty)$data[[1]]

    pxx <- p_empty +
      xlab("Fold Change") +
      ylab("Significance") +
      theme_minimal() +
      geom_point(aes(
        text = paste("Gene:", df[[label.col]], 
                     "<br>Log2FC:", df[[lfc_name]], 
                     "<br>P-value:", df[[sig_name]]),
        colour = as.character(plot_data$colour),
        fill = as.character(plot_data$colour)
      ), 
      shape = 21,
      size = 2,
      stroke = 0.1) + 
      scale_fill_identity()

    # Add interactive hover labels for the gene names
    interactive_plot <- ggplotly(pxx, tooltip = c("text"))

    # Store plots in the list with descriptive names
    plot_list[[paste0(group, "_static")]] <- p
    plot_list[[paste0(group, "_interactive")]] <- interactive_plot
  }

  # Combine original data with rank columns
  df.final <- cbind(df.orig, do.call(cbind, rank))

  # Return list with data first, then all plots
  result <- c(list(data = df.final), plot_list)
  return(result)
}