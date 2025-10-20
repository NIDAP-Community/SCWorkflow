Pseudobulk_LimmaStats <- function(Pseudobulk_AggregateCounts,BRCATrainingSubset_SampleMetadata) {
  
  library(nidapFunctions)
  libs <- c("limma","ggplot2","edgeR","stringr","gridExtra","reshape2","grid","gtable","dplyr","tibble","tidyr","plotly")
  nidapLoadPackages(libs)
  
  counts_matrix <- Pseudobulk_AggregateCounts
  sample_metadata <- BRCATrainingSubset_SampleMetadata 
  samples_to_include = c("E1","E2","E3","H1","H2","H3","T1","T2","T3")
  gene_names_column <- "Gene"
  contrast_variable_columns <- c("Subtype")
  contrasts <- c("HER2-ER","HER2-TNBC","ER-TNBC")
  donor_variable_column = c("My_Sample_Names")
  sample_names_column = "My_Renamed_Samples"
  
  covariate_columns = c()
  summarization_method = "mean" 
  
  #Advanced:
  return_matrix = TRUE
  fold_change_threshold = 1.2
  first_pvalue_threshold = 0.05
  second_pvalue_threshold = 0.01
  
  #Visualization:
  row_number <- 30  # cut number of rows
  column_number <- 10  # cut number of columns per table
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(!all(make.names(colnames(counts_matrix))==colnames(counts_matrix))){
    print("Error: The following counts matrix column names are not valid:\n")
    print(colnames(counts_matrix)[make.names(colnames(counts_matrix))!=colnames(counts_matrix)])
    
    print("Likely causes are columns starting with numbers or other special characters eg spaces.")
    stop("Bad column names.")
  }
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  # Function to split data frame rows
  split_data_frame_rows <- function(data, chunk_size) {
    data <- as.data.frame(data)
    split(data, ceiling(seq_along(1:nrow(data)) / chunk_size))
    #split(data, ceiling( nrow(data)/ chunk_size))
  }
  
  # Function to split data frame columns
  split_data_frame_cols <- function(data, chunk_size) {
    col_indices <- split(seq_len(ncol(data)), ceiling(seq_along(seq_len(ncol(data))) / chunk_size))
    lapply(col_indices, function(cols) data[, cols, drop = FALSE])
  }
  
  # Function to print chunks with titles
  print_chunk <- function(chunk, index, title) {
    for(i in seq_along(colnames(chunk))) {
      if(class(chunk[[i]]) == "list") {
        chunk[[i]] <- unlist(lapply(chunk[[i]], function(x) paste(x, collapse = ",")))
      }
    }
    
    # Create a grob with title
    title_grob <- textGrob(title, gp = gpar(fontsize = 15, fontface = "bold"))
    table_grob <- tableGrob(chunk, theme = ttheme_default(base_size = 10))
    
    # Arrange the title and table in a vertical layout
    combined_grob <- arrangeGrob(title_grob, table_grob, ncol = 1, heights = c(1, 10))
    
    grid.newpage()
    grid.draw(combined_grob)
    return(NULL)
  }
  
  # Function to split and print a table
  split_and_print_table <- function(data, row_chunk_size, col_chunk_size, tabtitle) {
    # Split rows into chunks
    row_chunks <- suppressWarnings(split_data_frame_rows(data, row_chunk_size))
    row_index <- 1
    
    for (row_chunk in row_chunks) {
      row_start <- row_index * row_chunk_size - row_chunk_size + 1
      row_stop <- row_start + nrow(row_chunk) - 1
      
      # Split columns into chunks for each row chunk
      col_chunks <- suppressWarnings(split_data_frame_cols(row_chunk, col_chunk_size))
      col_index <- 1
      
      for (col_chunk in col_chunks) {
        col_start <- col_index * col_chunk_size - col_chunk_size + 1
        col_stop <- col_start + ncol(col_chunk) - 1
        
        title <- paste(tabtitle,"Subset Rows", row_start, "-", row_stop, "Cols", col_start, "-", col_stop)
        print_chunk(col_chunk, paste0("Row ", row_start, "-", row_stop, ", Col ", col_start, "-", col_stop), title)
        
        col_index <- col_index + 1
      }
      
      row_index <- row_index + 1
    }
  }
  
  plot_sankey <- function(sample_metadata, contrast_variable_columns, covariate_columns = NULL) {
    
    # Prepare the data for the Sankey plot  
    all_columns <- c(contrast_variable_columns, covariate_columns)
    
    if (length(all_columns) == 1) {
      message("Cannot draw a Sankey plot with only one contrast variable column and no covariate columns.")
      return(NULL)
    } else {
      # Transform the data into a long format
      long_data <- sample_metadata %>%
        pivot_longer(cols = all_of(all_columns), names_to = "Category", values_to = "Value")
      
      # Create nodes
      nodes <- data.frame(name = unique(long_data$Value))
      
      # Create links dynamically
      links <- NULL
      for (i in 1:(length(all_columns) - 1)) {
        temp_links <- sample_metadata %>%
          group_by(across(all_of(all_columns[c(i, i + 1)]))) %>%
          summarise(value = n(), .groups = 'drop') %>%
          mutate(source = match(.data[[all_columns[i]]], nodes$name) - 1,
                 target = match(.data[[all_columns[i + 1]]], nodes$name) - 1) %>% select(source, target, value)
        links <- bind_rows(links, temp_links)
      }
      
      # Create the plotly Sankey plot
      sankey <- plot_ly( type = "sankey",
                         orientation = "h",
                         node = list(
                           label = nodes$name,
                           pad = 15,
                           thickness = 20,
                           line = list(color = "black", width = 0.5)
                         ),
                         link = list(source = links$source,
                                     target = links$target,
                                     value = links$value
                         )
      ) 
      return(sankey)
    }
  }
  
  plot_sankey_with_confounders <- function(data, contrast_columns, covariate_column, confounder_names) {
    # Transform data into long format
    long_data <- data %>%
      pivot_longer(cols = c(contrast_columns, covariate_column), names_to = "Category", values_to = "Value")
    
    # Create nodes
    nodes <- data.frame(name = unique(long_data$Value))
    
    # Create links
    long_data <- long_data %>%
      group_by(.data[[sample_names_column]]) %>%
      arrange(Category) %>%
      mutate(Source = lag(Value), Target = Value) %>%
      filter(!is.na(Source)) %>%
      ungroup()
    
    links <- long_data %>%
      group_by(Source, Target) %>%
      summarise(count = n(), .groups = 'drop') %>%
      mutate(source = match(Source, nodes$name) - 1,
             target = match(Target, nodes$name) - 1) %>%
      select(source, target, count)
    
    # Highlight confounded nodes
    if(length(confounder_names) > 0) {
      nodes$group <- ifelse(nodes$name %in% confounder_names, "confounded", "non_confounded")
      nodes$color <- ifelse(nodes$name %in% confounder_names, "red", "lightblue")
    } else { # give different colors to each name
      nodes$group <- nodes$name
      nodes$color <- scales::hue_pal()(length(unique(nodes$name)))
    }
    
    # Highlight confounded links
    links$group <- ifelse(nodes$name[links$source + 1] %in% confounder_names | nodes$name[links$target + 1] %in% confounder_names, "confounded", "non_confounded")
    links$color <- ifelse(nodes$name[links$source + 1] %in% confounder_names | nodes$name[links$target + 1] %in% confounder_names, "red", "lightblue")
    
    # Check if there are any confounders
    if (length(confounder_names) == 0) {
      # If no confounders, use grey color for all links
      links$color <- "#A9A9A9"
    }
    
    # Create the plotly Sankey plot
    sankey <- plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes$name,
        color = nodes$color,
        pad = 15,
        thickness = 20,
        line = list(color = "black", width = 0.5)
      ),
      link = list(
        source = links$source,
        target = links$target,
        value = links$count,
        color = links$color
      )
    ) %>%
      layout(
        title = "Sankey Diagram",
        font = list(size = 10)
      )
    
    return(sankey)
  }
  
  # Function to check for confounders without using design matrix
  check_confounders <- function(data, contrast_vars, covariate_column) {
    confounders <- list()
    confounder_names <- list()  
    for (conf_var in contrast_vars) {
      for (level in levels(data[[covariate_column]])) {
        contingency_table <- table(data[[conf_var]], data[[covariate_column]] == level)
        for (row_name in rownames(contingency_table)) {
          if (contingency_table[row_name, 1] == 0 & contingency_table[row_name, 2] > 0) {
            detailed_message <- paste(row_name, "is completely confounded with", covariate_column, level)
            confounders[[conf_var]] <- c(confounders[[conf_var]], detailed_message)
            confounder_names[[conf_var]] <- c(confounder_names[[conf_var]], row_name)
          }
        }
      }
    }
    return(list(confounders = confounders, confounder_names = confounder_names))
  }
  
  # Function to draw an error message with text wrapping and respecting newlines
  draw_error_message <- function(message, color, width = 120) {
    # Split the message by newlines
    segments <- str_split(message, "\n")[[1]]
    # Wrap each segment
    wrapped_segments <- lapply(segments, str_wrap, width = width)
    # Combine the wrapped segments with newlines
    wrapped_message <- paste(unlist(wrapped_segments), collapse = "\n")
    
    grid.newpage()
    grid.rect(gp = gpar(fill = color, col = NA))
    grid.text(wrapped_message, x = 0.5, y = 0.5, 
              gp = gpar(fontface = "italic", cex = 3, col = "black"), 
              just = "center")
  }
  
  #Function to get gene lists based on fold change and p-value thresholds
  getgenelists <- function(FClimit,pvallimit,pval){
    upreggenes <- list()
    downreggenes <- list()  
    for(i in 1:length(contrasts)){
      if(pval == "pval"){
        results %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvals)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
        results %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvals)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]]        
      } else {
        results %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvaladj)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
        results %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvaladj)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]] 
      }
    }
    names(upreggenes) <- contrasts
    names(downreggenes) <- contrasts
    allreggenes <- rbind(unlist(upreggenes),unlist(downreggenes))
    rownames(allreggenes) <- c(paste0("upreg>",FClimit, ", ",pval,"<",pvallimit),paste0("downreg<-",FClimit, ", ",pval,"<",pvallimit))
    return(allreggenes)
  }
  
  # Function to calculate fold change from group means
  calculate_fold_change <- function(group_means, contrast) {
    # Split the contrast string into individual group names
    groups <- strsplit(contrast, "-")[[1]]
    
    # Extract the two group names for the fold change calculation
    group1 <- groups[1]
    group2 <- groups[2]
    
    # Calculate fold change
    fold_change <- group_means[[group1]] / group_means[[group2]]
    fold_change[is.nan(fold_change)] <- 1
    fold_change[fold_change < 1] <- -1/fold_change
    return(fold_change)
  }
  
  # Function to calculate the difference in fold changes
  calculate_fold_change_difference <- function(group_means, contrast) {
    # Split the contrast string into individual contrasts
    contrasts <- strsplit(contrast, "\\)-\\(")[[1]]
    
    # Remove the residual parentheses
    contrasts <- gsub("[()]", "", contrasts)
    
    # Calculate fold changes for each contrast
    fc1 <- calculate_fold_change(group_means, contrasts[1])
    fc2 <- calculate_fold_change(group_means, contrasts[2])
    delta <- fc1/fc2
    delta[is.nan(delta)] <- 1
    delta[delta < 1] <- -1/delta  
    return(delta)
  }
  
  # Custom function to calculate standard error
  se <- function(x) {
    n <- length(x[!is.na(x)])  # Only non-NA values should be counted
    sd(x, na.rm = TRUE) / sqrt(n)
  }
  
  wraplines <- function(y){
    j = unlist(strsplit(y,"-"))
    k = strwrap(j, width = 10)
    l = paste(k,collapse="\n-")
    return(l)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  #Make sure that sample metadata contains all the sample names
  sample_metadata.filt <- sample_metadata %>% dplyr::filter(.data[[sample_names_column]] %in% samples_to_include)
  rownames(sample_metadata.filt) <- sample_metadata.filt[[sample_names_column]]
  
  ##Make sure the counts matrix matches sample names selected:
  counts_matrix.filt  <- counts_matrix %>% select(all_of(c(gene_names_column,samples_to_include))) 
  func <- match.fun(summarization_method)
  counts_matrix.filt <- counts_matrix.filt %>% group_by(.data[[gene_names_column]]) %>% summarise_all(func) 
  counts_matrix.final <- counts_matrix.filt[,sample_metadata.filt[[sample_names_column]]] %>% as.data.frame()
  rownames(counts_matrix.final) <- counts_matrix.filt[[gene_names_column]]
  counts_matrix.final <- as.matrix(counts_matrix.final)
  
  # Check that all values in counts_matrix.final are numeric
  if (!all(sapply(counts_matrix.final, is.numeric))) {
    stop("Error: counts_matrix.final contains non-numeric values. Please ensure all values are numeric.")
  }
  
  #Ensure that there is only 1 sample per condition per subject, and print out the groups:
  if(length(donor_variable_column) > 0){
    sample_metadata.filt.dedup <- sample_metadata.filt %>%
      group_by(across(all_of(c(donor_variable_column, contrast_variable_columns)))) %>%
      slice(1) %>%
      ungroup()
    #Print out individual groups:
    grouped_summary <- sample_metadata.filt.dedup %>%
      group_by(across(all_of(c(donor_variable_column, contrast_variable_columns)))) %>%
      summarise(Sample_Names = toString(.data[[sample_names_column]]),
                .groups = 'drop')
    wide_table <- grouped_summary %>%
      pivot_wider(names_from = contrast_variable_columns, values_from = Sample_Names)
    wide_table <- replace(wide_table, is.na(wide_table), "")
  } else {
    sample_metadata.filt.dedup <- sample_metadata.filt
    wide_table <- NULL
  }
  
  #Create another contrast column called contmerge if 2 factor:
  if(length(contrast_variable_columns)>1){
    sample_metadata.final <- sample_metadata.filt.dedup %>% dplyr::mutate(contmerge = paste0(.data[[contrast_variable_columns[1]]],".",.data[[contrast_variable_columns[2]]])) %>%
      as.data.frame()
  } else {
    sample_metadata.final <- sample_metadata.filt.dedup %>% dplyr::mutate(contmerge = .data[[contrast_variable_columns]]) %>% as.data.frame()
  }
  
  #Make design matrix:
  if(length(covariate_columns) > 0) {
    dm.formula <- as.formula(paste("~0 + contmerge + ", paste(covariate_columns, sep="+", collapse="+")))
  } else {
    dm.formula <- as.formula("~0 + contmerge")
  }
  
  cat("\n\nstatistical formula:\n")
  print(dm.formula)
  cat(sprintf("\n\n Contmerge factor created from %s \n\n", paste(contrast_variable_columns, collapse = ", ")))
  
  #Set rownames of sample metadata before creating design matrix:
  rownames(sample_metadata.final) = sample_metadata.final[[sample_names_column]]
  # Convert covariate column to factor
  sample_metadata.final <- sample_metadata.final %>%
    mutate(across(all_of(covariate_columns), as.factor))
  
  design <- model.matrix(dm.formula,sample_metadata.final)
  
  for(i in 1:length(contrast_variable_columns)){
    colnames(design) <- gsub(paste0("^", contrast_variable_columns[i]), "", colnames(design))
  }
  colnames(design) <- gsub(":",".",colnames(design))
  colnames(design) <- gsub("contmerge", "", colnames(design))
  
  # Identify non-estimable factors
  nonest0 <- colnames(design)[colSums(design) == 0]
  nonest1 <- colnames(design)[colSums(design) == 1]
  
  # Combine non-estimable factors
  nonest <- c(nonest0, nonest1)
  
  # Print out non-estimable factors
  if (length(nonest) > 0) {
    
    # Print initial dimensions
    print("Variables before removing non-estimable factors:")
    print(colnames(design))
    
    # Remove non-estimable factors from the design matrix
    design <- design[, !colnames(design) %in% nonest]
    
    #Remove uncategorized samples from sample metadata and counts matrix
    counts_matrix.final <- (counts_matrix.final[,rownames(design)]) %>% as.data.frame
    sample_metadata.final <- sample_metadata.final %>% filter(.data[[sample_names_column]] %in% rownames(design))
    
    # Print dimensions after removing non-estimable factors
    print("Variables after removing non-estimable factors:")
    print(colnames(design))
    
    message("The following non-estimable factors have been removed from the design matrix as they do not provide sufficient variability:\n",
            "Factors with all zero counts: ", paste(nonest0, collapse = ", "), "\n",
            "Factors with only one non-zero count: ", paste(nonest1, collapse = ", "), "\n")
  } else {
    message("All factors in the design matrix provide sufficient variability. No non-estimable factors were found.")
  }
  
  # Perform QR decomposition
  qr_check <- qr(design)
  
  # Extract the rank from the QR object
  qr_rank <- qr_check$rank
  
  # Get the number of columns in the design matrix
  num_columns <- ncol(design)
  
  # Check for collinearity
  if (qr_rank < num_columns) {
    cat("Collinearity detected: The design matrix does not have full rank.\n")
    confounded_covariates <- colnames(design)[qr_check$pivot[(qr_rank + 1):num_columns]]
    print(confounded_covariates)
  } else {
    cat("No collinearity detected: The design matrix has full rank.\n")
    confounded_covariates <- NULL
  }
  
  if(!is.null(confounded_covariates)){   
    # Identify confounded variables - column where confounding variables are found
    confounding_variable_column <- covariate_columns[grepl(covariate_columns,confounded_covariates)][1]
    confounder_info <- check_confounders(sample_metadata.final,  contrast_variable_columns, confounding_variable_column)
    confounders <- confounder_info$confounders
    confounder_names <- unique(unlist(confounder_info$confounder_names))
    
    # Draw Sankey Plot to show the experimental design
    sankeyplot <- plot_sankey_with_confounders(sample_metadata.final, 
                                               c(contrast_variable_columns), 
                                               covariate_columns, confounder_names)
    message <- paste0("Error: The following covariates are confounded:\n")#, "\n", paste(confounders, collapse = ", "), "\n\n")
    for (covariate in names(confounders)) {
      if (length(confounders[[covariate]]) > 0) {
        message <- paste0(message, "\n", covariate, "\n",
                          paste(confounders[[covariate]], collapse = "\n "), "\n")
      } else {
        message <- paste0(message, covariate, " (no specific confounding variable identified)\n")}
    }
    message <- paste0(message, "\n\nThis indicates that your experimental design contains confounding variables, which can bias your results.\n",
                      "Consider revising your experimental design to ensure that each covariate provides unique information.\n",
                      "You may need to collect more data or adjust the grouping of your samples.")
    draw_error_message(message, "tomato3")
    print(sankeyplot)
  } else {
    sankeyplot <- plot_sankey(sample_metadata.final, 
                              c(contrast_variable_columns), 
                              covariate_columns)
    message <- "No collinearity was detected in the design matrix."
    draw_error_message(message, "lightgreen")
    if(!is.null(sankeyplot)){
      print(sankeyplot)    
    } 
    #Prepare plot for distribution of values across samples
    counts.df <- as.data.frame(counts_matrix.final) %>% rownames_to_column("Features")
    
    long_counts <- counts.df %>% pivot_longer(
      cols = -Features, 
      names_to = "SampleID", 
      values_to = "Count"
    ) 
    
    # DrawBoxplot
    p <- ggplot(long_counts, 
                aes(x = SampleID, y = Count, fill = factor(SampleID))) + 
      geom_boxplot() +
      scale_fill_manual(values = rainbow(length(unique(long_counts$SampleID)))) +
      ggtitle("Distribution of Values Across Samples") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "none") +
      labs(x = "", y = "Value")
    
    print(p)
    
    # Splitting and printing the pivoted table for donor samples
    if(!is.null(wide_table )){
      split_and_print_table(wide_table, row_number, column_number, "Pivot Table of Different Groups by Donor")
    }
    
    #Run duplicate correlation and lmfit:
    if(length(donor_variable_column) > 0){
      cat("Running linear mixed model\n\n")
      donor_variable_column <- donor_variable_column[1]
      corfit <- duplicateCorrelation(counts_matrix.final,design,block=sample_metadata.final[[donor_variable_column]]) 
      fit <- lmFit(counts_matrix.final, design, block=sample_metadata.final[[donor_variable_column]], correlation=corfit$consensus.correlation)
    } else {
      cat("Running linear model\n\n")
      fit <- lmFit(counts_matrix.final, design)
    } 
    
    cm <- makeContrasts(contrasts = contrasts, levels=design)
    fit <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit)
    
    #Calculate Group Means:
    contgroups <- colnames(design)
    group_means <- list() 
    group_se <- list()
    for (group_name in contgroups) {
      group_indices <- which(design[, group_name] == 1)
      group_means[[group_name]] <- rowMeans(counts_matrix.final[, group_indices, drop = FALSE])
      group_se[[group_name]] <- apply(counts_matrix.final[, group_indices, drop = FALSE], 1, se)
    }
    
    #Put together the results dataframe:
    tstat <- fit2$t
    pvals <- fit2$p.value
    logFC <- fit2$coefficients    
    FC <- 2^fit2$coefficients
    FC <- apply(FC, MARGIN = c(1,2), function(x) ifelse(x < 1, -1/x, x))
    SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
    
    #Format the final results table 
    group_means <- do.call(data.frame, group_means)
    colnames(group_means) <- paste(colnames(group_means),"Mean",sep = "_") 
    group_se <- do.call(data.frame, group_se)
    colnames(group_se) <- paste(colnames(group_se),"SE",sep = "_")
    
    colnames(logFC) <- paste(colnames(logFC),"logFC", sep="_")
    colnames(FC) <- paste(colnames(FC),"FC", sep="_")
    colnames(SE) = paste(colnames(SE), "SE", sep="_")
    colnames(tstat) <- paste(colnames(tstat),"tstat", sep="_")
    colnames(pvals) <- paste(colnames(pvals),"pval",sep="_")
    pvaladj <- apply(pvals,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladj) <- paste(colnames(fit$coefficients),"adjpval",sep="_")
    results <- cbind(group_means,group_se,FC,logFC,SE,tstat,pvals,pvaladj)
    results <- as.data.frame(results) %>% rownames_to_column("Gene") 
    
    #Format group sample numbers:
    sampsize <- colSums(design)
    titleval <- "Please note Sample size:"
    titletext <- paste(names(sampsize), sampsize, sep = "=", collapse = " \n ")
    sampletitle <- textGrob(paste(titleval,"\n",titletext,"\n\n\n"),gp=gpar(fontsize=10))
    contrast <- colnames(cm)
    connames <- strsplit(contrast,"-")
    connames <- lapply(connames,function(x) {gsub("\\(","",gsub("\\)","",x))})
    contrastsize <- lapply(connames,function(x) sampsize[unlist(x)])
    contrasttext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n") 
    contrasttext <- textGrob(paste("\n\n\nContrasts:\n",contrasttext),gp=gpar(fontsize=10))
    
    grid.newpage()
    grid.draw(sampletitle)
    grid.newpage()
    grid.draw(contrasttext)
    
    #Format genelist sizes for typical cutoffs:
    FCpval1 <- getgenelists(FClimit = fold_change_threshold, pvallimit = first_pvalue_threshold,"pval")
    FCpval2 <- getgenelists(FClimit = fold_change_threshold, pvallimit = second_pvalue_threshold,"pval")
    FCadjpval1 <- getgenelists(FClimit = fold_change_threshold, pvallimit = first_pvalue_threshold,"adjpval")
    FCadjpval2 <- getgenelists(FClimit = fold_change_threshold, pvallimit = second_pvalue_threshold,"adjpval")
    pvaltab <- rbind(FCpval1,FCpval2,FCadjpval1,FCadjpval2)
    colnames(pvaltab) <- sapply(colnames(pvaltab), function(x) wraplines(x))
    
    #Splitting and printing the table for genelist sizes:
    split_and_print_table(pvaltab, row_number, column_number, "Number of contrasts")
    
    colnames(results) <- gsub(" - ","-",colnames(results))
    colnames(results) <- gsub("\\(","",gsub("\\)","",colnames(results)))
    
    if(return_matrix == TRUE){
      results = as.data.frame(cbind(results, counts_matrix.final))
    }    
    
    return(results)
  }
}