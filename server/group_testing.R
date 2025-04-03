# ----------------------------------------------------------------------------
# groupTesting
# ----------------------------------------------------------------------------
groupTestingUpdateText <- function(){
  output$groupTestingInfo <- renderText({
    if(is.null(input$groupTestingTissueSelect) & !isTRUE(input$groupTestingCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      if(is.null(input$groupTestingCellLineSelect) & !isTRUE(input$groupTestingCheckCellLineAll)){
        invisible("INFO: Please select the cell line(s) in the right panel! You can add additional filtering critera for gene dependency or expression.")
      }else{
        invisible("INFO: Click Load data!")
      }
    }
  })
}

#upon load display nothing
output$groupTestingPlot <- renderPlotly({
})
output$groupTestingDataTableResults <- renderDataTable({
})
output$groupTestingDataTableContrasts <- renderDataTable({
})

#query database and create dataframe
groupTestingResultsDataFrame <- reactive({
  
  presel_cell_line <- groupTestingCellLineList()
  
  presel_cell_line_rest <- groupTestingRestCellLineList()
  
  #get selected library
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    presel_library <- groupTestingLibraryList()
  }else{
    presel_library <- local(input$groupTestingLibrarySelect)
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    
    testing_data_target <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line, library_id %in% presel_library, 
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="target") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_data_rest <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line_rest, library_id %in% presel_library,
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="rest") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_data_all <- testing_data_target %>%
      bind_rows(testing_data_rest)
    
    if(input$groupTestingSearchRadio == "gene"){
      data <- read_tsv(file="essentiality_data/human_scaledLFC_geneLevel_HQscreens.tsv", col_names = T) %>%
        mutate(guide_id = entrez_id)
    }else{
      data <- read_tsv(file="essentiality_data/human_scaledLFC_guideLevel_HQscreens.tsv", col_names = T)
    }
    
    symbols <- data %>%
      select(symbol, entrez_id) %>%
      distinct()
    
    data <- data %>%
      dplyr::select(symbol, entrez_id, guide_id, starts_with(as.character(testing_data_all$contrast_id_QC))) %>%
      pivot_longer(matches(as.character(testing_data_all$contrast_id_QC)), names_to = "contrast_id", values_to="effect_essentialome") %>%
      inner_join(testing_data_all %>% select(contrast_id=contrast_id_QC, cellline_name, group)) %>%
      filter(!is.na(effect_essentialome))
    
    testing_data_target <- data %>%
      dplyr::filter(contrast_id %in% testing_data_target$contrast_id_QC) %>%
      dplyr::group_by(symbol, entrez_id, guide_id, cellline_name) %>%
      dplyr::summarize(effect_essentialome = median(effect_essentialome, na.rm=TRUE)) %>%
      dplyr::group_by(symbol, entrez_id, guide_id) %>%
      dplyr::filter(n() >= 2) %>%
      summarize(n=sum(!is.na(effect_essentialome)),
                LFCs=list(effect_essentialome),
                percentEssential=round(100*sum(effect_essentialome < -0.5, na.rm=TRUE)/n(), 2),
                avgLFC = round(mean(effect_essentialome, na.rm=TRUE), 2),
                medianLFC = round(median(effect_essentialome, na.rm=TRUE), 2),
                sdLFC = round(sd(effect_essentialome, na.rm=TRUE), 2)) %>%
      dplyr::filter(n >= 2) %>%
      dplyr::ungroup()
    
    testing_data_rest <- data %>%
      dplyr::filter(contrast_id %in% testing_data_rest$contrast_id_QC) %>%
      dplyr::group_by(symbol, entrez_id, guide_id) %>%
      dplyr::filter(n() >= 2) %>%
      summarize(n=sum(!is.na(effect_essentialome)),
                LFCs=list(effect_essentialome),
                percentEssential=100*sum(effect_essentialome < -0.5, na.rm=TRUE)/n(),
                avgLFC = mean(effect_essentialome, na.rm=TRUE),
                medianLFC = median(effect_essentialome, na.rm=TRUE),
                sdLFC = sd(effect_essentialome, na.rm=TRUE)) %>%
      dplyr::ungroup()
    
    testing_data_joined <- testing_data_target %>%
      inner_join(testing_data_rest, by=c(c("symbol", "entrez_id", "guide_id")))
    
    testing_data_joined <- testing_data_joined %>%
      dplyr::group_by(symbol, entrez_id, guide_id) %>%
      mutate(avgLFC_Differential = avgLFC.x - avgLFC.y, medianLFC_Differential = medianLFC.x - medianLFC.y,
             # p_value = format(t.test(unlist(LFCs.x), unlist(LFCs.y), alternative = "two.sided", paired = FALSE, var.equal = FALSE)$p.value,3),
             log10_BF= log10(extractBF(ttestBF(x = unlist(LFCs.x), y=unlist(LFCs.y)))$bf)) %>%
      dplyr::ungroup() %>%
      arrange(desc(log10_BF)) %>%
      select(-LFCs.x, -LFCs.y)
    
    # testing_data_joined$p_adjusted <- p.adjust(testing_data_joined$p_value, "fdr")
    testing_data_joined
    
    
  }else{
    con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm_tmm.db")
    
    celllines_expressionData <- con_expression %>%
      tbl("expression_data_meta_info") %>%
      distinct() %>%
      arrange(cell_line_name) %>%
      collect() %>%
      rename(cellline_name=cell_line_name)
    
    DBI::dbDisconnect(con_expression)
    
    testing_data_target <- celllines_expressionData %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             cellline_name %in% presel_cell_line) %>%
      mutate(group="target") %>%
      dplyr::rename(sample_name = sample_id)
    
    testing_data_rest <- celllines_expressionData %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             cellline_name %in% presel_cell_line_rest) %>%
      mutate(group="rest") %>%
      dplyr::rename(sample_name = sample_id)
    
    testing_data_all <- testing_data_target %>%
      bind_rows(testing_data_rest) 
    
    data <- read_tsv(file="expression_values_per_tissue/all_tissues_log2_TMM_rpkm_human_spread.tsv", col_names = T)
   
    symbols <- data %>%
      select(symbol, entrez_id) %>%
      distinct()
    
    data <- data %>%
      dplyr::select(symbol, entrez_id, starts_with(as.character(testing_data_all$sample_name))) %>%
      pivot_longer(starts_with(as.character(testing_data_all$sample_name)), names_to = "sample_name", values_to="log2_TMM_rpkm") %>%
      inner_join(testing_data_all %>% select(sample_name, cellline_name, group)) %>%
      filter(!is.na(log2_TMM_rpkm))
    
    testing_data_target <- data %>%
      dplyr::filter(sample_name %in% testing_data_target$sample_name) %>%
      dplyr::group_by(symbol, entrez_id, cellline_name) %>%
      dplyr::summarize(log2_TMM_rpkm = median(log2_TMM_rpkm, na.rm=TRUE)) %>%
      dplyr::group_by(symbol, entrez_id) %>%
      dplyr::filter(n() >= 2) %>%
      summarize(n=sum(!is.na(log2_TMM_rpkm)),
                log2_TMM_rpkms=list(log2_TMM_rpkm),
                avg_log2_TMM_rpkm = round(mean(log2_TMM_rpkm, na.rm=TRUE), 2),
                median_log2_TMM_rpkm = round(median(log2_TMM_rpkm, na.rm=TRUE), 2),
                sd_log2_TMM_rpkm = round(sd(log2_TMM_rpkm, na.rm=TRUE), 2)) %>%
      dplyr::filter(n >= 2) %>%
      dplyr::ungroup()
    
    testing_data_rest <- data %>%
      dplyr::filter(sample_name %in% testing_data_rest$sample_name) %>%
      dplyr::group_by(symbol, entrez_id) %>%
      dplyr::filter(n() >= 2) %>%
      summarize(n=sum(!is.na(log2_TMM_rpkm)),
                log2_TMM_rpkms=list(log2_TMM_rpkm),
                avg_log2_TMM_rpkm = mean(log2_TMM_rpkm, na.rm=TRUE),
                median_log2_TMM_rpkm = median(log2_TMM_rpkm, na.rm=TRUE),
                sd_log2_TMM_rpkm = sd(log2_TMM_rpkm, na.rm=TRUE)) %>%
      dplyr::ungroup()
    
    testing_data_joined <- testing_data_target %>%
      inner_join(testing_data_rest, by=c(c("symbol", "entrez_id")))
    
    testing_data_joined <- testing_data_joined %>%
      dplyr::group_by(symbol, entrez_id) %>%
      mutate(avg_log2_TMM_rpkm_Differential = avg_log2_TMM_rpkm.x - avg_log2_TMM_rpkm.y, median_log2_TMM_rpkm_Differential = median_log2_TMM_rpkm.x - median_log2_TMM_rpkm.y,
             # p_value = format(t.test(unlist(LFCs.x), unlist(LFCs.y), alternative = "two.sided", paired = FALSE, var.equal = FALSE)$p.value,3),
             log10_BF= ifelse(length(unique(unlist(log2_TMM_rpkms.x)))!=1 & length(unique(unlist(log2_TMM_rpkms.y)))!=1,
                              log10(extractBF(ttestBF(x = unlist(log2_TMM_rpkms.x), y=unlist(log2_TMM_rpkms.y)))$bf), NA)) %>%
      dplyr::ungroup() %>%
      arrange(desc(log10_BF)) %>%
      select(-log2_TMM_rpkms.x, -log2_TMM_rpkms.y)
    
    # testing_data_joined$p_adjusted <- p.adjust(testing_data_joined$p_value, "fdr")
    testing_data_joined 
  }
  
})

#query database and create dataframe
groupTestingContrastsDataFrame <- reactive({

  presel_cell_line <- groupTestingCellLineList()
  
  presel_cell_line_rest <- groupTestingRestCellLineList()
  
  #get selected library
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    presel_library <- groupTestingLibraryList()
  }else{
    presel_library <- local(input$groupTestingLibrarySelect)
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    
    testing_data_target <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line, library_id %in% presel_library, 
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="target") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_data_rest <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line_rest, library_id %in% presel_library,
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="rest") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_data_all <- testing_data_target %>%
      bind_rows(testing_data_rest)
  }else{
    con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm_tmm.db")
    
    celllines_expressionData <- con_expression %>%
      tbl("expression_data_meta_info") %>%
      distinct() %>%
      arrange(cell_line_name) %>%
      collect() %>%
      rename(cellline_name=cell_line_name)
    
    DBI::dbDisconnect(con_expression)
    
    testing_data_target <- celllines_expressionData %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             cellline_name %in% presel_cell_line) %>%
      mutate(group="target")
    
    testing_data_rest <- celllines_expressionData %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             cellline_name %in% presel_cell_line_rest) %>%
      mutate(group="rest")
    
    testing_data_all <- testing_data_target %>%
      bind_rows(testing_data_rest)
  }
  testing_data_all
})

#create plot out of dataframe
groupTestingResultsPlot <- eventReactive(input$groupTestingLoadButton,{
  df <- groupTestingResultsDataFrame()
  if (nrow(df) > 0) {
    output$groupTestingInfo <- renderText({"INFO: Loading completed!"})
    
    if(input$groupTestingDatasetSelect == "dependency"){
      values <- df$avgLFC_Differential
      df$values <- values
      delta <- "delta LFC"
      xaxis <- "avgLFC_Differential"
    }else{
      values <- df$avg_log2_TMM_rpkm_Differential
      df$values <- values
      delta <- "delta log2-TMM-RPKM"
      xaxis <- "avg_log2_TMM_rpkm"
    }
    
    df$text <- paste("Gene: ", df$symbol, "\n",
                     "BF: ", round(df$log10_BF,2), "\n",
                     paste0(delta,": ", round(values,2)))
    
    df$log10_BF<- 10^(-df$log10_BF)
    
    
    # Create a grouping variable using the same cutoffs as in EnhancedVolcano
    df$group <- with(df, ifelse((abs(values) >= 0.25) & (log10_BF < 0.1), paste0("BF & ", delta),
                                ifelse((abs(values) >= 0.25), delta,
                                       ifelse((log10_BF < 0.1), "BF", "Not sig."))))
    
    p<- EnhancedVolcano(df,
                        lab = df$symbol,
                        x = "values",
                        y = "log10_BF",
                        pCutoff = 0.1,
                        FCcutoff = 0.25,
                        pointSize = 1.0,
                        labSize = 2,
                        ylim = c(min(-log10(df$log10_BF)), max(-log10(df$log10_BF))),
                        xlim = c(min(values), max(values))) + 
      xlab(delta) + 
      ylab("log10 Bayes Factor")
    
    p<-ggplotly(p)
    
    legend_names <- c("Not sig.", delta, "BF", paste0("BF & ", delta))
    # Loop over the Plotly traces to update their legend names.
    trace_index <- 1
    for(i in seq_along(p$x$data)) {
      # Check if the trace has a non-empty name (these become the legend labels)
      if (!is.null(p$x$data[[i]]$name) && p$x$data[[i]]$name != "") {
        # Update the trace's name with our desired label
        p$x$data[[i]]$name <- legend_names[trace_index]
        trace_index <- trace_index + 1
        # If we've updated all desired legend names, exit the loop.
        if(trace_index > length(legend_names)) break
      }
    }
    
    # Now, update the hover text for each trace based on the group membership
    for (i in seq_along(p$x$data)) {
      trace_name <- p$x$data[[i]]$name
      if (trace_name %in% legend_names) {
        # Find the indices in df corresponding to this group
        idx <- which(df$group == trace_name)
        # Assign hover text only for these points
        p$x$data[[i]]$hovertext <- df$text[idx]
        p$x$data[[i]]$hoverinfo <- "text"
      }
    }

    p <- p %>% 
      style(textposition = "top left") %>%
      layout(
        height = 800,
        title = list(font = list(size = 12)),
        xaxis = list(title = list(font = list(size = 10)), tickfont = list(size = 10)),
        yaxis = list(title = list(font = list(size = 10)), tickfont = list(size = 10)),
        legend = list(
          title = list(text = "Legend"),
          font = list(size = 12),
          orientation = "v",
          x = 1.05,  # Adjust horizontal position if needed
          y = 0.5    # Adjust vertical position if needed
        )
      )
    
    # Assume p is your existing Plotly object.
    for(i in seq_along(p$x$data)) {
      # Check if the trace is a scatter trace (or has mode markers)
      if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "scatter") {
        p$x$data[[i]]$type <- "scattergl"
      }
      if (!is.null(p$x$data[[i]]$hoveron)) {
        p$x$data[[i]]$hoveron <- NULL
      }
    }
    p
    
  }else{
    NULL
  }
})

#create datatable out of dataframe
groupTestingResultsDataTable <- eventReactive(input$groupTestingLoadButton,{
  df <- groupTestingResultsDataFrame() %>%
    mutate_if(is.numeric, round, 2) %>%
    mutate(entrez_id = round(entrez_id))
  
  if(input$groupTestingSearchRadio == "gene" & input$groupTestingDatasetSelect == "dependency"){
    df <- df %>%
      select(-guide_id)
  }
  
  if (nrow(df) > 0) {
    output$groupTestingInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                class = "display nowrap",
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  # fixedColumns = list(leftColumns = 3),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    NULL
  }
})


#create datatable out of dataframe
groupTestingContrastsDataTable <- eventReactive(input$groupTestingLoadButton,{
  df <- groupTestingContrastsDataFrame()
  
  if (nrow(df) > 0) {
    output$groupTestingInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                class = "display nowrap",
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  # fixedColumns = list(leftColumns = 3),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    NULL
  }
})


#----------------------------------------------------------------------------
# Selectbox Lists
#----------------------------------------------------------------------------
groupTestingLibraryList <- reactive({
  
  contrasts %>%
    filter(species == local(input$groupTestingSpeciesSelect), type == "dropout", reference_type == "reference") %>%
    .$library_id %>%
    unique
  
})

groupTestingTissueList <- reactive({
  
  if(input$groupTestingDatasetSelect == "dependency"){
    contrasts_dropout <- contrasts %>%
      filter(species == local(input$groupTestingSpeciesSelect), type == "dropout", reference_type == "reference")
    
    if(isTRUE(input$groupTestingCheckLibraryAll)){
      library <- contrasts_dropout %>%
        .$library_id %>%
        unique
    }else{
      library <- local(input$groupTestingLibrarySelect)
    }
    tissues <- contrasts_dropout %>%
      filter(library_id %in% library)
  }else{
    tissues <- cellline_list_expressionData %>%
      filter(species == local(input$groupTestingSpeciesSelect))
  }
  
  tissues %>% 
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

groupTestingGeneDependencyList <- reactive({
  #get selected tissue
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    library <- contrasts %>%
      filter(species == local(input$groupTestingSpeciesSelect), type == "dropout", reference_type == "reference") %>%
      .$library_id %>%
      unique
  }else{
    library <- local(input$groupTestingLibrarySelect)
  }
  
  gene_list_screens %>%
    left_join(libraries %>% dplyr::select(library_id, species) %>% distinct) %>%
    dplyr::filter(library_id %in% local(library), species == local(input$groupTestingSpeciesSelect)) %>%
    dplyr::select(symbol) %>%
    distinct() %>%
    arrange(symbol) %>%
    .$symbol
  
})

groupTestingGeneExpressionList <- reactive({
  gene_list_expressionData %>%
    dplyr::filter(species == local(input$groupTestingSpeciesSelect)) %>%
    dplyr::select(symbol) %>%
    distinct() %>%
    arrange(symbol) %>%
    .$symbol
  
})

############# group testing target group ##############

groupTestingCellLineList <- reactive({
  
  preselLibrary <- groupTestingLibraryList()
  if(!isTRUE(input$groupTestingCheckLibraryAll) & !is.null(input$groupTestingLibrarySelect)){
    preselLibrary <- input$groupTestingLibrarySelect
  }
  
  preselTissue <- groupTestingTissueList()
  if(!isTRUE(input$groupTestingCheckTissueAll) & !is.null(input$groupTestingTissueSelect)){
    preselTissue <- input$groupTestingTissueSelect
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    cell_lines <- contrasts %>%
      filter(library_id %in% preselLibrary) %>%
      dplyr::rename(cell_line_name = cellline_name)
  }else{
    cell_lines <- cellline_list_expressionData
  }
  
  cell_lines <- cell_lines %>%
    dplyr::filter(tissue_name %in% preselTissue, species == local(input$groupTestingSpeciesSelect)) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  #get cell lines with requested gene dependency
  if(!is.null(input$groupTestingGeneDependencySelect) & length(cell_lines) > 0){
    con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
    
    gene_ids <- gene_list_screens %>%
      dplyr::filter(symbol %in% local(input$groupTestingGeneDependencySelect)) %>%
      select(gene_id) %>%
      distinct() %>%
      .$gene_id
    
    contrasts_filtered <- contrasts %>%
      dplyr::filter(type == "dropout", reference_type == "reference", dynamic_range <= -1.5, auc > 0.9, cellline_name %in% local(cell_lines)) %>%
      distinct()
    
    contrast_ids <- con %>%
      tbl("gene_stats") %>%
      dplyr::filter(contrast_id %in% local(contrasts_filtered$contrast_id), gene_id %in% local(gene_ids), effect_essentialome <= local(input$groupTestingGeneDependencySlider)) %>%
      dplyr::select(contrast_id, effect_essentialome) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(contrast_id) %>%
      summarize(n_check=sum(effect_essentialome <= local(input$groupTestingGeneDependencySlider)), n_genes=length(local(input$groupTestingGeneDependencySelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$contrast_id %>%
      unique
    
    cellline_names <- contrasts_filtered %>%
      dplyr::filter(contrast_id %in% contrast_ids) %>%
      .$cellline_name %>%
      unique
    
    cell_lines <- cell_lines[cell_lines %in% cellline_names]
    DBI::dbDisconnect(con)
    
    if(length(cell_lines) == 0){
      
      showModal(modalDialog(
        title = "WARNING!", 
        paste0("WARNING: The specified gene dependency filtering criteria do not match any cell lines. Try removing/changing some criteria to find matching cell lines."),
        footer = tagList(
          modalButton("OK"),
        )
      ))
      
    }
  }
  
  #get cell lines with requested gene expression
  if(!is.null(input$groupTestingGeneExpressionSelect) & length(cell_lines) > 0){
    con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm_tmm.db")
    
    sample_ids_filtered <- cellline_list_expressionData %>%
      dplyr::filter(cell_line_name %in% local(cell_lines)) %>%
      dplyr::select(sample_id) %>%
      distinct() %>%
      collect() %>%
      .$sample_id
    
    sample_ids <- con_expression %>%
      tbl("expression_data_values") %>%
      dplyr::filter(sample_id %in% local(sample_ids_filtered), symbol %in% local(input$groupTestingGeneExpressionSelect), log2_tpm >= local(input$groupTestingGeneExpressionSlider)) %>%
      dplyr::select(sample_id, log2_tpm) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(sample_id) %>%
      summarize(n_check=sum(log2_tpm >= local(input$groupTestingGeneExpressionSlider)), n_genes=length(local(input$groupTestingGeneExpressionSelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$sample_id %>%
      unique
    
    cellline_names <- cellline_list_expressionData %>%
      dplyr::filter(sample_id %in% sample_ids) %>%
      .$cell_line_name %>%
      unique
    
    cell_lines <- cell_lines[cell_lines %in% cellline_names]
    DBI::dbDisconnect(con_expression)
    
    if(length(cell_lines) == 0){
      
      showModal(modalDialog(
        title = "WARNING!", 
        paste0("WARNING: The specified gene expression filtering criteria do not match any cell lines. Try removing/changing some criteria to find matching cell lines."),
        footer = tagList(
          modalButton("OK"),
        )
      ))
      
    }
  }
  
  cell_lines %>% unique
})

############# group testing rest group ##############

groupTestingRestCellLineList <- reactive({
  
  preselLibrary <- groupTestingLibraryList()
  if(!isTRUE(input$groupTestingCheckLibraryAll) & !is.null(input$groupTestingLibrarySelect)){
    preselLibrary <- input$groupTestingLibrarySelect
  }

  preselTissue <- groupTestingTissueList()
  if(!isTRUE(input$groupTestingRestCheckTissueAll) & !is.null(input$groupTestingRestTissueSelect)){
    preselTissue <- input$groupTestingRestTissueSelect
  }
  
  preselCellLinesTarget <- c()
  if(!isTRUE(input$groupTestingCheckCellLineAll) & !is.null(input$groupTestingCellLineSelect)){
    preselCellLinesTarget <- input$groupTestingCellLineSelect
  }else{
    if(isTRUE(input$groupTestingCheckCellLineAll) & is.null(input$groupTestingCellLineSelect)){
      preselCellLinesTarget <- groupTestingCellLineList()
    }
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    cell_lines <- contrasts %>%
      filter(library_id %in% preselLibrary) %>%
      dplyr::rename(cell_line_name = cellline_name)
  }else{
    cell_lines <- cellline_list_expressionData
  }
  
  cell_lines <- cell_lines %>%
    dplyr::filter(tissue_name %in% preselTissue, species == local(input$groupTestingSpeciesSelect), !cell_line_name %in% preselCellLinesTarget) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name %>%
    unique
  
  #get cell lines with requested gene dependency
  if(!is.null(input$groupTestingRestGeneDependencySelect) & length(cell_lines) > 0){
    con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
    
    gene_ids <- gene_list_screens %>%
      dplyr::filter(symbol %in% local(input$groupTestingRestGeneDependencySelect)) %>%
      select(gene_id) %>%
      distinct() %>%
      .$gene_id
    
    contrasts_filtered <- contrasts %>%
      dplyr::filter(type == "dropout", reference_type == "reference", dynamic_range <= -1.5, auc > 0.9, cellline_name %in% local(cell_lines)) %>%
      distinct()
    
    contrast_ids <- con %>%
      tbl("gene_stats") %>%
      dplyr::filter(contrast_id %in% local(contrasts_filtered$contrast_id), gene_id %in% local(gene_ids), effect_essentialome >= local(input$groupTestingRestGeneDependencySlider)) %>%
      dplyr::select(contrast_id, effect_essentialome) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(contrast_id) %>%
      summarize(n_check=sum(effect_essentialome >= local(input$groupTestingRestGeneDependencySlider)), n_genes=length(local(input$groupTestingRestGeneDependencySelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$contrast_id %>%
      unique
    
    cellline_names <- contrasts_filtered %>%
      dplyr::filter(contrast_id %in% contrast_ids) %>%
      .$cellline_name %>%
      unique
    
    cell_lines <- cell_lines[cell_lines %in% cellline_names]
    DBI::dbDisconnect(con)
    
    if(length(cell_lines) == 0){
      
      showModal(modalDialog(
        title = "WARNING!", 
        paste0("WARNING: The specified gene dependency filtering criteria do not match any cell lines. Try removing/changing some criteria to find matching cell lines."),
        footer = tagList(
          modalButton("OK"),
        )
      ))
      
    }
  }
  
  #get cell lines with requested gene expression
  if(!is.null(input$groupTestingGeneExpressionSelect) & length(cell_lines) > 0){
    con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm_tmm.db")
    
    sample_ids_filtered <- cellline_list_expressionData %>%
      dplyr::filter(cell_line_name %in% local(cell_lines)) %>%
      dplyr::select(sample_id) %>%
      distinct() %>%
      collect() %>%
      .$sample_id
    
    sample_ids <- con_expression %>%
      tbl("expression_data_values") %>%
      dplyr::filter(sample_id %in% local(sample_ids_filtered), symbol %in% local(input$groupTestingGeneExpressionSelect), log2_tpm <= local(input$groupTestingGeneExpressionSlider)) %>%
      dplyr::select(sample_id, log2_tpm) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(sample_id) %>%
      summarize(n_check=sum(log2_tpm <= local(input$groupTestingGeneExpressionSlider)), n_genes=length(local(input$groupTestingGeneExpressionSelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$sample_id %>%
      unique
    
    cellline_names <- cellline_list_expressionData %>%
      dplyr::filter(sample_id %in% sample_ids) %>%
      .$cell_line_name %>%
      unique
    
    cell_lines <- cell_lines[cell_lines %in% cellline_names]
    DBI::dbDisconnect(con_expression)
    
    if(length(cell_lines) == 0){
      
      showModal(modalDialog(
        title = "WARNING!", 
        paste0("WARNING: The specified gene expression filtering criteria do not match any cell lines. Try removing/changing some criteria to find matching cell lines."),
        footer = tagList(
          modalButton("OK"),
        )
      ))
      
    }
    
  }
  
  cell_lines %>% 
    unique
})


#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
observeEvent(input$groupTestingLoadButton, {
  output$groupTestingPlot <- renderPlotly({
    p <- groupTestingResultsPlot()
    if(!is.null(p)){
      p
    }
  })
  output$groupTestingDataTableResults <- renderDataTable({
    dt <- groupTestingResultsDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  output$groupTestingDataTableContrasts <- renderDataTable({
    dt <- groupTestingContrastsDataTable()
    if(!is.null(dt)){
      dt
    }
  })
})

observeEvent(input$groupTestingSpeciesSelect, {
  #update library selectbox
  updateSelectizeInput(session, 'groupTestingLibrarySelect', choices = groupTestingLibraryList(), server = TRUE)
  #select checbox library
  updateCheckboxInput(session, 'groupTestingCheckLibraryAll', value = TRUE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update gene dependency selectbox
  updateSelectizeInput(session, 'groupTestingGeneDependencySelect', choices = groupTestingGeneDependencyList(), server = TRUE)
  #update gene expression selectbox
  updateSelectizeInput(session, 'groupTestingGeneExpressionSelect', choices = groupTestingGeneExpressionList(), server = TRUE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  #update gene dependency selectbox
  updateSelectizeInput(session, 'groupTestingRestGeneDependencySelect', choices = groupTestingGeneDependencyList(), server = TRUE)
  #update gene expression selectbox
  updateSelectizeInput(session, 'groupTestingRestGeneExpressionSelect', choices = groupTestingGeneExpressionList(), server = TRUE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingDatasetSelect, {
  
  if(input$groupTestingDatasetSelect == "dependency"){
    enable("groupTestingSearchRadio")
    enable("groupTestingLibrarySelect")
    enable("groupTestingCheckLibraryAll")
  }else{
    disable("groupTestingSearchRadio")
    disable("groupTestingLibrarySelect")
    disable("groupTestingCheckLibraryAll")
  }
  
  #update library selectbox
  updateSelectizeInput(session, 'groupTestingLibrarySelect', choices = groupTestingLibraryList(), server = TRUE)
  #select checbox library
  updateCheckboxInput(session, 'groupTestingCheckLibraryAll', value = TRUE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingSearchRadio, {
  if(!is.null(input$groupTestingTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = TRUE)
  }
  
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingLibrarySelect, {
  if(!is.null(input$groupTestingLibrarySelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'groupTestingCheckLibraryAll', value = FALSE)
  }else{
    updateCheckboxInput(session, 'groupTestingCheckLibraryAll', value = TRUE)
  }
  
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select tissue checkbox
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCheckLibraryAll, {
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'groupTestingLibrarySelect', choices = groupTestingLibraryList(), server = TRUE)
  }
  
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select tissue checkbox
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

############# group testing target group ##############

observeEvent(input$groupTestingTissueSelect, {
  if(!is.null(input$groupTestingTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = TRUE)
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$groupTestingCheckTissueAll) | (!is.null(input$groupTestingTissueSelect))){
    enable("groupTestingCellLineSelect")
    enable("groupTestingCheckCellLineAll")
    enable("groupTestingGeneDependencySelect")
    enable("groupTestingGeneExpressionSelect")
  }else{
    disable("groupTestingCellLineSelect")
    disable("groupTestingCheckCellLineAll")
    disable("groupTestingGeneDependencySelect")
    disable("groupTestingGeneExpressionSelect")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCheckTissueAll, {
  if(isTRUE(input$groupTestingCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = TRUE)
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$groupTestingCheckTissueAll) | (!is.null(input$groupTestingTissueSelect))){
    enable("groupTestingCellLineSelect")
    enable("groupTestingCheckCellLineAll")
    enable("groupTestingGeneDependencySelect")
    enable("groupTestingGeneExpressionSelect")
  }else{
    disable("groupTestingCellLineSelect")
    disable("groupTestingCheckCellLineAll")
    disable("groupTestingGeneDependencySelect")
    disable("groupTestingGeneExpressionSelect")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCellLineSelect, {
  if(!is.null(input$groupTestingCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = TRUE)
  
  if(isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)){
    enable("groupTestingRest")
    if(isTRUE(input$groupTestingRestCheckCellLineAll) | !is.null(input$groupTestingRestCellLineSelect)){
      enable("groupTestingLoadButton")
    }
  }else{
    disable("groupTestingRest")
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCheckCellLineAll, {
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = TRUE)
  
  if(isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)){
    enable("groupTestingRest")
    if(length(groupTestingCellLineList()) > 0 & length(groupTestingTissueList()) > 0 & (isTRUE(input$groupTestingRestCheckCellLineAll) | !is.null(input$groupTestingRestCellLineSelect))){
      enable("groupTestingLoadButton")
    }
  }else{
    disable("groupTestingRest")
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingGeneDependencySelect, {
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingGeneExpressionSelect, {
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingGeneDependencySlider, {
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingGeneExpressionSlider, {
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
}, ignoreNULL = FALSE)

############# groupt testing rest group ##############
observeEvent(input$groupTestingRestTissueSelect, {
  if(!is.null(input$groupTestingRestTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)){
    if(isTRUE(input$groupTestingRestCheckTissueAll) | (!is.null(input$groupTestingRestTissueSelect))){
      enable("groupTestingRestCellLineSelect")
      enable("groupTestingRestCheckCellLineAll")
      enable("groupTestingRestGeneDependencySelect")
      enable("groupTestingRestGeneExpressionSelect")
      enable("groupTestingRestGeneDependencySlider")
      enable("groupTestingRestGeneExpressionSlider")
    }else{
      disable("groupTestingRestCellLineSelect")
      disable("groupTestingRestCheckCellLineAll")
      disable("groupTestingRestGeneDependencySelect")
      disable("groupTestingRestGeneExpressionSelect")
      disable("groupTestingRestGeneDependencySlider")
      disable("groupTestingRestGeneExpressionSlider")
    }
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCheckTissueAll, {
  if(isTRUE(input$groupTestingRestCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  #update cell line selectb
  if(isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)){
    if(isTRUE(input$groupTestingRestCheckTissueAll) | (!is.null(input$groupTestingRestTissueSelect))){
      enable("groupTestingRestCellLineSelect")
      enable("groupTestingRestCheckCellLineAll")
      enable("groupTestingRestGeneDependencySelect")
      enable("groupTestingRestGeneExpressionSelect")
      enable("groupTestingRestGeneDependencySlider")
      enable("groupTestingRestGeneExpressionSlider")
      
    }else{
      disable("groupTestingRestCellLineSelect")
      disable("groupTestingRestCheckCellLineAll")
      disable("groupTestingRestGeneDependencySelect")
      disable("groupTestingRestGeneExpressionSelect")
      disable("groupTestingRestGeneDependencySlider")
      disable("groupTestingRestGeneExpressionSlider")
    }
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCellLineSelect, {
  if(!is.null(input$groupTestingRestCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  }
  
  if((isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)) & 
     (isTRUE(input$groupTestingRestCheckCellLineAll) | !is.null(input$groupTestingRestCellLineSelect))){
    enable("groupTestingLoadButton")
  }else{
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCheckCellLineAll, {
  if(isTRUE(input$groupTestingRestCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  }
  
  if(length(groupTestingCellLineList()) > 0 & length(groupTestingTissueList()) > 0 & 
     (isTRUE(input$groupTestingCheckCellLineAll) | !is.null(input$groupTestingCellLineSelect)) & 
     (isTRUE(input$groupTestingRestCheckCellLineAll) | !is.null(input$groupTestingRestCellLineSelect))){
    enable("groupTestingLoadButton")
  }else{
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestGeneDependencySelect, {
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestGeneExpressionSelect, {
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestGeneDependencySlider, {
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestGeneExpressionSlider, {
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
}, ignoreNULL = FALSE)


#############################################

observeEvent(input$groupTestingCancelModal, {
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  removeModal()
})

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------
output$groupTestingButtonDownload <- downloadHandler(
  filename = function() {
    paste0("group_testing_results.zip")
  },
  content = function(file) {
    results_table <- groupTestingResultsDataTable()$x$data
    contrast_table <- groupTestingContrastsDataTable()$x$data
    df <- groupTestingResultsDataFrame()
    
    if(input$groupTestingDatasetSelect == "dependency"){
      values <- df$avgLFC_Differential
      df$values <- values
      delta <- "delta LFC"
      xaxis <- "avgLFC_Differential"
    }else{
      values <- df$avg_log2_TMM_rpkm_Differential
      df$values <- values
      delta <- "delta log2-TMM-RPKM"
      xaxis <- "avg_log2_TMM_rpkm"
    }
    
    #plot
    df$log10_BF <- 10^(-df$log10_BF)
    p<-EnhancedVolcano(df,
                        lab = df$symbol,
                        x = "values",
                        y = "log10_BF",
                        pCutoff = 0.1,
                        FCcutoff = 0.25,
                        pointSize = 1.0,
                        labSize = 1,
                        max.overlaps = 5,
                        maxoverlapsConnectors = 5,
                        drawConnectors = TRUE,
                        widthConnectors = 0.25,
                        arrowheads = F,
                        legendLabels = c("Not sig.", delta, "BF", paste0("BF & ", delta)),
                        legendPosition = 'right',
                        legendLabSize = 10,
                        legendIconSize = 3.0,
                        ylim = c(min(-log10(df$log10_BF)), max(-log10(df$log10_BF))),
                        xlim = c(min(df$values), max(df$values))) + 
      xlab(delta) + 
      ylab("Bayes Factor (log10)")
    
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    files <- NULL;
    
    #write each sheet to a csv file, save the name
    fileName <- "group_testing_results.pdf"
    ggsave(fileName, p, "pdf", width=12, height=9)
    files <- c(fileName,files)
    fileName <- "group_testing_results.txt"
    write_tsv(results_table,fileName)
    files <- c(fileName,files)
    fileName <- "group_testing_metadata.txt"
    write_tsv(contrast_table,fileName)
    files <- c(fileName,files)
    
    if(!is.null(files)){
      #create the zip file
      zip(file,files, compression_level = 2)
    }else{
      NULL
    }
  }
)