# ----------------------------------------------------------------------------
# groupTesting
# ----------------------------------------------------------------------------
groupTestingUpdateText <- function(){
  output$groupTestingInfo <- renderText({
    if(is.null(input$groupTestingTissueSelect) & !isTRUE(input$groupTestingCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      if(is.null(input$groupTestingCellLineSelect) & !isTRUE(input$groupTestingCheckCellLineAll)){
        invisible("INFO: Please select the cell line(s) in the right panel!")
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
  
  #get selected tissue
  if(isTRUE(input$groupTestingCheckTissueAll)){
    presel_tissue <- groupTestingTissueList()
  }else{
    presel_tissue <- local(input$groupTestingTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    presel_cell_line <- groupTestingCellLineList()
  }else{
    presel_cell_line <- local(input$groupTestingCellLineSelect)
  }
  
  #get selected cell lines rest
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    presel_cell_line_rest <- groupTestingRestCellLineList()
  }else{
    presel_cell_line_rest <- local(input$groupTestingRestCellLineSelect)
  }
  
  #get selected library
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    presel_library <- groupTestingLibraryList()
  }else{
    presel_library <- local(input$groupTestingLibrarySelect)
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    
    testing_contrasts_target <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line, library_id %in% presel_library, 
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="target") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_contrasts_rest <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line_rest, library_id %in% presel_library,
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="rest") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_contrasts_all <- testing_contrasts_target %>%
      bind_rows(testing_contrasts_rest)
    
    if(input$groupTestingSearchRadio == "gene"){
      dropout_data <- read_tsv(file="essentiality_data/human_scaledLFC_geneLevel_HQscreens.tsv", col_names = T) %>%
        mutate(guide_id = entrez_id)
    }else{
      dropout_data <- read_tsv(file="essentiality_data/human_scaledLFC_guideLevel_HQscreens.tsv", col_names = T)
    }
    symbols <- dropout_data %>%
      select(symbol, entrez_id) %>%
      distinct()
    
    dropout_data <- dropout_data %>%
      dplyr::select(symbol, entrez_id, guide_id, starts_with(as.character(testing_contrasts_all$contrast_id_QC))) %>%
      pivot_longer(matches(as.character(testing_contrasts_all$contrast_id_QC)), names_to = "contrast_id", values_to="effect_essentialome") %>%
      inner_join(testing_contrasts_all %>% select(contrast_id=contrast_id_QC, cellline_name, group)) %>%
      filter(!is.na(effect_essentialome))
    
    testing_data_target <- dropout_data %>%
      dplyr::filter(contrast_id %in% testing_contrasts_target$contrast_id_QC) %>%
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
    
    testing_data_rest <- dropout_data %>%
      dplyr::filter(contrast_id %in% testing_contrasts_rest$contrast_id_QC) %>%
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
      # filter(abs(avgLFC.x - avgLFC.y) > 0.05, sdLFC.x > 0 | sdLFC.y > 0) %>%
      mutate(avgLFC_Differential = avgLFC.x - avgLFC.y, medianLFC_Differential = medianLFC.x - medianLFC.y,
             # p_value = format(t.test(unlist(LFCs.x), unlist(LFCs.y), alternative = "two.sided", paired = FALSE, var.equal = FALSE)$p.value,3),
             log10_BF= log10(extractBF(ttestBF(x = unlist(LFCs.x), y=unlist(LFCs.y)))$bf)) %>%
      dplyr::ungroup() %>%
      arrange(desc(log10_BF)) %>%
      select(-LFCs.x, -LFCs.y)

      # testing_data_joined$p_adjusted <- p.adjust(testing_data_joined$p_value, "fdr")
      
      testing_data_joined
    
  }else{
    
    
  }
  testing_data_joined
})

#query database and create dataframe
groupTestingContrastsDataFrame <- reactive({
  
  #get selected tissue
  if(isTRUE(input$groupTestingCheckTissueAll)){
    presel_tissue <- groupTestingTissueList()
  }else{
    presel_tissue <- local(input$groupTestingTissueSelect)
  }
  
  #get selected cell lines
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    presel_cell_line <- groupTestingCellLineList()
  }else{
    presel_cell_line <- local(input$groupTestingCellLineSelect)
  }
  
  #get selected cell lines rest
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    presel_cell_line_rest <- groupTestingRestCellLineList()
  }else{
    presel_cell_line_rest <- local(input$groupTestingRestCellLineSelect)
  }
  
  #get selected library
  if(isTRUE(input$groupTestingCheckLibraryAll)){
    presel_library <- groupTestingLibraryList()
  }else{
    presel_library <- local(input$groupTestingLibrarySelect)
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    
    testing_contrasts_target <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line, library_id %in% presel_library, 
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="target") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_contrasts_rest <- contrasts %>% 
      filter(species == local(input$groupTestingSpeciesSelect), 
             type == "dropout", reference_type == "reference", 
             cellline_name %in% presel_cell_line_rest, library_id %in% presel_library,
             auc > 0.95, dynamic_range < -1.5) %>%
      mutate(group="rest") %>%
      rowwise() %>%
      mutate(contrast_id_QC = paste0(contrast_id, "_", abs(dynamic_range))) %>%
      ungroup()
    
    testing_contrasts_all <- testing_contrasts_target %>%
      bind_rows(testing_contrasts_rest)
  }
})

#create plot out of dataframe
groupTestingResultsPlot <- eventReactive(input$groupTestingLoadButton,{
  df <- groupTestingResultsDataFrame()
  if (nrow(df) > 0) {
    output$groupTestingInfo <- renderText({"INFO: Loading completed!"})
    
    df$text <- paste("Gene: ", df$symbol, "\n",
                     "BF: ", round(df$log10_BF,2), "\n",
                     "delta LFC: ", round(df$avgLFC_Differential,2))
    
    df$log10_BF<- 10^(-df$log10_BF)
    
    p<- EnhancedVolcano(df,
                        lab = df$symbol,
                        x = "avgLFC_Differential",
                        y = "log10_BF",
                        pCutoff = 0.1,
                        FCcutoff = 0.25,
                        pointSize = 1.0,
                        labSize = 2,
                        ylim = c(min(-log10(df$log10_BF)), max(-log10(df$log10_BF))),
                        xlim = c(min(df$avgLFC_Differential), max(df$avgLFC_Differential))) + 
      xlab("delta LFC") + 
      ylab("log10 Bayes Factor")
    
    p<-ggplotly(p)
    
    legend_names <- c('Not sig.', 'delta LFC', 'BF', 'BF & delta LFC')
    p$x$data[[1]]$hoverinfo <- "none"
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

        p <- p %>% 
      style(textposition = "top left", hovertext = df$text, hoverinfo = "text", traces = seq_along(p$x$data)) %>%
      layout(
        height = 800,
        title = list(
          font = list(size = 12)
        ),
        xaxis = list(
          title = list(
            font = list(size = 10)
          ),
          tickfont = list(size = 10)
        ),
        yaxis = list(
          title = list(
            font = list(size = 10)
          ),
          tickfont = list(size = 10)
        ),
        legend = list(
          title = list(text = "Legend"),
          font = list(size = 12),
          orientation = "v",
          x = 1.05,  # Adjust horizontal position if needed
          y = 0.5    # Adjust vertical position if needed
        )
      )
    p$x$data[[1]]$hoverinfo <- "none"
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
  
  if(input$groupTestingSearchRadio == "gene"){
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

############# groupt testing target group ##############

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
  
  cell_lines %>%
    dplyr::filter(tissue_name %in% preselTissue, species == local(input$groupTestingSpeciesSelect)) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
})

############# groupt testing rest group ##############

groupTestingRestTissueList <- reactive({
  
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
  
  preselTissueTarget <- groupTestingTissueList()
  if(!isTRUE(input$groupTestingCheckTissueAll) & !is.null(input$groupTestingTissueSelect)){
    preselTissueTarget <- input$groupTestingTissueSelect
  }
  
  tissues %>% 
    filter(!tissue_name %in% preselTissueTarget) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

groupTestingRestCellLineList <- reactive({
  
  preselLibrary <- groupTestingLibraryList()
  if(!isTRUE(input$groupTestingCheckLibraryAll) & !is.null(input$groupTestingLibrarySelect)){
    preselLibrary <- input$groupTestingLibrarySelect
  }
  
  preselTissue <- groupTestingRestTissueList()
  if(!isTRUE(input$groupTestingRestCheckTissueAll) & !is.null(input$groupTestingRestTissueSelect)){
    preselTissue <- input$groupTestingRestTissueSelect
  }
  
  if(input$groupTestingDatasetSelect == "dependency"){
    cell_lines <- contrasts %>%
      filter(library_id %in% preselLibrary) %>%
      dplyr::rename(cell_line_name = cellline_name)
  }else{
    cell_lines <- cellline_list_expressionData
  }
  
  cell_lines %>%
    dplyr::filter(tissue_name %in% preselTissue, species == local(input$groupTestingSpeciesSelect)) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
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
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
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
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
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
    updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  }
  
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
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
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
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
  updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
  #select checbox tissue
  updateCheckboxInput(session, 'groupTestingRestCheckTissueAll', value = FALSE)
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

############# groupt testing target group ##############

observeEvent(input$groupTestingTissueSelect, {
  if(!is.null(input$groupTestingTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'groupTestingCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
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
  }else{
    disable("groupTestingCellLineSelect")
    disable("groupTestingCheckCellLineAll")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCheckTissueAll, {
  if(isTRUE(input$groupTestingCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'groupTestingTissueSelect', choices = groupTestingTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
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
  }else{
    disable("groupTestingCellLineSelect")
    disable("groupTestingCheckCellLineAll")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCellLineSelect, {
  if(!is.null(input$groupTestingCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'groupTestingCheckCellLineAll', value = FALSE)
  }
  
  if(isTRUE(input$groupTestingCheckCellLineAll) | (!is.null(input$groupTestingCellLineSelect))){
    enable("groupTestingLoadButton")
  }else{
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingCheckCellLineAll, {
  if(isTRUE(input$groupTestingCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'groupTestingCellLineSelect', choices = groupTestingCellLineList(), server = TRUE)
  }
  
  if(isTRUE(input$groupTestingCheckCellLineAll) | (!is.null(input$groupTestingCellLineSelect))){
    enable("groupTestingLoadButton")
  }else{
    disable("groupTestingLoadButton")
  }
  
  groupTestingUpdateText()
  
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
  if(isTRUE(input$groupTestingRestCheckTissueAll) | (!is.null(input$groupTestingRestTissueSelect))){
    enable("groupTestingRestCellLineSelect")
    enable("groupTestingRestCheckCellLineAll")
  }else{
    disable("groupTestingRestCellLineSelect")
    disable("groupTestingRestCheckCellLineAll")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCheckTissueAll, {
  if(isTRUE(input$groupTestingRestCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'groupTestingRestTissueSelect', choices = groupTestingRestTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = TRUE)
  #update cell line selectb
  if(isTRUE(input$groupTestingRestCheckTissueAll) | (!is.null(input$groupTestingRestTissueSelect))){
    enable("groupTestingRestCellLineSelect")
    enable("groupTestingRestCheckCellLineAll")
  }else{
    disable("groupTestingRestCellLineSelect")
    disable("groupTestingRestCheckCellLineAll")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCellLineSelect, {
  if(!is.null(input$groupTestingRestCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'groupTestingRestCheckCellLineAll', value = FALSE)
  }
  
  if(isTRUE(input$groupTestingRestCheckCellLineAll) | (!is.null(input$groupTestingRestCellLineSelect))){
    enable("groupTestingRestLoadButton")
  }else{
    disable("groupTestingRestLoadButton")
  }
  
  groupTestingUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$groupTestingRestCheckCellLineAll, {
  if(isTRUE(input$groupTestingRestCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'groupTestingRestCellLineSelect', choices = groupTestingRestCellLineList(), server = TRUE)
  }
  
  if(isTRUE(input$groupTestingRestCheckCellLineAll) | (!is.null(input$groupTestingRestCellLineSelect))){
    enable("groupTestingRestLoadButton")
  }else{
    disable("groupTestingRestLoadButton")
  }
  
  groupTestingUpdateText()
  
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
    
    #plot
    df$log10_BF <- 10^(-df$log10_BF)
    p<-EnhancedVolcano(df,
                        lab = df$symbol,
                        x = "avgLFC_Differential",
                        y = "log10_BF",
                        pCutoff = 0.1,
                        FCcutoff = 0.25,
                        pointSize = 1.0,
                        labSize = 2,
                        max.overlaps = Inf,
                        maxoverlapsConnectors = Inf,
                        drawConnectors = TRUE,
                        widthConnectors = 0.25,
                        arrowheads = F,
                        legendLabels = c('Not sig.', 'delta LFC', 'BF', 'BF & delta LFC'),
                        legendPosition = 'right',
                        legendLabSize = 16,
                        legendIconSize = 5.0,
                        ylim = c(min(-log10(df$log10_BF)), max(-log10(df$log10_BF))),
                        xlim = c(min(df$avgLFC_Differential), max(df$avgLFC_Differential))) + 
      xlab("delta LFC") + 
      ylab("Bayes Factor (log10)")
    
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    files <- NULL;
    
    #write each sheet to a csv file, save the name
    fileName <- "group_testing_results.pdf"
    ggsave(fileName, p, "pdf", width=9, height=9)
    files <- c(fileName,files)
    fileName <- "group_testing_results.txt"
    write_tsv(results_table,fileName)
    files <- c(fileName,files)
    fileName <- "group_testing_contrasts.txt"
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