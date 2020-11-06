# ----------------------------------------------------------------------------
# ExpressionData
# ----------------------------------------------------------------------------

expressionDataUpdateText <- function(){
  output$expressionDataInfo <- renderText({
    if(is.null(input$expressionDataTissueSelect) & !isTRUE(input$expressionDataCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      if(is.null(input$expressionDataCellLineSelect) & !isTRUE(input$expressionDataCheckCellLineAll)){
        invisible("INFO: Please select the cell line(s) in the right panel!")
      }else{
        if(is.null(input$expressionDataGeneSelect) & !isTRUE(input$expressionDataCheckGeneAll)){
          invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
        }else{
          invisible("INFO: Click Load data!")
        }
      }
    }
  })
}

#upon load display nothing
output$expressionDataTable <- renderDataTable({
})

#query database and create dataframe
expressionDataDataFrame <- reactive({
  #get selected species
  if(input$expressionDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$expressionDataSpeciesSelect
  }
  
  #get selected tissue
  if(isTRUE(input$expressionDataCheckTissueAll)){
    presel_tissue <- expressionDataTissueList()
  }else{
    presel_tissue <- local(input$expressionDataTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$expressionDataCheckCellLineAll)){
    presel_cell_line <- expressionDataCellLineList()
  }else{
    presel_cell_line <- local(input$expressionDataCellLineSelect)
  }
  
  if(isTRUE(input$expressionDataCheckGeneAll)){
    presel_genes <- expressionDataGeneList()
  }else{
    presel_genes <- local(input$expressionDataGeneSelect)
  }
  
  presel_genes<- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  sample_ids <- cellline_list_expressionData %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::filter(tissue_name %in%  presel_tissue) %>%
    dplyr::filter(cell_line_name %in% presel_cell_line) %>%
    dplyr::select(sample_id) %>%
    .$sample_id
  
  df <- con_expression %>%
    tbl("expression_data_values") %>%
    dplyr::select(sample_id, gene_symbol, entrez_id, expression_value) %>%
    dplyr::filter(entrez_id %in% presel_gene_entrez, gene_symbol %in% presel_gene_symbol, sample_id %in% sample_ids) %>%
    distinct() %>%
    collect()
  
  df <- df %>%
    left_join(cellline_list_expressionData %>% dplyr::select(sample_id, species, unit))
  
  if(input$expressionDataSpeciesSelect == "all"){
    df_human <- df %>%
      dplyr::filter(species == "human") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_mouse) %>% dplyr::filter(!is.na(Symbol_human)), by=c("gene_symbol" = "Symbol_human")) %>%
      dplyr::rename(Symbol_human = gene_symbol, EntrezID_human = entrez_id)
    
    df_mouse <- df %>%
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_human) %>% dplyr::filter(!is.na(Symbol_mouse)), by=c("gene_symbol" = "Symbol_mouse")) %>%
      dplyr::rename(Symbol_mouse = gene_symbol, EntrezID_mouse = entrez_id)
    
    df <- df_human %>% rbind(df_mouse)
  }
  
  df
  
})

#create datatable out of dataframe
expressionDataDataTable <- eventReactive(input$expressionDataLoadButton,{
  
  df <- expressionDataDataFrame() %>%
    dplyr::distinct()
  
  if (nrow(df) > 0) {
    
    brks <- seq(0, 20, length.out = 40)
    clrs <- round(seq(255, 5, length.out = (length(brks) + 1)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}

    nfreezeColumns <- 2
    
    if(input$expressionDataSpeciesSelect == "all"){
      
      dt <- df %>%
        dplyr::select(sample_id, expression_value, Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse) %>%
        spread(sample_id, expression_value) %>%
        arrange(Symbol_human, Symbol_mouse) %>%
        dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, everything())
      
      nfreezeColumns <- nfreezeColumns + 2
    }else{
      dt <- df %>%
        dplyr::select(sample_id, gene_symbol, entrez_id, expression_value) %>%
        spread(sample_id, expression_value) %>%
        arrange(gene_symbol)
    }
    
    #tooltips
    colnames_dt <- colnames(dt)
    # 
    # contrast_ids <- df$contrast_id %>% unique
    # tooltip <- ''
    # 
    # if(input$gwsBrowseScreenDatasetSelect %in% c("dropout")){
    #   for(i in 1:length(colnames_dt)){
    #     if(i < length(colnames_dt)){
    #       tooltip <- paste0(tooltip, "'", colnames_dt[i], "'",  ", " )
    #     }else{
    #       tooltip <- paste0(tooltip, "'", colnames_dt[i], "'")
    #     }
    #     if(colnames_dt[i] %in% contrast_ids){
    #       colnames_dt[i] <- contrasts %>% select(contrast_id, contrast_id_QC) %>% filter(contrast_id == colnames_dt[i]) %>% .$contrast_id_QC
    #     }
    #   }
    #   colnames(dt) <- colnames_dt
    # }
    
    # headerCallback <- c(
    #   "function(thead, data, start, end, display){",
    #   "  var $ths = $(thead).find('th');",
    #   "  $ths.css({'vertical-align': 'bottom', 'white-space': 'nowrap'});",
    #   "  var betterCells = [];",
    #   "  $ths.each(function(){",
    #   "    var cell = $(this);",
    #   "    var newDiv = $('<div>', {height: 'auto', width: 'auto'});",
    #   "    var newInnerDiv = $('<div>', {text: cell.text()});",
    #   "    newDiv.css({margin: 'auto'});",
    #   "    newInnerDiv.css({",
    #   "      transform: 'rotate(180deg)',",
    #   "      'writing-mode': 'tb-rl',",
    #   "      'white-space': 'nowrap'",
    #   "    });",
    #   "    newDiv.append(newInnerDiv);",
    #   "    betterCells.push(newDiv);",
    #   "  });",
    #   "  $ths.each(function(i){",
    #   "    $(this).html(betterCells[i]);",
    #   "  });",
    #   paste0("  var tooltips = [", tooltip, "];"),
    #   paste0("  for(var i=0; i<", length(colnames_dt), "; i++){"),
    #   "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
    #   "  }",
    #   "}")
    
    dt <- dt %>%
      DT::datatable(extensions = c('FixedColumns','FixedHeader'),
                    options = list(
                      autoWidth = FALSE,
                      headerCallback = JS(headerCallback),
                      scrollX=TRUE,
                      fixedColumns = list(leftColumns = nfreezeColumns),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
                      searchHighlight = TRUE
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
    formatStyle(seq(nfreezeColumns+1, length(colnames_dt),1),
                backgroundColor = styleInterval(brks, clrs))
    
    if(!is.null(input$expressionDataGeneSelect) | isTRUE(input$expressionDataCheckGeneAll)){
      output$expressionDataInfo <- renderText({
        "Info: Loading completed! Table shows log2-transformed TPM values (+1 pseudocount)"
      })
    }
    #display datatable
    dt
  }else{
    if(!is.null(input$expressionDataGeneSelect) | isTRUE(input$expressionDataCheckGeneAll)){
      output$expressionDataInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      expressionDataUpdateText()
    }
  }
})

#----------------------------------------------------------------------------
# Selectbox Lists
#----------------------------------------------------------------------------

expressionDataTissueList <- reactive({
  
  if(class(cellline_list_expressionData)[1] == "tbl_SQLiteConnection" & loadExpressionDataTissueList){
    cellline_list_expressionData <<- cellline_list_expressionData %>%
      collect()
    gene_list_expressionData <<- gene_list_expressionData %>%
      collect()
  }
  
  if(input$expressionDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$expressionDataSpeciesSelect
  }
  
  cellline_list_expressionData %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

expressionDataCellLineList <- reactive({
  if(input$expressionDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$expressionDataSpeciesSelect
  }
  
  preselTissue = cellline_list_expressionData %>%
    dplyr::filter(species %in%  speciesList) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    .$tissue_name
  
  if(!isTRUE(input$expressionDataCheckTissueAll) & !is.null(input$expressionDataTissueSelect)){
    preselTissue = cellline_list_expressionData %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::filter(tissue_name %in% input$expressionDataTissueSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
  }
  
  cellline_list_expressionData %>%
    dplyr::filter(species %in% speciesList, tissue_name %in% preselTissue) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
})

expressionDataGeneList <- reactive({
  if(!isTRUE(input$expressionDataCheckCellLineAll) & (is.null(input$expressionDataCellLineSelect))){
    NULL
  }else{
    if(input$expressionDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$expressionDataSpeciesSelect
    }
    
    gene_list_expressionData %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::select(gene_symbol, entrez_id) %>%
      collect() %>%
      dplyr::mutate(gene = ifelse(is.na(gene_symbol), paste0("No symbol found (", entrez_id, ")"), paste0(gene_symbol , " (", entrez_id, ")"))) %>%
      arrange(gene) %>%
      .$gene
  }
})


#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
observe(
  if(input$tabs == "expressionDataSidebar"){
    loadExpressionDataTissueList <<- T
    if(input$expressionDataSpeciesSelect == ""){
      select = "human"
    }else{
      select = input$expressionDataSpeciesSelect
    }
    updateRadioButtons(session, 'expressionDataSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = select, inline = T)
  }
)

observeEvent(input$expressionDataLoadButton, {
  output$expressionDataTable <- renderDataTable({
    expressionDataDataTable()
  })
})

observeEvent(input$expressionDataSpeciesSelect, {
  #select checkbox tissue
  updateCheckboxInput(session, 'expressionDataCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'expressionDataTissueSelect', choices = expressionDataTissueList(), server = TRUE)
  #update cell linne selectbox
  updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  #disable laod button
  disable("expressionDataLoadButton")
  expressionDataUpdateText()
})

observeEvent(input$expressionDataTissueSelect, {
  if(!is.null(input$expressionDataTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'expressionDataCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect))){
    enable("expressionDataCellLineSelect")
    enable("expressionDataCheckCellLineAll")
  }else{
    disable("expressionDataCellLineSelect")
    disable("expressionDataCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$expressionDataCheckTissueAll, {
  if(isTRUE(input$expressionDataCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'expressionDataTissueSelect', choices = expressionDataTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect))){
    enable("expressionDataCellLineSelect")
    enable("expressionDataCheckCellLineAll")
  }else{
    disable("expressionDataCellLineSelect")
    disable("expressionDataCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$expressionDataCellLineSelect, {
  if(!is.null(input$expressionDataCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$expressionDataCheckCellLineAll) | (!is.null(input$expressionDataCellLineSelect))) & (isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect)))) {
    enable("expressionDataGeneSelect")
    enable("expressionDataCheckGeneAll")
  }else{
    disable("expressionDataGeneSelect")
    disable("expressionDataCheckGeneAll")
  }
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$expressionDataCheckCellLineAll, {
  if(isTRUE(input$expressionDataCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$expressionDataCheckCellLineAll) | (!is.null(input$expressionDataCellLineSelect))) & (isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect)))) {
    enable("expressionDataGeneSelect")
    enable("expressionDataCheckGeneAll")
  }else{
    disable("expressionDataGeneSelect")
    disable("expressionDataCheckGeneAll")
  }
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$expressionDataGeneSelect, {
  #unselect library checkbox
  if(!is.null(input$expressionDataGeneSelect)){
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
  }
  if(isTRUE(input$expressionDataCheckGeneAll) | (!is.null(input$expressionDataGeneSelect))){
    enable("expressionDataLoadButton")
  }else{
    disable("expressionDataLoadButton")
  }
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$expressionDataCheckGeneAll, {
  #unselect library selectbox
  if(isTRUE(input$expressionDataCheckGeneAll)){
    geneList <- expressionDataGeneList()
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = geneList, server = TRUE)
  }
  if(isTRUE(input$expressionDataCheckGeneAll) | (!is.null(input$expressionDataGeneSelect))){
    enable("expressionDataLoadButton")
  }else{
    disable("expressionDataLoadButton")
  }
  expressionDataUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$expressionDataButtonDownload <- downloadHandler(
  filename = function() {
    if(isTRUE(input$expressionDataCheckGeneAll)){
      paste0(paste(c("expresssion_all_genes_", local(input$expressionDataTissueSelect), local(input$expressionDataCellLineSelect)),  collapse="_"), ".txt")
    }else{
      str_replace_all(string = paste0("expression_", paste(local(input$expressionDataGeneSelect),collapse="_"), ".txt"), pattern = " ", replacement = "_")
    }
  },
  content = function(file) {
    df <- expressionDataDataFrame()
    
    if (nrow(df) > 0) {
      dt <- df %>%
        dplyr::select(sample_id, gene_symbol, entrez_id, expression_value, unit) %>%
        spread(sample_id, expression_value) %>%
        arrange(gene_symbol) %>% write_tsv(file)
    }
  }
)