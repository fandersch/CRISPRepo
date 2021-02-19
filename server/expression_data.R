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
  
  sample_ids <- cellline_list_expressionData %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::filter(tissue_name %in%  presel_tissue) %>%
    dplyr::filter(cell_line_name %in% presel_cell_line) %>%
    dplyr::select(sample_id) %>%
    .$sample_id
  
  if(isTRUE(input$expressionDataCheckGeneAll)){
    presel_genes <- expressionDataGeneList()
  }else{
    presel_genes <- local(input$expressionDataGeneSelect)
  }
  
  presel_genes<- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws() 
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  if(length(sample_ids)>=900){
    sample_ids_filter_string <- c(paste(paste("sample_id", paste0("'", sample_ids[1:899], "'"), sep="="), collapse=" OR "))
    i<-900
    while(i <= length(sample_ids)){
      end<-i+899
      if(end>length(sample_ids)){
        end<-length(sample_ids)
      }
      sample_ids_filter_string <- c(sample_ids_filter_string, paste(paste("sample_id", paste0("'", sample_ids[i:end], "'"), sep="="), collapse=" OR "))
      i<-i+end
    }
  }else{
    sample_ids_filter_string <- paste(paste("sample_id", paste0("'", sample_ids, "'"), sep="="), collapse=" OR ")
  }

  query <- paste0("SELECT sample_id, gene_symbol, entrez_id, ensembl_id, expression_value FROM expression_data_values ",
                  "WHERE (", sample_ids_filter_string, ") ")
  
  if(!isTRUE(input$expressionDataCheckGeneAll)){
    gene_filter_str <- paste("entrez_id", paste0(presel_gene_entrez), sep="=")
    for(i in 1:length(presel_gene_entrez)){
      if(i==1){
        query_final <- paste0(query, "AND (", gene_filter_str[i], ") ")
      }else{
        query_final <- c(query_final, paste0(query, "AND (", gene_filter_str[i], ") "))
      }
    }
  }
  
  if(!isTRUE(input$expressionDataCheckCellLineAll) | length(presel_gene_entrez)<=50){
    print("query")
    df<- NULL
    for(z in 1:length(query_final)){
      # a chunk at a time
      res <- dbSendQuery(con_expression, query_final[z])
      i<-1
      while(!dbHasCompleted(res)){
        chunk <- dbFetch(res, n = 5000000)
        if(is.null(df)){
          df <- chunk
        }else{
          df <- df %>% rbind(chunk)
        }
        if(i %% 10==0){
          gc()
        }
        i<-i+1
      }
      dbClearResult(res)
    }
    df <- df %>%
      left_join(cellline_list_expressionData %>% dplyr::select(sample_id, species))
  }else{
    print("rsd")
    if(isTRUE(input$expressionDataCheckTissueAll)){
      if(input$expressionDataSpeciesSelect == "all"){
        df <- readRDS(file = paste0("expression_values_per_tissue/ALL_TISSUES.rds"))
      }
      if(input$expressionDataSpeciesSelect == "human"){
        df <- readRDS(file = paste0("expression_values_per_tissue/ALL_TISSUES_HUMAN.rds"))
      }
      if(input$expressionDataSpeciesSelect == "mouse"){
        df <- readRDS(file = paste0("expression_values_per_tissue/ALL_TISSUES_MOUSE.rds"))
      }
    }else{
      for(i in 1:length(presel_tissue)){
        if(i==1){
          df <- readRDS(file = paste0("expression_values_per_tissue/", presel_tissue[i], ".rds")) %>%
            filter(species %in% speciesList)
        }else{
          df <- df %>% 
            rbind(readRDS(file = paste0("expression_values_per_tissue/", presel_tissue[i], ".rds"))) %>%
            filter(species %in% speciesList)
        }
      }
    }
    if(!isTRUE(input$expressionDataCheckGeneAll)){
      df <- df %>%
        filter(entrez_id %in% presel_gene_entrez)
    }
  }
  
  df <- df %>% 
    distinct
  
  if(input$expressionDataSpeciesSelect == "all" & !is.null(df)){
    
    dict_joined <- dict_joined %>%
      left_join(df %>% 
                  dplyr::filter(species == "human") %>%
                  dplyr::select(gene_symbol, entrez_id, EnsemblID_human = ensembl_id) %>% distinct, by=c("EntrezID_human" = "entrez_id")) %>%
      left_join(df %>% 
                  dplyr::filter(species == "mouse") %>%
                  dplyr::select(gene_symbol, entrez_id, EnsemblID_mouse = ensembl_id) %>% distinct, by=c("EntrezID_mouse" = "entrez_id")) %>%
      distinct
    
    df_human <- df %>%
      dplyr::filter(species == "human") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_mouse, EnsemblID_mouse) %>% 
                  dplyr::filter(!is.na(Symbol_human)), by=c("gene_symbol" = "Symbol_human")) %>%
      dplyr::rename(Symbol_human = gene_symbol, EntrezID_human = entrez_id, EnsemblID_human = ensembl_id)
    
    df_mouse <- df %>%
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_human, EnsemblID_human) %>% 
                  dplyr::filter(!is.na(Symbol_mouse)), by=c("gene_symbol" = "Symbol_mouse")) %>%
      dplyr::rename(Symbol_mouse = gene_symbol, EntrezID_mouse = entrez_id, EnsemblID_mouse = ensembl_id)
    
    df <- df_human %>% rbind(df_mouse) %>% distinct
    
  }
  
  df
  
})

#create datatable out of dataframe
expressionDataDataTable <- eventReactive(input$expressionDataLoadButton,{
  
  df <- expressionDataDataFrame()
  
  if(!is.null(df)){
    
    brks <- seq(0, 20, length.out = 40)
    clrs <- round(seq(255, 5, length.out = (length(brks) + 1)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}

    nfreezeColumns <- 3
    
    if(input$expressionDataSpeciesSelect == "all"){

      dt <- df %>%
        dplyr::select(sample_id, expression_value, Symbol_human, EntrezID_human, EnsemblID_human, Symbol_mouse, EntrezID_mouse, EnsemblID_mouse) %>%
        pivot_wider(names_from=sample_id, values_from=expression_value) %>%
        arrange(Symbol_human, Symbol_mouse) %>%
        dplyr::select(Symbol_human, EntrezID_human, EnsemblID_human, Symbol_mouse, EntrezID_mouse, EnsemblID_mouse, everything())
      
      nfreezeColumns <- nfreezeColumns + 3
    }else{
      dt <- df %>%
        dplyr::select(sample_id, gene_symbol, entrez_id, ensembl_id, expression_value) %>%
        pivot_wider(names_from=sample_id, values_from=expression_value) %>%
        arrange(gene_symbol)
    }
    
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
    formatStyle(seq(nfreezeColumns+1, length(colnames(dt)),1),
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
    df
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
      arrange(entrez_id) %>%
      dplyr::mutate(gene = ifelse(is.na(gene_symbol), paste0("No symbol found (", entrez_id, ")"), paste0(gene_symbol , " (", entrez_id, ")"))) %>%
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
    dt <- expressionDataDataTable()
    if(!is.null(dt)){
      dt
    }
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
    
    sample_ids <- cellline_list_expressionData %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::filter(tissue_name %in%  presel_tissue) %>%
      dplyr::filter(cell_line_name %in% presel_cell_line) %>%
      dplyr::select(sample_id) %>%
      .$sample_id
    
    showModal(modalDialog(
      title = "WARNING!", 
      paste0("WARNING: You have selected ", 
             length(sample_ids), 
             " samples. Are you sure you want to load all samples for this selection?"),
      footer = tagList(
        modalButton("OK"),
        actionButton("expressionDataCancelModal", "Cancel")
      )
    ))
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

observeEvent(input$expressionDataCancelModal, {
  updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
  removeModal()
})

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
      
      if(input$expressionDataSpeciesSelect == "all"){
        
        dt <- df %>%
          dplyr::select(sample_id, expression_value, Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse) %>%
          pivot_wider(names_from=sample_id, values_from=expression_value) %>%
          arrange(Symbol_human, Symbol_mouse) %>%
          dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, everything())
        
        nfreezeColumns <- nfreezeColumns + 3
      }else{
        dt <- df %>%
          dplyr::select(sample_id, gene_symbol, entrez_id, ensembl_id, expression_value) %>%
          pivot_wider(names_from=sample_id, values_from=expression_value) %>%
          arrange(gene_symbol)
      }
      
      dt %>% write_tsv(file)
    }
  }
)