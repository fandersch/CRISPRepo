# ----------------------------------------------------------------------------
# slamseqData
# ----------------------------------------------------------------------------
slamseqData_geneInputFile <- reactiveValues(data = NULL)


slamseqDataUpdateText <- function(){
  output$slamseqDataInfo <- renderText({
    if(is.null(input$slamseqDataTissueSelect) & !isTRUE(input$slamseqDataCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      if(is.null(input$slamseqDataCellLineSelect) & !isTRUE(input$slamseqDataCheckCellLineAll)){
        invisible("INFO: Please select the cell line(s) in the right panel!")
      }else{
        if(is.null(input$slamseqDataGeneSelect) & !isTRUE(input$slamseqDataCheckGeneAll) & is.null(slamseqData_geneInputFile$data)){
          invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                    " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                    " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
        }else{
          invisible("INFO: Click Load data!")
        }
      }
    }
  })
}

#upon load display nothing
output$slamseqDataTable <- renderDataTable({
})

output$slamseqSampleMetaTable <- renderDataTable({
})

#query database and create dataframe
slamseqDataDataFrame <- reactive({
  #get selected species
  if(input$slamseqDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$slamseqDataSpeciesSelect
  }
  
  #get selected tissue
  if(isTRUE(input$slamseqDataCheckTissueAll)){
    presel_tissue <- slamseqDataTissueList()
  }else{
    presel_tissue <- local(input$slamseqDataTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$slamseqDataCheckCellLineAll)){
    presel_cell_line <- slamseqDataCellLineList()
  }else{
    presel_cell_line <- local(input$slamseqDataCellLineSelect)
  }
  
  sample_ids <- cellline_list_slamseq %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::filter(tissue_name %in%  presel_tissue) %>%
    dplyr::filter(cell_line_name %in% presel_cell_line) %>%
    dplyr::select(sample_id) %>%
    .$sample_id
  
  if(!is.null(slamseqData_geneInputFile$data)){
    presel_genes_buff <- slamseqDataGeneList()
    genes_fileUpload <- c(paste0("\\(", (slamseqData_geneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (slamseqData_geneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$slamseqDataCheckGeneAll)){
      presel_genes <- slamseqDataGeneList()
    }else{
      presel_genes <- local(input$slamseqDataGeneSelect)
    }
  }
  
  presel_genes<- presel_genes %>% strsplit(split="\\(|\\)")
  presel_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws() 
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  
  con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
  
  df <- con_slamseq %>%
    tbl("slamseq_counts")
  
  if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
    df <- df %>%
      dplyr::filter(sample_id %in% sample_ids)
  }
  
  if(!input$slamseqDataCheckGeneAll){
    df <- df %>%
      dplyr::filter(entrez_id %in% presel_gene_entrez)
  }
  
  df <- df %>%
    dplyr::select(sample_id, entrez_id, symbol, local(input$slamseqDataUnitSelect)) %>%
    distinct() %>%
    collect() %>%
    left_join(cellline_list_slamseq %>% dplyr::select(sample_id, sample_name, species) %>% distinct) %>%
    dplyr::rename(expression_value = input$slamseqDataUnitSelect)
  
  DBI::dbDisconnect(con_slamseq)
  
  if(input$slamseqDataSpeciesSelect == "all" & !is.null(df)){
    
    dict_joined <- dict_joined %>%
      left_join(df %>% 
                  dplyr::filter(species == "human") %>%
                  dplyr::select(symbol, entrez_id) %>% distinct, by=c("EntrezID_human" = "entrez_id")) %>%
      left_join(df %>% 
                  dplyr::filter(species == "mouse") %>%
                  dplyr::select(symbol, entrez_id) %>% distinct, by=c("EntrezID_mouse" = "entrez_id")) %>%
      distinct
    
    df_human <- df %>%
      dplyr::filter(species == "human") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_mouse) %>% 
                  dplyr::filter(!is.na(Symbol_human)), by=c("symbol" = "Symbol_human")) %>%
      dplyr::rename(Symbol_human = symbol, EntrezID_human = entrez_id)
    
    df_mouse <- df %>%
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined %>% dplyr::select(Symbol_human, Symbol_mouse, EntrezID_human) %>% 
                  dplyr::filter(!is.na(Symbol_mouse)), by=c("symbol" = "Symbol_mouse")) %>%
      dplyr::rename(Symbol_mouse = symbol, EntrezID_mouse = entrez_id)
    
    df <- df_human %>% rbind(df_mouse) %>% distinct
    #clean up
    df_human <- NULL
    df_mouse <- NULL
  }
  gc()
  df
  
})

#create datatable out of dataframe
slamseqDataDataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqDataDataFrame()
  
  if(!is.null(df)){

    nfreezeColumns <- 2
    
    if(input$slamseqDataSpeciesSelect == "all"){

      dt <- df %>%
        dplyr::select(sample_name, expression_value, Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse) %>%
        pivot_wider(names_from=sample_name, values_from=expression_value) %>%
        arrange(Symbol_human, Symbol_mouse) %>%
        dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, everything())
      
      nfreezeColumns <- nfreezeColumns + 2
    }else{
      dt <- df %>%
        dplyr::select(sample_name, symbol, entrez_id, expression_value) %>%
        pivot_wider(names_from=sample_name, values_from=expression_value) %>%
        arrange(symbol)
    }

    values<-df$expression_value
    max_value <- max(values, na.rm=T)
    min_value <- min(values, na.rm=T)

    brks <- exp(seq(0, log2(max_value), length.out = 40))
    clrs <- round(seq(255, 5, length.out = (length(brks) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
    
    df<-NULL
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
    
    if(!is.null(input$slamseqDataGeneSelect) | isTRUE(input$slamseqDataCheckGeneAll | !is.null(input$slamseqData_inputFile))){
      if(input$slamseqDataUnitSelect == "RPMu"){
        output$slamseqDataInfo <- renderText({
          "Info: Loading completed! Table shows RPMu values (expression dynamic).\n
          
          RPMus are calculated accordingly:\n
            nonTcReadCount = readCount - tcReadCount\n
            RPMu = (tcReadCount / sum(nonTcReadCount)) * 10^6"
        })
      }
      if(input$slamseqDataUnitSelect == "RPM"){
        output$slamseqDataInfo <- renderText({
          "Info: Loading completed! Table shows RPM values (expression level)."
        })
      }
      if(input$slamseqDataUnitSelect == "count"){
        output$slamseqDataInfo <- renderText({
          "Info: Loading completed! Table shows raw read count values (expression dynamic)."
        })
      }
      if(input$slamseqDataUnitSelect == "TCcount"){
        output$slamseqDataInfo <- renderText({
          "Info: Loading completed! Table shows raw T>C counts (expression level)."
        })
      }
    }
    gc()
    #display datatable
    dt
  }else{
    if(!is.null(input$slamseqDataGeneSelect) | isTRUE(input$slamseqDataCheckGeneAll)){
      output$slamseqDataInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      slamseqDataUpdateText()
    }
  }
})

#query database and create dataframe
slamseqSampleMetaDataFrame <- reactive({
  
  #get selected species
  if(input$slamseqDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$slamseqDataSpeciesSelect
  }
  
  #get selected tissue
  if(isTRUE(input$slamseqDataCheckTissueAll)){
    presel_tissue <- slamseqDataTissueList()
  }else{
    presel_tissue <- local(input$slamseqDataTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$slamseqDataCheckCellLineAll)){
    presel_cell_line <- slamseqDataCellLineList()
  }else{
    presel_cell_line <- local(input$slamseqDataCellLineSelect)
  }
  
  df <- cellline_list_slamseq %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::filter(tissue_name %in%  presel_tissue) %>%
    dplyr::filter(cell_line_name %in% presel_cell_line)
  
  df
})

#create datatable out of dataframe
slamseqSampleMetaDataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqSampleMetaDataFrame()
  if (nrow(df) > 0) {
    nfreezeColumns <- 1
    
    dt <- df %>%
      DT::datatable(extensions = c('FixedColumns','FixedHeader'), 
                    options = list(
                      autoWidth = FALSE,
                      headerCallback = JS(headerCallback),
                      scrollX=TRUE,
                      fixedColumns = list(leftColumns = nfreezeColumns),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      pageLength = 10,
                      lengthMenu = c(10, 25, 50, 100, 200),
                      searchHighlight = TRUE
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
      formatStyle(names(df),"white-space"="nowrap")
    
    
    output$slamseqDataInfo <- renderText({
      "Info: Loading completed!"
    })
    
    #display datatable
    dt 
  }else{
    output$slamseqDataInfo <- renderText({
      "WARNING: No data found!"
    })
  }
})

#----------------------------------------------------------------------------
# Selectbox Lists
#----------------------------------------------------------------------------

slamseqDataTissueList <- reactive({
 
  if(input$slamseqDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$slamseqDataSpeciesSelect
  }
  
  cellline_list_slamseq %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

slamseqDataCellLineList <- reactive({
  if(input$slamseqDataSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$slamseqDataSpeciesSelect
  }
  
  preselTissue = cellline_list_slamseq %>%
    dplyr::filter(species %in%  speciesList) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    .$tissue_name
  
  if(!isTRUE(input$slamseqDataCheckTissueAll) & !is.null(input$slamseqDataTissueSelect)){
    preselTissue = cellline_list_slamseq %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::filter(tissue_name %in% input$slamseqDataTissueSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
  }
  
  cellline_list_slamseq %>%
    dplyr::filter(species %in% speciesList, tissue_name %in% preselTissue) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
})

slamseqDataGeneList <- reactive({
  if(!isTRUE(input$slamseqDataCheckCellLineAll) & (is.null(input$slamseqDataCellLineSelect))){
    NULL
  }else{
    if(input$slamseqDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$slamseqDataSpeciesSelect
    }
    
    gene_list_slamseq %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::select(symbol, entrez_id) %>%
      arrange(entrez_id) %>%
      dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
      .$gene
  }
})


#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
observe(
  if(input$tabs == "slamseqDataSidebar"){
    loadslamseqDataTissueList <<- T
    if(input$slamseqDataSpeciesSelect == ""){
      select = "human"
    }else{
      select = input$slamseqDataSpeciesSelect
    }
    if(input$slamseqDataUnitSelect == ""){
      select_unit = "RPMu"
    }else{
      select_unit = input$slamseqDataUnitSelect
    }
    updateRadioButtons(session, 'slamseqDataUnitSelect', choices = list("RPMus" = "RPMu", "RPMs" = "RPM", "TCcounts" = "TCcount", "Counts" = "count"), selected = select_unit, inline = T)
    updateRadioButtons(session, 'slamseqDataSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = select, inline = T)
  }
)

observeEvent(input$slamseqDataLoadButton, {
  output$slamseqDataTable <- renderDataTable({
    dt <- slamseqDataDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$slamseqSampleMetaTable <- renderDataTable({
    dt <- slamseqSampleMetaDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
})

observeEvent(input$slamseqDataSpeciesSelect, {
  #update other species selects
  updateSpecies(input$slamseqDataSpeciesSelect)
  #select checkbox tissue
  updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'slamseqDataTissueSelect', choices = slamseqDataTissueList(), server = TRUE)
  #update cell linne selectbox
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = slamseqDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  #disable laod button
  disable("slamseqDataLoadButton")
  slamseqDataUpdateText()
})

observeEvent(input$slamseqDataTissueSelect, {
  if(!is.null(input$slamseqDataTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = slamseqDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$slamseqDataCheckTissueAll) | (!is.null(input$slamseqDataTissueSelect))){
    enable("slamseqDataCellLineSelect")
    enable("slamseqDataCheckCellLineAll")
  }else{
    disable("slamseqDataCellLineSelect")
    disable("slamseqDataCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCheckTissueAll, {
  if(isTRUE(input$slamseqDataCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'slamseqDataTissueSelect', choices = slamseqDataTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = slamseqDataCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$slamseqDataCheckTissueAll) | (!is.null(input$slamseqDataTissueSelect))){
    enable("slamseqDataCellLineSelect")
    enable("slamseqDataCheckCellLineAll")
  }else{
    disable("slamseqDataCellLineSelect")
    disable("slamseqDataCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCellLineSelect, {
  if(!is.null(input$slamseqDataCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$slamseqDataCheckCellLineAll) | (!is.null(input$slamseqDataCellLineSelect))) & (isTRUE(input$slamseqDataCheckTissueAll) | (!is.null(input$slamseqDataTissueSelect)))) {
    enable("slamseqDataGeneSelect")
    enable("slamseqData_inputFile")
    enable("slamseqDataCheckGeneAll")
  }else{
    disable("slamseqDataGeneSelect")
    disable("slamseqData_inputFile")
    disable("slamseqDataCheckGeneAll")
  }
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCheckCellLineAll, {
  if(isTRUE(input$slamseqDataCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = slamseqDataCellLineList(), server = TRUE)
    
    #get selected species
    if(input$slamseqDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$slamseqDataSpeciesSelect
    }
    
    #get selected tissue
    if(isTRUE(input$slamseqDataCheckTissueAll)){
      presel_tissue <- slamseqDataTissueList()
    }else{
      presel_tissue <- local(input$slamseqDataTissueSelect)
    }
    
    #get selected cell line
    if(isTRUE(input$slamseqDataCheckCellLineAll)){
      presel_cell_line <- slamseqDataCellLineList()
    }else{
      presel_cell_line <- local(input$slamseqDataCellLineSelect)
    }
    
    sample_ids <- cellline_list_slamseq %>%
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
        actionButton("slamseqDataCancelModal", "Cancel")
      )
    ))
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$slamseqDataCheckCellLineAll) | (!is.null(input$slamseqDataCellLineSelect))) & (isTRUE(input$slamseqDataCheckTissueAll) | (!is.null(input$slamseqDataTissueSelect)))) {
    enable("slamseqDataGeneSelect")
    enable("slamseqData_inputFile")
    enable("slamseqDataCheckGeneAll")
  }else{
    disable("slamseqDataGeneSelect")
    disable("slamseqData_inputFile")
    disable("slamseqDataCheckGeneAll")
  }
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCancelModal, {
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  removeModal()
})

observeEvent(input$slamseqDataGeneSelect, {
  if(!is.null(input$slamseqDataGeneSelect)){
    updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
    reset('slamseqData_inputFile')
    slamseqData_geneInputFile$data <- NULL
  }
  if(isTRUE(input$slamseqDataCheckGeneAll) | (!is.null(input$slamseqDataGeneSelect)) | !is.null(slamseqData_geneInputFile$data)){
    enable("slamseqDataLoadButton")
  }else{
    disable("slamseqDataLoadButton")
  }
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqData_inputFile, {
  if(!is.null(input$slamseqData_inputFile)){
    updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
    updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
    req(input$slamseqData_inputFile)
    slamseqData_geneInputFile$data <- read_tsv(input$slamseqData_inputFile$datapath, col_names = F)
  }else{
    slamseqData_geneInputFile$data <- NULL
  }
  if(isTRUE(input$slamseqDataCheckGeneAll) | (!is.null(input$slamseqDataGeneSelect)) | !is.null(slamseqData_geneInputFile$data)){
    enable("slamseqDataLoadButton")
  }else{
    disable("slamseqDataLoadButton")
  }
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCheckGeneAll, {
  if(isTRUE(input$slamseqDataCheckGeneAll)){
    geneList <- slamseqDataGeneList()
    updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = geneList, server = TRUE)
    reset('slamseqData_inputFile')
    slamseqData_geneInputFile$data <- NULL
  }
  if(isTRUE(input$slamseqDataCheckGeneAll) | (!is.null(input$slamseqDataGeneSelect)) | !is.null(slamseqData_geneInputFile$data)){
    enable("slamseqDataLoadButton")
  }else{
    disable("slamseqDataLoadButton")
  }
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$slamseqDataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_data_", local(input$slamseqDataUnitSelect), ".txt")
  },
  content = function(file) {
    df <- slamseqDataDataFrame()
    
    if (nrow(df) > 0) {
      
      if(input$slamseqDataSpeciesSelect == "all"){
        
        dt <- df %>%
          dplyr::select(sample_name, expression_value, Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse) %>%
          pivot_wider(names_from=sample_name, values_from=expression_value) %>%
          arrange(Symbol_human, Symbol_mouse) %>%
          dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, everything())
        
      }else{
        dt <- df %>%
          dplyr::select(sample_name, symbol, entrez_id, expression_value) %>%
          pivot_wider(names_from=sample_name, values_from=expression_value) %>%
          arrange(symbol)
      }
      #clean up
      df<-NULL
      gc()
      
      dt %>% write_tsv(file)
    }
  }
)

# output$slamseqDataButtonDownloadPrimaryTables <- downloadHandler(
#   filename = function() {
#     "expression_data.zip"
#   },
#   content = function(file) {
#     shiny::withProgress(
#       message = paste0("Downloading", input$dataset, " Data"),
#       value = 0,
#       {
#         files <- NULL;
#         if("Human" %in% input$slamseqDataDownloadPrimaryTablesCheck){
#           if("read_count" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_counts_human_spread.tsv"
#           }
#           if("tpm" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_tpms_human_spread.tsv"
#           }
#           if("log2_TMM" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_log2_TMM_human_spread.tsv"
#           }
#           
#           files <- c(fileName,files)
#         }
#         if("Mouse" %in% input$slamseqDataDownloadPrimaryTablesCheck){
#           if("read_count" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_counts_mouse_spread.tsv"
#           }
#           if("tpm" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_tpms_mouse_spread.tsv"
#           }
#           if("log2_TMM" %in% input$slamseqDataUnitSelect){
#             fileName <- "expression_values_per_tissue/all_tissues_log2_TMM_mouse_spread.tsv"
#           }
#           files <- c(fileName,files)
#         }
#         shiny::incProgress(1/2)
#         if(!is.null(files)){
#           #create the zip file
#           zip(file,files, compression_level = 2)
#         }else{
#           NULL
#         }
#       }
#     )
#   }
# )