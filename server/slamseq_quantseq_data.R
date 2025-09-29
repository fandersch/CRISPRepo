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
output$slamseqTotalDataTable <- renderDataTable({
})

output$slamseqDynamicsDataTable <- renderDataTable({
})

output$slamseqSampleMetaTable <- renderDataTable({
})

output$slamseqDeseq2Table <- renderDataTable({
})

output$slamseqDeseq2MetaTable <- renderDataTable({
})

#query database and create dataframe
slamseqDataTotalDataFrame <- reactive({

  sample_ids <- slamseqSampleMetaDataFrame() %>%
    .$sample_id %>%
    unique
  
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
  
  df_quant <- data.frame()
  if("quant" %in% input$slamseqDataIncludeDataset){
    
    con_quantseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/quantseq.db")
    
    df_quant <- con_quantseq %>%
      tbl("quantseq_sample_counts")
    
    if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
      df_quant <- df_quant %>%
        dplyr::filter(sample_id %in% sample_ids)
    }
    
    if(!input$slamseqDataCheckGeneAll){
      df_quant <- df_quant %>%
        dplyr::filter(entrez_id %in% presel_gene_entrez)
    }
    
    unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "TPM", "count")
    
    df_quant <- df_quant %>%
      dplyr::select(sample_id, entrez_id, symbol, ends_with(unit)) %>%
      distinct() %>%
      collect() %>%
      mutate(sample_id = as.character(sample_id)) %>%
      left_join(sample_list_quantseq %>% dplyr::select(sample_id, sample_name, species) %>% distinct) %>%
      dplyr::rename(expression_value = ends_with(unit))
    
    DBI::dbDisconnect(con_quantseq)
  }
  
  df_slam <- data.frame()
  if("slam" %in% input$slamseqDataIncludeDataset){
  
    con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
    
    df_slam <- con_slamseq %>%
      tbl("slamseq_sample_counts")
    
    if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
      df_slam <- df_slam %>%
        dplyr::filter(sample_id %in% sample_ids)
    }
    
    if(!input$slamseqDataCheckGeneAll){
      df_slam <- df_slam %>%
        dplyr::filter(entrez_id %in% presel_gene_entrez)
    }
    
    unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "RPM", "TCcount")
    
    df_slam <- df_slam %>%
      dplyr::select(sample_id, entrez_id, symbol, ends_with(unit)) %>%
      distinct() %>%
      collect() %>%
      mutate(sample_id = as.character(sample_id)) %>%
      left_join(sample_list_slamseq %>% dplyr::select(sample_id, sample_name, species) %>% distinct) %>%
      dplyr::rename(expression_value = ends_with(unit))
    
    DBI::dbDisconnect(con_slamseq)
  }
  
  df <- df_quant %>%
    bind_rows(df_slam)
  
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
slamseqDataTotalDataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqDataTotalDataFrame()
  
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
    
    max_value_quant <- 0
    if("quant" %in% input$slamseqDataIncludeDataset){
      unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "TPM", "count")
      #get value between max and average
      con_quantseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/quantseq.db")
      values <- con_quantseq %>%
        tbl("quantseq_sample_counts") %>%
        summarise(mean = mean(!!sym(unit), na.rm = T),
                  max = max(!!sym(unit), na.rm = T)) %>%
        collect()
      
      max_value_quant <- (as.numeric(values$max) + as.numeric(values$mean))/2
      
      DBI::dbDisconnect(con_quantseq)
      
    }
    
    max_value_slam <- 0
    if("slam" %in% input$slamseqDataIncludeDataset){
      unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "RPM", "TCcount")
      #get value between max and average
      con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
      values <- con_slamseq %>%
        tbl("slamseq_sample_counts") %>%
        summarise(mean = mean(!!sym(unit), na.rm = T),
                  max = max(!!sym(unit), na.rm = T)) %>%
        collect()
      
      max_value_slam <- (as.numeric(values$max) + as.numeric(values$mean))/2
      
      DBI::dbDisconnect(con_quantseq)
    }
    
    max_value <- ifelse(max_value_quant >= max_value_slam, max_value_quant, max_value_slam)

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
      if(input$slamseqDataUnitSelect == "rpm"){
        output$slamseqDataInfo <- renderText({
          paste0(
            "Info: Loading completed!", HTML('<br/>'), HTML('<br/>'), 
            "Expression Levels Table shows RPM values.", HTML('<br/>'), HTML('<br/>'),
            "Expression Dynamics Table shows RPMu values.", HTML('<br/>'), HTML('<br/>'),
            "RPMus are calculated accordingly:", HTML('<br/>'),
            "nonTcReadCount = readCount - tcReadCount", HTML('<br/>'),
            "RPMu = (tcReadCount / sum(nonTcReadCount)) * 10^6"
          )
        })
      }
      if(input$slamseqDataUnitSelect == "count"){
        output$slamseqDataInfo <- renderText({
          paste0(
            "Info: Loading completed!", HTML('<br/>'), HTML('<br/>'), 
            "Expression Levels Table shows raw count values.", HTML('<br/>'), HTML('<br/>'),
            "Expression Dynamics Table shows raw T>C counts values.", HTML('<br/>'), HTML('<br/>')
          )
        })
      }
    }
    gc()
    #display datatable
    dt
  }else{
    data.frame() %>% 
      datatable()
  }
})

#query database and create dataframe
slamseqDataDynamicsDataFrame <- reactive({
  
  df <- data.frame()
  if("slam" %in% input$slamseqDataIncludeDataset){
    
    sample_ids <- slamseqSampleMetaDataFrame() %>%
      .$sample_id %>%
      unique
    
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
      tbl("slamseq_sample_counts")
    
    if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
      df <- df %>%
        dplyr::filter(sample_id %in% sample_ids)
    }
    
    if(!input$slamseqDataCheckGeneAll){
      df <- df %>%
        dplyr::filter(entrez_id %in% presel_gene_entrez)
    }
    
    unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "RPMu", "TCcount")
    
    df <- df %>%
      dplyr::select(sample_id, entrez_id, symbol, matches(unit)) %>%
      distinct() %>%
      collect() %>%
      mutate(sample_id = as.character(sample_id)) %>%
      left_join(sample_list_slamseq %>% dplyr::select(sample_id, sample_name, species) %>% distinct) %>%
      dplyr::rename(expression_value =  matches(unit))
    
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
  }
  df
  
})

#create datatable out of dataframe
slamseqDataDynamicsDataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqDataDynamicsDataFrame()
  
  if("slam" %in% input$slamseqDataIncludeDataset & !is.null(df) & nrow(df) > 0){
    
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
    
    #get value between max and averag
    unit <- ifelse(local(input$slamseqDataUnitSelect) == "rpm", "RPMu", "TCcount")
    
    con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
    values <- con_slamseq %>%
      tbl("slamseq_sample_counts") %>%
      summarise(mean = mean(!!sym(unit), na.rm = T),
                max = max(!!sym(unit), na.rm = T)) %>%
      collect()
    
    max_value <- (as.numeric(values$max) + as.numeric(values$mean))/2
    
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
      if(input$slamseqDataUnitSelect == "rpm"){
        output$slamseqDataInfo <- renderText({
          paste0(
            "Info: Loading completed!", HTML('<br/>'), HTML('<br/>'), 
            "Expression Levels Table shows RPM values.", HTML('<br/>'), HTML('<br/>'),
            "Expression Dynamics Table shows RPMu values.", HTML('<br/>'), HTML('<br/>'),
            "RPMus are calculated accordingly:", HTML('<br/>'),
            "nonTcReadCount = readCount - tcReadCount", HTML('<br/>'),
            "RPMu = (tcReadCount / sum(nonTcReadCount)) * 10^6"
          )
        })
      }
      if(input$slamseqDataUnitSelect == "count"){
        output$slamseqDataInfo <- renderText({
          paste0(
            "Info: Loading completed!", HTML('<br/>'), HTML('<br/>'), 
            "Expression Levels Table shows raw count values.", HTML('<br/>'), HTML('<br/>'),
            "Expression Dynamics Table shows raw T>C counts values.", HTML('<br/>'), HTML('<br/>')
          )
        })
      }
    }
    gc()
    #display datatable
    dt
  }else{
    data.frame() %>% 
      datatable()
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
    presel_cell_line <- slamseqDataCellLineList() %>%
      .$cell_line_name %>% 
      unique
  }else{
    presel_cell_line <- local(input$slamseqDataCellLineSelect)
  }
  
  df <- slamseqDataCellLineList() %>%
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
    
    #display datatable
    dt 
  }else{
    data.frame() %>% 
      datatable()
  }
})


#query database and create dataframe
slamseqDeseq2DataFrame <- reactive({
  
  contrast_ids <- slamseqDeseq2MetaDataFrame() %>%
    .$contrast_id %>%
    unique
  
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
  
  df_quant <- data.frame()
  if("quant" %in% input$slamseqDataIncludeDataset){
    con_quantseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/quantseq.db")
    
    df_quant <- con_quantseq %>%
      tbl("quantseq_deseq2_results")
    
    if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
      df_quant <- df_quant %>%
        dplyr::filter(contrast_id %in% contrast_ids)
    }
    
    if(!input$slamseqDataCheckGeneAll){
      df_quant <- df_quant %>%
        dplyr::filter(entrez_id %in% presel_gene_entrez)
    }
    
    df_quant <- df_quant %>%
      distinct() %>%
      collect() %>%
      mutate(pvalue = as.numeric(formatC(pvalue, format = "e", digits = 2)),
             padj = as.numeric(formatC(padj, format = "e", digits = 2))) %>%
      left_join(slamseqDeseq2MetaDataFrame() %>% dplyr::select(contrast_id, species) %>% distinct)
    
    DBI::dbDisconnect(con_quantseq)
  }
  
  df_slam <- data.frame()
  if("slam" %in% input$slamseqDataIncludeDataset){
    con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
    
    df_slam <- con_slamseq %>%
      tbl("slamseq_deseq2_results")
    
    if(!input$slamseqDataCheckTissueAll | !input$slamseqDataCheckCellLineAll){
      df_slam <- df_slam %>%
        dplyr::filter(contrast_id %in% contrast_ids)
    }
    
    if(!input$slamseqDataCheckGeneAll){
      df_slam <- df_slam %>%
        dplyr::filter(entrez_id %in% presel_gene_entrez)
    }
    
    df_slam <- df_slam %>%
      distinct() %>%
      collect() %>%
      mutate(pvalue = as.numeric(formatC(pvalue, format = "e", digits = 2)),
             padj = as.numeric(formatC(padj, format = "e", digits = 2))) %>%
      left_join(slamseqDeseq2MetaDataFrame() %>% dplyr::select(contrast_id, species) %>% distinct)
    
    DBI::dbDisconnect(con_slamseq)
  }
  
  df <- df_quant %>%
    bind_rows(df_slam)
  
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
  df %>% dplyr::select(-species)
  
})


#create datatable out of dataframe
slamseqDeseq2DataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqDeseq2DataFrame()
  
  if (nrow(df) > 0) {
    
    nfreezeColumns <- 2
    
    if(input$slamseqDataSpeciesSelect == "all"){
      nfreezeColumns <- nfreezeColumns + 2
    }
    
    columns_wider <- colnames(df)
    columns_wider <- columns_wider[!columns_wider %in% c("contrast_id", "symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse")]
    
    df <- df %>%
      pivot_wider(names_from=contrast_id, values_from=c(columns_wider), names_glue = "{contrast_id}_{.value}") %>%
      arrange(dplyr::across(starts_with("symbol")))
    
    #order columns
    colnames(df) <- colnames(df) %>%
      str_replace("log2FoldChange$", "1log2FoldChange") %>%
      str_replace("pvalue$", "2pvalue") %>%
      str_replace("padj$", "3padj") %>%
      str_replace("lfcSE$", "4lfcSE") %>%
      str_replace("baseMean$","5baseMean") %>%
      str_replace("treatment_mean$", "6treatment_mean") %>%
      str_replace("control_mean$", "7control_mean")
    
    df <- df %>%
      dplyr::select(order(colnames(df))) %>%
      dplyr::select(any_of(c("symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse")), everything())
    
    colnames(df) <- colnames(df) %>%
      str_replace("1log2FoldChange$", "log2FoldChange") %>%
      str_replace("2pvalue$", "pvalue") %>%
      str_replace("3padj$", "padj") %>%
      str_replace("4lfcSE$", "lfcSE") %>%
      str_replace("5baseMean$","baseMean") %>%
      str_replace("6treatment_mean$", "treatment_mean") %>%
      str_replace("7control_mean$", "control_mean")
    
    #get value between max and averag
    con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
    values <- con_slamseq %>%
      tbl("slamseq_deseq2_results") %>%
      summarise(mean = mean(!!sym("log2FoldChange"), na.rm = T),
                max = max(!!sym("log2FoldChange"), na.rm = T),
                min = min(!!sym("log2FoldChange"), na.rm = T)) %>%
      collect()
    
    max_value <- (as.numeric(values$max) + as.numeric(values$mean))/2
    min_value <- (as.numeric(values$min) + as.numeric(values$mean))/2
    
    brks_smaller <- seq(min_value, 0, .05)
    brks_bigger <- seq(0, max_value, .05)
    
    clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
    clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
      {paste0("rgb(", ., ",", ., ",255)")}
    
    brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
    clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
    
    colorInterval <- length(columns_wider)
    
    dt <- df %>%
      DT::datatable(class = "display nowrap",
                    extensions = c('FixedColumns','FixedHeader'),
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
      formatStyle(seq(nfreezeColumns+1, length(colnames(df)),colorInterval),
                backgroundColor = styleInterval(brks, clrs))
    
    gc()
    #display datatable
    dt
  }else{
    data.frame() %>% 
      datatable()
  }
})



#query database and create dataframe
slamseqDeseq2MetaDataFrame <- reactive({
  
  sample_metadata <- slamseqSampleMetaDataFrame()
  
  sample_ids <- sample_metadata %>%
    .$sample_id
  
  contrast_ids <- deseq2_list_slamseq %>%
    bind_rows(deseq2_list_quantseq) %>%
    separate_longer_delim(cols = treatment_sample_ids, delim = ";") %>%
    separate_longer_delim(cols = control_sample_ids, delim = ";") %>%
    filter(treatment_sample_ids %in% sample_ids | control_sample_ids %in% sample_ids) %>%
    .$contrast_id %>%
    unique
  
  deseq2_list_slamseq %>%
    bind_rows(deseq2_list_quantseq) %>%
    filter(contrast_id %in% contrast_ids)
})

#create datatable out of dataframe
slamseqDeseq2MetaDataTable <- eventReactive(input$slamseqDataLoadButton,{
  
  df <- slamseqDeseq2MetaDataFrame()
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
    
    #display datatable
    dt 
  }else{
    data.frame() %>% 
      datatable()
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
  
  sample_list <- sample_list_slamseq %>%
    bind_rows(sample_list_quantseq)
  
  if(!("quant" %in% input$slamseqDataIncludeDataset)){
    sample_list <- sample_list %>%
      anti_join(sample_list_quantseq)
  }
  
  if(!("slam" %in% input$slamseqDataIncludeDataset)){
    sample_list <- sample_list %>%
      anti_join(sample_list_slamseq)
  }
  
  if(isTRUE(input$slamseqDataCheckOnlyIncludeBaseline)){
    sample_list <- sample_list %>%
      filter(is_baseline_reference == "yes")
  }
  
  sample_list %>%
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
  
  sample_list <- sample_list_slamseq %>%
    bind_rows(sample_list_quantseq)
  
  if(!("quant" %in% input$slamseqDataIncludeDataset)){
    sample_list <- sample_list %>%
      anti_join(sample_list_quantseq)
  }
  
  if(!("slam" %in% input$slamseqDataIncludeDataset)){
    sample_list <- sample_list %>%
      anti_join(sample_list_slamseq)
  }
  
  if(isTRUE(input$slamseqDataCheckOnlyIncludeBaseline)){
    sample_list <- sample_list %>%
      filter(is_baseline_reference == "yes")
  }
  
  preselTissue = sample_list %>%
    dplyr::filter(species %in%  speciesList) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    .$tissue_name
  
  if(!isTRUE(input$slamseqDataCheckTissueAll) & !is.null(input$slamseqDataTissueSelect)){
    preselTissue = sample_list %>%
      dplyr::filter(species %in% speciesList) %>%
      dplyr::filter(tissue_name %in% input$slamseqDataTissueSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
  }
  
  sample_list %>%
    dplyr::filter(species %in% speciesList, tissue_name %in% preselTissue)
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
      bind_rows(gene_list_quantseq) %>%
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
      select_unit = c("rpm")
    }else{
      select_unit = input$slamseqDataUnitSelect
    }
    updateRadioButtons(session, 'slamseqDataUnitSelect', choices = list("RPMs/RPMus" = "rpm", "Counts/TCcounts" = "count"), selected = select_unit, inline = T)
    updateRadioButtons(session, 'slamseqDataSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = select, inline = T)
  }
)

observeEvent(input$slamseqDataLoadButton, {
  output$slamseqTotalDataTable <- renderDataTable({
    dt <- slamseqDataTotalDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$slamseqDynamicsDataTable <- renderDataTable({
    dt <- slamseqDataDynamicsDataTable()
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
  
  output$slamseqDeseq2Table <- renderDataTable({
    dt <- slamseqDeseq2DataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$slamseqDeseq2MetaTable <- renderDataTable({
    dt <- slamseqDeseq2MetaDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
})

observeEvent(input$slamseqDataSpeciesSelect, {
  #select checkbox tissue
  updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'slamseqDataTissueSelect', choices = slamseqDataTissueList(), server = TRUE)
  if(length(slamseqDataTissueList()) >0){
    enable("slamseqDataCheckTissueAll")
  }else{
    disable("slamseqDataCheckTissueAll")
  }
  #update cell linne selectbox
  celllines <- slamseqDataCellLineList() %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
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

observeEvent(input$slamseqDataIncludeDataset, {
  #select checkbox tissue
  updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'slamseqDataTissueSelect', choices = slamseqDataTissueList(), server = TRUE)
  if(length(slamseqDataTissueList()) >0){
    enable("slamseqDataCheckTissueAll")
  }else{
    disable("slamseqDataCheckTissueAll")
  }
  #update cell linne selectbox
  celllines <- slamseqDataCellLineList() %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  #disable laod button
  disable("slamseqDataLoadButton")
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataCheckOnlyIncludeBaseline, {
  #select checkbox tissue
  updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'slamseqDataTissueSelect', choices = slamseqDataTissueList(), server = TRUE)
  if(length(slamseqDataTissueList()) >0){
    enable("slamseqDataCheckTissueAll")
  }else{
    disable("slamseqDataCheckTissueAll")
  }
  #update cell linne selectbox
  celllines <- slamseqDataCellLineList() %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'slamseqDataCheckCellLineAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'slamseqDataGeneSelect', choices = slamseqDataGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'slamseqDataCheckGeneAll', value = FALSE)
  #disable laod button
  disable("slamseqDataLoadButton")
  slamseqDataUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$slamseqDataTissueSelect, {
  if(!is.null(input$slamseqDataTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'slamseqDataCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  celllines <- slamseqDataCellLineList() %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
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
  celllines <- slamseqDataCellLineList() %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
  
  updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
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
    celllines <- slamseqDataCellLineList() %>%
      dplyr::select(cell_line_name) %>%
      arrange(cell_line_name) %>%
      .$cell_line_name
    
    updateSelectizeInput(session, 'slamseqDataCellLineSelect', choices = celllines, server = TRUE)
    
    sample_ids <- slamseqDataCellLineList() %>%
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

output$slamseqExpressionDataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_expression_data_", local(input$slamseqDataUnitSelect), ".txt")
  },
  content = function(file) {
    df <- slamseqDataTotalDataTable()$x$data
    
    if (nrow(df) > 0) {
      df %>% write_tsv(file)
    }
  }
)

output$slamseqExpressionDynamicsDataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_expression_data_", local(input$slamseqDataUnitSelect), ".txt")
  },
  content = function(file) {
    df <- slamseqDataDynamicsDataTable()$x$data
    
    if (nrow(df) > 0) {
      df %>% write_tsv(file)
    }
  }
)

output$slamseqSampleMetaDataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_sample_metadata.txt")
  },
  content = function(file) {
    df <- slamseqSampleMetaDataTable()$x$data
    
    if (nrow(df) > 0) {
      df %>% write_tsv(file)
    }
  }
)

output$slamseqDeseq2DataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_DEA_data.txt")
  },
  content = function(file) {
    df <- slamseqDeseq2DataTable()$x$data
    
    if (nrow(df) > 0) {
      df %>% write_tsv(file)
    }
  }
)

output$slamseqDeseq2MetaDataButtonDownload <- downloadHandler(
  filename = function() {
    paste0("crisprepo_slamseq_DEA_metadata_.txt")
  },
  content = function(file) {
    df <- slamseqDeseq2MetaDataTable()$x$data
  
    if (nrow(df) > 0) {
      df %>% write_tsv(file)
    }
  }
)