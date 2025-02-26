# ----------------------------------------------------------------------------
# patientMutation
# ----------------------------------------------------------------------------
patientMutationGeneInputFile <- reactiveValues(data = NULL)


patientMutationUpdateText <- function(){
  output$patientMutationInfo <- renderText({
    if(is.null(input$patientMutationCancerTypeSelect) & !isTRUE(input$patientMutationCheckCancerTypeAll)){
      invisible("INFO: Please select the cancer type(s) in the right panel!")
    }else{
      if(is.null(input$patientMutationGeneSelect) & !isTRUE(input$patientMutationCheckGeneAll) & is.null(patientMutationGeneInputFile$data)){
        invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                  " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                  " The file must have one gene per row (one single column): gene-symbol (no header)!"))
      }else{
        invisible("INFO: Click Load data!")
      }
    }
  })
}

#upon load display nothing
output$patientMutationDataTableAlterations <- renderDataTable({
})
output$patientMutationDataTableMutations <- renderDataTable({
})
output$patientMutationDataTableFusions <- renderDataTable({
})
output$patientMutationDataTableCNVs <- renderDataTable({
})

patientMutationAlterationsDataFrame <- reactive({
  
  #get selected CancerType
  if(isTRUE(input$patientMutationCheckCancerTypeAll)){
    presel_CancerType <- patient_cancer_types
  }else{
    presel_CancerType <- local(input$patientMutationCancerTypeSelect)
  }
  
  if(!is.null(patientMutationGeneInputFile$data)){
    presel_genes_buff <- patientMutationGeneList()
    genes_fileUpload <- c(patientMutationGeneInputFile$data$X1 %>% as.character)
    presel_genes_both <- grep(paste(paste0(" ", genes_fileUpload, " "),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$patientMutationCheckGeneAll)){
      presel_genes_both <- patientMutationGeneList()
    }else{
      presel_genes_both <- local(input$patientMutationGeneSelect)
    }
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  con_patient_data <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cBioportal_mutations_CNAs_fusions.db")
  
  #number patients per cancer type and dataset
  patients_per_cancerType <- con_patient_data %>%
    tbl("curated_set_non_redundant_sample_info") %>%
    dplyr::select(cancer_type, study_id, patient_id, dataset) %>%
    distinct %>%
    dplyr::filter(cancer_type %in% presel_CancerType) %>%
    collect()
  
  n_patients_per_cancerType_total <- patients_per_cancerType %>%
    dplyr::select(cancer_type, study_id, patient_id) %>%
    distinct %>%
    group_by(cancer_type) %>%
    summarize(n_patients_total=n()) %>%
    ungroup
    
  n_patients_per_cancerType_mutations <- patients_per_cancerType %>%
    filter(dataset == "mutations") %>%
    dplyr::select(cancer_type, study_id, patient_id) %>%
    distinct %>%
    group_by(cancer_type) %>%
    summarize(n_patients_total=n()) %>%
    ungroup
  
  n_patients_per_cancerType_structuralVariants <- patients_per_cancerType %>%
    filter(dataset == "structuralVariants") %>%
    dplyr::select(cancer_type, study_id, patient_id) %>%
    distinct %>%
    group_by(cancer_type) %>%
    summarize(n_patients_total=n()) %>%
    ungroup
  
  n_patients_per_cancerType_CNAs <- patients_per_cancerType %>%
    filter(dataset == "CNAs") %>%
    dplyr::select(cancer_type, study_id, patient_id) %>%
    distinct %>%
    group_by(cancer_type) %>%
    summarize(n_patients_total=n()) %>%
    ungroup
  
  #mutations
  patientMutation_mutation_data <- con_patient_data %>%
    tbl("mutations") %>%
    dplyr::filter(cancer_type %in% presel_CancerType, symbol %in% presel_genes | entrez_id %in% presel_entrez) %>%
    dplyr::select(symbol, entrez_id, cancer_type, patient_id, proteinChange) %>%
    distinct() %>%
    collect()
  
  patientMutation_mutation_data_per_cancerType <- patientMutation_mutation_data %>%
    dplyr::select(symbol, entrez_id, cancer_type, patient_id) %>%
    distinct() %>%
    inner_join(n_patients_per_cancerType_mutations) %>%
    group_by(symbol, entrez_id, cancer_type, n_patients_total_mutations=n_patients_total) %>%
    summarize(n_patients_gene_mutation=n(), percent_patients_gene_mutation=round(100*(n()/unique(n_patients_total_mutations)),2)) %>%
    ungroup()
  
  #fusions
  patientMutation_fusions_data <- con_patient_data %>%
    tbl("structuralVariants") %>%
    dplyr::filter(cancer_type %in% presel_CancerType, (symbol_1 %in% presel_genes | entrez_id_1 %in% presel_entrez) | (symbol_2 %in% presel_genes | entrez_id_2 %in% presel_entrez))  %>%
    dplyr::select(symbol_1,entrez_id_1, symbol_2, entrez_id_2, cancer_type, patient_id, eventInfo) %>%
    distinct() %>%
    collect()
  
  patientMutation_fusions_data_per_cancerType <- patientMutation_fusions_data %>% 
    dplyr::rename(entrez_id=entrez_id_1, symbol=symbol_1) %>% 
    bind_rows(patientMutation_fusions_data %>%
                dplyr::rename(entrez_id=entrez_id_2, symbol=symbol_2)) %>%
    distinct() %>%
    dplyr::select(-entrez_id_1, -symbol_1, -entrez_id_2, -symbol_2, -eventInfo) %>%
    dplyr::filter((symbol %in% presel_genes | entrez_id %in% presel_entrez)) %>%
    distinct() %>%
    inner_join(n_patients_per_cancerType_structuralVariants) %>%
    group_by(symbol, entrez_id, cancer_type, n_patients_total_fusions=n_patients_total) %>%
    summarize(n_patients_gene_fusion = n(), percent_patients_gene_fusion = round(100*(n()/unique(n_patients_total_fusions)),2)) %>%
    ungroup()
    
  #cna
  patientMutation_cna_data <- con_patient_data %>%
    tbl("CNAs") %>%
    dplyr::filter(symbol %in% presel_genes | entrez_id %in% presel_entrez) %>%
    dplyr::select(symbol, entrez_id, cancer_type, patient_id, alteration) %>%
    distinct() %>%
    collect()
  
  patientMutation_cna_data_per_cancerType <- patientMutation_cna_data %>%
    inner_join(n_patients_per_cancerType_CNAs) %>%
    group_by(symbol, entrez_id, cancer_type, n_patients_total_cna=n_patients_total) %>%
    summarize(n_patients_gene_cna = n(), percent_patients_gene_cna = round(100*(n()/unique(n_patients_total_cna)),2),
              n_patients_gene_cna_deletion=sum(alteration == -2), percent_patients_gene_cna_deletion=round(100*(n_patients_gene_cna_deletion/unique(n_patients_total_cna)),2),
              n_patients_gene_cna_amplification=sum(alteration == 2), percent_patients_gene_cna_amplification=round(100*(n_patients_gene_cna_amplification/unique(n_patients_total_cna)),2)) %>%
    ungroup()
  
  patientMutation_alteration_data <- patientMutation_mutation_data %>%
    full_join(patientMutation_cna_data) %>%
    full_join(patientMutation_fusions_data) %>%
    dplyr::select(symbol, entrez_id, cancer_type, patient_id) %>%
    distinct() %>%
    inner_join(n_patients_per_cancerType_total) %>%
    group_by(symbol, entrez_id, cancer_type, n_patients_total) %>%
    summarize(n_patients_gene_alteration = n(), percent_patients_gene_alteration = round(100*(n()/unique(n_patients_total)),2)) %>%
    ungroup()
  
  patientMutation_alteration_data_details <- patientMutation_alteration_data %>% 
    full_join(patientMutation_mutation_data_per_cancerType) %>%
    full_join(patientMutation_cna_data_per_cancerType) %>%
    full_join(patientMutation_fusions_data_per_cancerType)
  
  patientMutation_alteration_data <- patientMutation_alteration_data %>% 
    mutate(cancer_type = paste0(cancer_type, " (n=", n_patients_total, ")")) %>%
    filter(!is.na(symbol)) %>%
    filter(n_patients_total >= input$patientMutationMinPatients) %>%
    pivot_wider(id_cols = "symbol", names_from = "cancer_type", values_from = "percent_patients_gene_alteration", values_fill = 0) %>%
    column_to_rownames(var = "symbol") %>%
    filter(rowSums(.[,] >= 0.1) > 0) %>%
    select_if(colSums(. >= 0.1) > 0)
  
  DBI::dbDisconnect(con_patient_data)
  
  list(patientMutation_alteration_data, patientMutation_alteration_data_details)
})

patientMutationMutationsDataFrame <- reactive({
  
  #get selected CancerType
  if(isTRUE(input$patientMutationCheckCancerTypeAll)){
    presel_CancerType <- patient_cancer_types
  }else{
    presel_CancerType <- local(input$patientMutationCancerTypeSelect)
  }
  
  if(!is.null(patientMutationGeneInputFile$data)){
    presel_genes_buff <- patientMutationGeneList()
    genes_fileUpload <- c(patientMutationGeneInputFile$data$X1 %>% as.character)
    presel_genes_both <- grep(paste(paste0(" ", genes_fileUpload, " "),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$patientMutationCheckGeneAll)){
      presel_genes_both <- patientMutationGeneList()
    }else{
      presel_genes_both <- local(input$patientMutationGeneSelect)
    }
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  con_patient_data <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cBioportal_mutations_CNAs_fusions.db")
  
  patientMutation_mutation_data <- con_patient_data %>%
    tbl("mutations") %>%
    dplyr::filter(cancer_type %in% presel_CancerType, symbol %in% presel_genes | entrez_id %in% presel_entrez) %>%
    collect()
  
  DBI::dbDisconnect(con_patient_data)
  
  patientMutation_mutation_data
})

patientMutationFusionsDataFrame <- reactive({
  
  #get selected CancerType
  if(isTRUE(input$patientMutationCheckCancerTypeAll)){
    presel_CancerType <- patient_cancer_types
  }else{
    presel_CancerType <- local(input$patientMutationCancerTypeSelect)
  }
  
  if(!is.null(patientMutationGeneInputFile$data)){
    presel_genes_buff <- patientMutationGeneList()
    genes_fileUpload <- c(patientMutationGeneInputFile$data$X1 %>% as.character)
    presel_genes_both <- grep(paste(paste0(" ", genes_fileUpload, " "),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$patientMutationCheckGeneAll)){
      presel_genes_both <- patientMutationGeneList()
    }else{
      presel_genes_both <- local(input$patientMutationGeneSelect)
    }
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  con_patient_data <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cBioportal_mutations_CNAs_fusions.db")
  
  patientMutation_fusions_data <- con_patient_data %>%
    tbl("structuralVariants") %>%
    dplyr::filter(cancer_type %in% presel_CancerType, (symbol_1 %in% presel_genes | entrez_id_1 %in% presel_entrez) | (symbol_2 %in% presel_genes | entrez_id_2 %in% presel_entrez))  %>%
    collect()
  
  DBI::dbDisconnect(con_patient_data)
  
  patientMutation_fusions_data
})

patientMutationCNVsDataFrame <- reactive({
  
  #get selected CancerType
  if(isTRUE(input$patientMutationCheckCancerTypeAll)){
    presel_CancerType <- patient_cancer_types
  }else{
    presel_CancerType <- local(input$patientMutationCancerTypeSelect)
  }
  
  if(!is.null(patientMutationGeneInputFile$data)){
    presel_genes_buff <- patientMutationGeneList()
    genes_fileUpload <- c(patientMutationGeneInputFile$data$X1 %>% as.character)
    presel_genes_both <- grep(paste(paste0(" ", genes_fileUpload, " "),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$patientMutationCheckGeneAll)){
      presel_genes_both <- patientMutationGeneList()
    }else{
      presel_genes_both <- local(input$patientMutationGeneSelect)
    }
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  con_patient_data <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cBioportal_mutations_CNAs_fusions.db")
  
  patientMutation_cna_data <- con_patient_data %>%
    tbl("CNAs") %>%
    dplyr::filter(cancer_type %in% presel_CancerType, symbol %in% presel_genes | entrez_id %in% presel_entrez) %>%
    collect()
  
  DBI::dbDisconnect(con_patient_data)
  
  patientMutation_cna_data
})

#create plot out of dataframe
patientMutationHeatmapAlterations <- eventReactive(input$patientMutationLoadButton,{
  heatmap_matrix <- patientMutationAlterationsDataFrame()[[1]]
  
  if (nrow(heatmap_matrix) > 0) {
    output$patientMutationInfo <- renderText({"INFO: Loading completed!"})
    
    nrows <- nrow(heatmap_matrix)
    ncols <- ncol(heatmap_matrix)
    clust_rows<-ifelse(nrows>1,T,F)
    clust_cols<-ifelse(ncols>1,T,F)
    max_font <- 12
    min_font <- 8
    fontsize_row <- max(min(max_font, 50 / nrows), min_font)
    fontsize_col <- max(min(max_font, 50 / ncols), min_font)
    
    max_percent_color <- input$patientMutationMaxHeatmapColour
    breaks<- seq(0.0001, max_percent_color, max_percent_color/20)

    #plot heatmap
    p<-pheatmap(heatmap_matrix, 
                          cluster_rows=clust_rows,
                          cluster_cols=clust_cols,
                          color=colorRampPalette(c("white", "firebrick3" ))(22), 
                          breaks = c(0,breaks,max(heatmap_matrix)), 
                          legend_breaks = c(0,max_percent_color,floor(max(heatmap_matrix))),
                          fontsize_col=fontsize_col, 
                          fontsize_row = fontsize_row, 
                          angle_col = 90, main="Gene alterations per cancer type",
                          silent=T
               )
    list(p, nrow(heatmap_matrix))
  }else{
    NULL
  }
})


#create datatable out of dataframe
patientMutationDataTableAlterations <- eventReactive(input$patientMutationLoadButton,{
  
  df <- patientMutationAlterationsDataFrame()[[2]]
  
  if (nrow(df) > 0) {
    output$patientMutationInfo <- renderText({
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

patientMutationDataTableMutations <- eventReactive(input$patientMutationLoadButton,{
  
  df <- patientMutationMutationsDataFrame()
  
  if (nrow(df) > 0) {
    output$patientMutationInfo <- renderText({
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

patientMutationDataTableFusions <- eventReactive(input$patientMutationLoadButton,{
  
  df <- patientMutationFusionsDataFrame()
  
  if (nrow(df) > 0) {
    output$patientMutationInfo <- renderText({
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

patientMutationDataTableCNVs <- eventReactive(input$patientMutationLoadButton,{
  
  df <- patientMutationCNVsDataFrame()
  
  if (nrow(df) > 0) {
    output$patientMutationInfo <- renderText({
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

patientMutationGeneList <- reactive({
  presel_cancer_type <- patient_cancer_types
    
  if(!is.null(input$patientMutationCancerTypeSelect)){
    presel_cancer_type <- input$patientMutationCancerTypeSelect
  }
  
  patient_genes_all %>%
    filter(cancer_type %in% presel_cancer_type) %>%
    dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (" , entrez_id, " )"), paste0(" ", symbol , " ( ", entrez_id, " )"))) %>%
    dplyr::select(gene) %>%
    distinct %>%
    arrange(gene) %>%
    .$gene
})

#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
output$patientMutationHeatmapAlterations <- renderImage({
  p <- patientMutationHeatmapAlterations()[[1]]
  nrow <-  patientMutationHeatmapAlterations()[[2]]
  
  if(!is.null(dt)){
    container_width <- session$clientData$output_patientMutationHeatmapAlterations_width
    # if (is.null(container_width) || container_width == 0) {
    
    height <- (nrow*10) + 400 
    session$sendCustomMessage("updateImageHeight", list(height = height))
    
    # A temp file to save the output.
    outfile <- tempfile(fileext='.png')
    
    png(outfile, width=container_width*7, height=height*7, res=600)
    print(p)
    dev.off()
    # Return a list containing the filename
    list(src = outfile,
         width = container_width,
         height = height,
         alt = "This is alternate text")
  }
}, deleteFile = TRUE)

observeEvent(input$patientMutationLoadButton, {
  
  output$patientMutationDataTableAlterations <- renderDataTable({
    dt <- patientMutationDataTableAlterations()
    if(!is.null(dt)){
      dt
    }
  })

  output$patientMutationDataTableMutations <- renderDataTable({
    dt <- patientMutationDataTableMutations()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$patientMutationDataTableFusions <- renderDataTable({
    dt <- patientMutationDataTableFusions()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$patientMutationDataTableCNVs <- renderDataTable({
    dt <- patientMutationDataTableCNVs()
    if(!is.null(dt)){
      dt
    }
  })
})

observeEvent(input$patientMutationCancerTypeSelect, {
  if(!is.null(input$patientMutationCancerTypeSelect)){
    #unselect checkbox CancerType
    updateCheckboxInput(session, 'patientMutationCheckCancerTypeAll', value = FALSE)
    enable("patientMutationGeneSelect")
    enable("patientMutationGeneInputFile")
  }else{
    updateCheckboxInput(session, 'patientMutationCheckCancerTypeAll', value = TRUE)
  }
  
  updateSelectizeInput(session, 'patientMutationGeneSelect', choices = patientMutationGeneList(), server = TRUE)
  
  patientMutationUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$patientMutationCheckCancerTypeAll, {
  if(isTRUE(input$patientMutationCheckCancerTypeAll)){
    #reset CancerType selectbox
    updateSelectizeInput(session, 'patientMutationCancerTypeSelect', choices = patient_cancer_types, server = TRUE)
    enable("patientMutationGeneSelect")
    enable("patientMutationGeneInputFile")
  }else{
    if(is.null(input$patientMutationCancerTypeSelect)){
      disable("patientMutationGeneSelect")
      disable("patientMutationGeneInputFile")
    }
  }
  
  updateSelectizeInput(session, 'patientMutationGeneSelect', choices = patientMutationGeneList(), server = TRUE)
  
}, ignoreNULL = FALSE)

observeEvent(input$patientMutationCancelModal, {
  updateCheckboxInput(session, 'patientMutationCheckpatientMutationAll', value = FALSE)
  removeModal()
})

observeEvent(input$patientMutationGeneSelect, {
  if(!is.null(input$patientMutationGeneSelect)){
    reset('patientMutationGeneInputFile')
    patientMutationGeneInputFile$data <- NULL
  }
  if((!is.null(input$patientMutationGeneSelect)) | !is.null(patientMutationGeneInputFile$data)){
    enable("patientMutationLoadButton")
  }else{
    disable("patientMutationLoadButton")
  }
  patientMutationUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$patientMutationGeneInputFile, {
  if(!is.null(input$patientMutationGeneInputFile)){
    updateSelectizeInput(session, 'patientMutationGeneSelect', choices = patientMutationGeneList(), server = TRUE)
    req(input$patientMutationGeneInputFile)
    patientMutationGeneInputFile$data <- read_tsv(input$patientMutationGeneInputFile$datapath, col_names = F)
  }else{
    patientMutationGeneInputFile$data <- NULL
  }
  if((!is.null(input$patientMutationGeneSelect)) | !is.null(patientMutationGeneInputFile$data)){
    enable("patientMutationLoadButton")
  }else{
    disable("patientMutationLoadButton")
  }
  patientMutationUpdateText()
  
}, ignoreNULL = FALSE)


# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$patientMutationButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_patientMutation_meta_data.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading data"),
      value = 0,
      {
        heatmap <- patientMutationHeatmapAlterations()[[1]]
        nrow <- patientMutationHeatmapAlterations()[[2]]
        heatmap_matrix <- patientMutationAlterationsDataFrame()[[1]]
        alterations <- patientMutationAlterationsDataFrame()[[2]]
        mutations <- patientMutationMutationsDataFrame()
        fusions <- patientMutationFusionsDataFrame()
        cnvs <- patientMutationCNVsDataFrame()
        
        print("here")
        
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;
        
        if("Gene alterations heatmap" %in% input$patientMutationDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          plot <- heatmap
          fileName <- "patient_gene_alteration_heatmap.pdf"
          height <- (6+nrow/10)
          ggsave(fileName, plot, "pdf", width=9, height=height)
          files <- c(fileName,files)
          table <- heatmap_matrix
          fileName <- "patient_gene_alteration_heatmap_matrix.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(1/5)
        print("here2")
        if("Gene alterations" %in% input$patientMutationDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          table <- alterations
          fileName <- "patient_alterations.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(2/5)
        print("here3")
        if("Gene mutations" %in% input$patientMutationDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          table <- mutations
          fileName <- "patient_mutations.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(3/5)
        print("here4")
        if("Gene fusions" %in% input$patientMutationDownloadCheck & nrow(fusions)>0){
          #write each sheet to a csv file, save the name
          table <- fusions
          fileName <- "patient_fusions.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(4/5)
        print("here5")
        if("Gene CNVs" %in% input$patientMutationDownloadCheck & nrow(cnvs)>0){
          #write each sheet to a csv file, save the name
          table <- cnvs
          fileName <- "patient_cnvs.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(5/5)
        print("here6")
        if(!is.null(files)){
          #create the zip file
          zip(file,files, compression_level = 2)
        }else{
          NULL
        }
      }
    )
  }
)