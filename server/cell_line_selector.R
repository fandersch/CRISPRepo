# ----------------------------------------------------------------------------
# cellLineSelector
# ----------------------------------------------------------------------------
cellLineSelector_geneInputFile <- reactiveValues(data = NULL)


cellLineSelectorUpdateText <- function(){
  output$cellLineSelectorInfo <- renderText({
    if(is.null(input$cellLineSelectorTissueSelect) & !isTRUE(input$cellLineSelectorCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      invisible("INFO: Please select further cell line attributes in the right panel! Click the 'Load data!'-button if you are satisfied with your selection.")
    }
  })
}

#upon load display nothing
output$cellLineSelectorDataTableMeta <- renderDataTable({
})
output$cellLineSelectorDataTableScreens <- renderDataTable({
})
output$cellLineSelectorDataTableExpression <- renderDataTable({
})
output$cellLineSelectorDataTableMutations <- renderDataTable({
})
output$cellLineSelectorDataTableFusions <- renderDataTable({
})
output$cellLineSelectorDataTableCNVs <- renderDataTable({
})
output$cellLineSelectorDataTableHLAs <- renderDataTable({
})

#--------------------------------------------------------
# query database and create dataframe
#--------------------------------------------------------

#
cellLineSelectorMetaDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLineSelector_meta_data <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(sanger_model_ID %in% local(c(model_ids))) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)

  cellLineSelector_meta_data
})

#
cellLineSelectorScreenDataFrame <- reactive({
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
  model_ids <- cellLineSelectorModelList()
  
  cellLineSelector_meta_data_cellline <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(sanger_model_ID %in% local(c(model_ids))) %>%
    collect() %>%
    .$cell_line_name
  
  contrast_ids <- contrasts %>%
    dplyr::filter(type == "dropout", dynamic_range <= -1.5, auc > 0.9, cellline_name %in% local(cellLineSelector_meta_data_cellline)) %>%
    distinct() %>%
    .$contrast_id
  
  library_ids <- contrasts %>%
    dplyr::filter(type == "dropout", dynamic_range <= -1.5, auc > 0.9, cellline_name %in% local(cellLineSelector_meta_data_cellline)) %>%
    distinct() %>%
    .$library_id
  
  genes <- features %>%
    filter(library_id %in% local(library_ids)) %>%
    dplyr::select(gene_id, symbol, entrez_id) %>%
    distinct()
  
  if(length(contrast_ids)>=900){
    contrasts_ids_filter_string <- c(paste(paste("contrast_id", paste0("'", contrast_ids[1:899], "'"), sep="="), collapse=" OR "))
    i<-900
    while(i <= length(contrast_ids)){
      end<-i+899
      if(end>length(contrast_ids)){
        end<-length(contrast_ids)
      }
      contrasts_ids_filter_string <- c(sample_ids_filter_string, paste(paste("contrast_id", paste0("'", contrast_ids[i:end], "'"), sep="="), collapse=" OR "))
      i<-i+end
    }
  }else{
    contrasts_ids_filter_string <- paste(paste("contrast_id", paste0("'", contrast_ids, "'"), sep="="), collapse=" OR ")
  }
  
  query <- paste0("SELECT contrast_id, gene_id, adjusted_effect_essentialome FROM gene_stats ",
                  "WHERE (", contrasts_ids_filter_string, ") ")
  
  
  df<- NULL
  for(z in 1:length(query)){
    # a chunk at a time
    res <- dbSendQuery(con, query[z])
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
    chunk<-NULL
    dbClearResult(res)
  }
  
  DBI::dbDisconnect(con)
  DBI::dbDisconnect(con_cell_lines)

  if(!is.null(df)){
    df %>% 
      dplyr::distinct() %>%
      inner_join(genes) %>%
      dplyr::select(-gene_id) %>%
      pivot_wider(names_from = "contrast_id", values_from = "adjusted_effect_essentialome") %>%
      arrange(symbol)
  }else{
    data.frame()
  }
})

#
cellLineSelectorExpressionDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm.db")
  
  cellLineSelector_meta_data_cellline <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(sanger_model_ID %in% local(c(model_ids))) %>%
    collect() %>%
    .$cell_line_name
  
  sample_ids <- cellline_list_expressionData %>%
    dplyr::filter(cell_line_name %in% cellLineSelector_meta_data_cellline) %>%
    dplyr::select(sample_id) %>%
    .$sample_id
  
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
  
  query <- paste0("SELECT sample_id, symbol, entrez_id, tpm FROM expression_data_values ",
                  "WHERE (", sample_ids_filter_string, ") ")
  
  df<- NULL
  for(z in 1:length(query)){
    # a chunk at a time
    res <- dbSendQuery(con_expression, query[z])
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
    chunk<-NULL
    dbClearResult(res)
  }
  
  DBI::dbDisconnect(con_cell_lines)
  DBI::dbDisconnect(con_expression)
  
  if(!is.null(df)){
    df %>% 
      distinct %>%
      dplyr::select(sample_id, symbol, entrez_id, tpm) %>%
      pivot_wider(names_from=sample_id, values_from=tpm) %>%
      arrange(symbol)
  }else{
    data.frame()
  }
  
})

#
cellLineSelectorMutationsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLineSelector_mutation_data <- con_cell_lines %>%
    tbl("cell_line_gene_mutations") %>%
    dplyr::filter(model_id %in% local(model_ids))
  
  if(!is.null(input$cellLineSelectorGeneMutationSelect)){
    cellLineSelector_mutation_data <- cellLineSelector_mutation_data %>%
      dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect))
  }
  if(!is.null(input$cellLineSelectorGeneMutationProteinSelect)){
    cellLineSelector_mutation_data <- cellLineSelector_mutation_data %>%
      dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect), protein_mutation %in% local(input$cellLineSelectorGeneMutationProteinSelect)) 
  }
  
  cellLineSelector_mutation_data <- cellLineSelector_mutation_data %>%
    dplyr::distinct() %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)

  cellLineSelector_mutation_data
})

#
cellLineSelectorFusionsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLineSelector_fusion_data <- con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% local(c(model_ids)))
  
  if(!is.null(input$cellLineSelectorGeneFusion3primeSelect)){
    cellLineSelector_fusion_data <- cellLineSelector_fusion_data %>%
      dplyr::filter(gene_symbol_3prime %in% local(input$cellLineSelectorGeneFusion3primeSelect))
  }
    
  if(!is.null(input$cellLineSelectorGeneFusion5primeSelect)){
    cellLineSelector_fusion_data <- cellLineSelector_fusion_data %>%
      dplyr::filter(gene_symbol_5prime %in% local(input$cellLineSelectorGeneFusion5primeSelect))
  }
  
  cellLineSelector_fusion_data <- cellLineSelector_fusion_data %>%
    dplyr::distinct() %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLineSelector_fusion_data
    
})

#
cellLineSelectorCNVsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLineSelector_cnv_data <- con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% local(c(model_ids)))
  
  if(!is.null(input$cellLineSelectorGeneCNVSelect)){
    cellLineSelector_cnv_data <- cellLineSelector_cnv_data %>%
      dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneCNVSelect))
  }
  if(!is.null(input$cellLineSelectorGeneCNVCategorySelect)){
    cellLineSelector_cnv_data <- cellLineSelector_cnv_data %>%
      dplyr::filter(cn_category %in% local(input$cellLineSelectorGeneCNVCategorySelect))
  }
  
  cellLineSelector_cnv_data <- cellLineSelector_cnv_data %>%
    dplyr::distinct() %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLineSelector_cnv_data
  
})

#
cellLineSelectorHLAsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLineSelector_meta_data_cellline <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(sanger_model_ID %in% local(c(model_ids))) %>%
    collect() %>%
    .$cell_line_name
  
  cellLineSelector_HLA_data <- con_cell_lines %>%
    tbl("cell_line_hla_type") %>%
    dplyr::filter(cell_line_name %in% local(c(cellLineSelector_meta_data_cellline))) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLineSelector_HLA_data
  
})

#---------------------------------------------
# create datatable out of dataframe
#---------------------------------------------
cellLineSelectorDataTableMeta <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorMetaDataFrame()

  if (nrow(df) > 0) {

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

cellLineSelectorDataTableScreens <- eventReactive(input$cellLineSelectorLoadButton,{
  
  df <- cellLineSelectorScreenDataFrame()
  
  if (nrow(df) > 0) {
    
    vals<-df[,3:ncol(df)]
    values_min <- vals %>% min(na.rm = TRUE)
    values_max <- vals %>% max(na.rm = TRUE)
    
    brks_smaller <- seq(values_min, 0, .05)
    brks_bigger <- seq(0, values_max, .05)
    
    clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
    clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
      {paste0("rgb(", ., ",", ., ",255)")}
    
    brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
    clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
    
    df %>%
      DT::datatable(extensions = c('FixedColumns','FixedHeader'), 
                    options = list(
                      autoWidth = FALSE,
                      headerCallback = JS(headerCallback),
                      scrollX=TRUE,
                      fixedColumns = list(leftColumns = 2),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
                      searchHighlight = TRUE
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
      formatStyle(seq(3, ncol(df),1),
                  backgroundColor = styleInterval(brks, clrs))
  }else{
    NULL
  }
})

cellLineSelectorDataTableExpression <- eventReactive(input$cellLineSelectorLoadButton,{
  
  df <- cellLineSelectorExpressionDataFrame()
  
  if (nrow(df) > 0) {
    
    nfreezeColumns <- 2
    
    values<- df[,3:ncol(df)]
    max_value <- max(values, na.rm=T)
    min_value <- min(values, na.rm=T)
    
    brks <- exp(seq(0, log(max_value), length.out = 40))
    
    clrs <- round(seq(255, 5, length.out = (length(brks) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
    
    df %>%
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
      formatStyle(seq(nfreezeColumns+1, length(colnames(df)),1),
                  backgroundColor = styleInterval(brks, clrs))
    
  }else{
    NULL
  }
})

cellLineSelectorDataTableMutations <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorMutationsDataFrame()

  if (nrow(df) > 0) {

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

cellLineSelectorDataTableFusions <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorFusionsDataFrame()

  if (nrow(df) > 0) {

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

cellLineSelectorDataTableCNVs <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorCNVsDataFrame()

  if (nrow(df) > 0) {
    
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

cellLineSelectorDataTableHLAs <- eventReactive(input$cellLineSelectorLoadButton,{
  
  df <- cellLineSelectorHLAsDataFrame()
  
  if (nrow(df) > 0) {
    
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
# Reactive cell line filtering
#----------------------------------------------------------------------------

cellLineSelectorModelList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm.db")
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$sanger_model_ID
  
  #get cell lines with requested mutations
  cell_lines_mutation<-c()
  if(!is.null(input$cellLineSelectorGeneMutationSelect) & length(model_ids) > 0){
    if(is.null(input$cellLineSelectorGeneMutationProteinSelect)){
      cell_lines_mutation <- con_cell_lines %>%
        tbl("cell_line_gene_mutations") %>%
        dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
      
    }else{
      
      cell_lines_mutation <- con_cell_lines %>%
        tbl("cell_line_gene_mutations") %>%
        dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect), protein_mutation %in% local(input$cellLineSelectorGeneMutationProteinSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
    }
    model_ids <- model_ids[model_ids %in% cell_lines_mutation]
  }
  
  #get cell lines with requested gene fusions
  cell_lines_fusion<-c()
  if(!is.null(input$cellLineSelectorGeneFusion3primeSelect) & length(model_ids) > 0){
    
    if(!is.null(input$cellLineSelectorGeneFusion5primeSelect)){
      
      cell_lines_fusion <- con_cell_lines %>%
        tbl("cell_line_gene_fusions") %>%
        dplyr::filter(gene_symbol_3prime %in% local(input$cellLineSelectorGeneFusion3primeSelect), gene_symbol_5prime %in% local(input$cellLineSelectorGeneFusion5primeSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
      
    }else{
      
      cell_lines_fusion <- con_cell_lines %>%
        tbl("cell_line_gene_fusions") %>%
        dplyr::filter(gene_symbol_3prime %in% local(input$cellLineSelectorGeneFusion3primeSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
    }
    model_ids <- model_ids[model_ids %in% cell_lines_fusion]
  }else{
    if(!is.null(input$cellLineSelectorGeneFusion5primeSelect) & length(model_ids) > 0){
      
      cell_lines_fusion <- con_cell_lines %>%
        tbl("cell_line_gene_fusions") %>%
        dplyr::filter(gene_symbol_5prime %in% local(input$cellLineSelectorGeneFusion5primeSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
      
      model_ids <- model_ids[model_ids %in% cell_lines_fusion]
    }
  }
  
  #get cell lines with requested gene CNV
  cell_lines_cnv<-c()
  if(!is.null(input$cellLineSelectorGeneCNVSelect) & length(model_ids) > 0){
    if(is.null(input$cellLineSelectorGeneCNVCategorySelect)){
      
      cell_lines_cnv <- con_cell_lines %>%
        tbl("cell_line_gene_cnv") %>%
        dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneCNVSelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
      
    }else{
      
      cell_lines_cnv <- con_cell_lines %>%
        tbl("cell_line_gene_cnv") %>%
        dplyr::filter(gene_symbol %in% local(input$cellLineSelectorGeneCNVSelect), cn_category %in% local(input$cellLineSelectorGeneCNVCategorySelect)) %>%
        dplyr::select(model_id) %>%
        dplyr::distinct() %>%
        collect() %>%
        .$model_id
    }
    
    model_ids <- model_ids[model_ids %in% cell_lines_cnv]

  }
  
  #get cell lines with requested gene HLA type
  cell_lines_hla<-c()
  if(local(input$cellLineSelectorHLAtypeSelect) != "none" & length(model_ids) > 0){
    if(is.null(input$cellLineSelectorHLAalleleSelect)){
      
      type <- local(input$cellLineSelectorHLAtypeSelect) %>%
        str_replace("HLA_", "RPKM_")
      
      cellline_names <- con_cell_lines %>%
        tbl("cell_line_hla_type") %>%
        collect() %>%
        dplyr::filter(if_all(starts_with(type), ~ . >= local(input$cellLineSelectorHLAExpressionSlider))) %>%
        dplyr::select(cell_line_name) %>%
        dplyr::distinct() %>%
        .$cell_line_name
      
    }else{
      type <- local(input$cellLineSelectorHLAtypeSelect) %>%
        str_replace("HLA_", "RPKM_")

      allele <- c(paste0(local(input$cellLineSelectorHLAtypeSelect), "_A1"),
                   paste0(local(input$cellLineSelectorHLAtypeSelect), "_A2"))
      
      cellline_names <- con_cell_lines %>%
        tbl("cell_line_hla_type") %>%
        collect() %>%
        dplyr::filter(if_any(starts_with(allele), ~ . %in% local(input$cellLineSelectorHLAalleleSelect))) %>%
        dplyr::filter(if_all(starts_with(type), ~ . >= local(input$cellLineSelectorHLAExpressionSlider))) %>%
        dplyr::select(cell_line_name) %>%
        dplyr::distinct() %>%
        .$cell_line_name
    }
    
    cell_lines_hla <- con_cell_lines %>%
      tbl("cell_line_meta") %>%
      dplyr::filter(cell_line_name %in% local(cellline_names)) %>%
      dplyr::select(sanger_model_ID) %>%
      dplyr::distinct() %>%
      collect() %>%
      .$sanger_model_ID
    
    model_ids <- model_ids[model_ids %in% cell_lines_hla]
    
  }
  
  #get cell lines with requested gene dependency
  cell_lines_dependency<-c()
  if(!is.null(input$cellLineSelectorGeneDependencySelect) & length(model_ids) > 0){
    
    gene_ids <- features %>%
      dplyr::filter(symbol %in% local(input$cellLineSelectorGeneDependencySelect)) %>%
      select(gene_id) %>%
      distinct() %>%
      .$gene_id
    
    contrasts_filtered <- contrasts %>%
      dplyr::filter(type == "dropout", dynamic_range <= -1.5, auc > 0.9) %>%
      distinct()
      
    contrast_ids <- con %>%
      tbl("gene_stats") %>%
      dplyr::filter(contrast_id %in% local(contrasts_filtered$contrast_id), gene_id %in% local(gene_ids), adjusted_effect_essentialome <= local(input$cellLineSelectorGeneDependencySlider)) %>%
      dplyr::select(contrast_id, adjusted_effect_essentialome) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(contrast_id) %>%
      summarize(n_check=sum(adjusted_effect_essentialome <= local(input$cellLineSelectorGeneDependencySlider)), n_genes=length(local(input$cellLineSelectorGeneDependencySelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$contrast_id %>%
      unique
    
    cellline_names <- contrasts_filtered %>%
      dplyr::filter(contrast_id %in% contrast_ids) %>%
      .$cellline_name %>%
      unique
    
    cell_lines_dependency <- con_cell_lines %>%
      tbl("cell_line_meta") %>%
      dplyr::filter(cell_line_name %in% local(cellline_names)) %>%
      dplyr::select(sanger_model_ID) %>%
      dplyr::distinct() %>%
      collect() %>%
      .$sanger_model_ID
    
    model_ids <- model_ids[model_ids %in% cell_lines_dependency]
  }
  
  #get cell lines with requested gene expression
  cell_lines_expression<-c()
  if(!is.null(input$cellLineSelectorGeneExpressionSelect) & length(model_ids) > 0){
    
    sample_ids_filtered <- cellline_list_expressionData %>%
      dplyr::filter(tissue_name %in% presel_tissue) %>%
      dplyr::select(sample_id) %>%
      distinct() %>%
      collect() %>%
      .$sample_id
    
    sample_ids <- con_expression %>%
      tbl("expression_data_values") %>%
      dplyr::filter(sample_id %in% local(sample_ids_filtered), symbol %in% local(input$cellLineSelectorGeneExpressionSelect), tpm >= local(input$cellLineSelectorGeneExpressionSlider)) %>%
      dplyr::select(sample_id, tpm) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(sample_id) %>%
      summarize(n_check=sum(tpm >= local(input$cellLineSelectorGeneExpressionSlider)), n_genes=length(local(input$cellLineSelectorGeneExpressionSelect))) %>%
      dplyr::filter(n_check==n_genes) %>%
      .$sample_id %>%
      unique
    
    
    cellline_names <- cellline_list_expressionData %>%
      dplyr::filter(sample_id %in% sample_ids) %>%
      .$cell_line_name %>%
      unique
    
    cell_lines_expression <- con_cell_lines %>%
      tbl("cell_line_meta") %>%
      dplyr::filter(cell_line_name %in% local(cellline_names)) %>%
      dplyr::select(sanger_model_ID) %>%
      dplyr::distinct() %>%
      collect() %>%
      .$sanger_model_ID
    
    model_ids <- model_ids[model_ids %in% cell_lines_expression]
  }
  
  #close connections
  DBI::dbDisconnect(con_cell_lines)
  DBI::dbDisconnect(con)
  DBI::dbDisconnect(con_expression)
  
  model_ids %>% 
    unique
})


#----------------------------------------------------------------------------
# Selectbox Lists
#----------------------------------------------------------------------------

cellLineSelectorGeneMutationProteinList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    collect() %>%
    .$sanger_model_ID
  
  protein_mutation <- con_cell_lines %>%
    tbl("cell_line_gene_mutations") %>%
    dplyr::filter(model_id %in% local(model_ids), gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect)) %>%
    dplyr::select(protein_mutation) %>%
    distinct() %>%
    collect() %>%
    arrange(protein_mutation) %>%
    .$protein_mutation
  
  DBI::dbDisconnect(con_cell_lines)
  
  protein_mutation
    
})

cellLineSelectorGeneDependencyList <- reactive({
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  model_names <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(cell_line_name) %>%
    collect() %>%
    .$cell_line_name
  
  DBI::dbDisconnect(con_cell_lines)
  
  library_ids <- libraries %>%
    dplyr::filter(tissue_name %in% presel_tissue, cellline_name %in% local(model_names),  species == "human") %>%
    dplyr::select(library_id) %>%
    distinct %>%
    .$library_id
  
  gene_list_screens %>%
    dplyr::filter(library_id %in% local(library_ids)) %>%
    dplyr::select(symbol) %>%
    distinct() %>%
    arrange(symbol) %>%
    .$symbol
  
})

cellLineSelectorGeneExpressionList <- reactive({

  gene_list_expressionData %>%
    filter(species == "human") %>%
    arrange(symbol) %>%
    .$symbol
  
})

cellLineSelectorGeneFusion3primeList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")

  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$sanger_model_ID
  
  gene_symbol_3prime  <- con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol_3prime) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol_3prime) %>%
    .$gene_symbol_3prime
  
  DBI::dbDisconnect(con_cell_lines)
  
  gene_symbol_3prime

})

cellLineSelectorGeneFusion5primeList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$sanger_model_ID
  
  gene_symbol_5prime <- con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol_5prime) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol_5prime) %>%
    .$gene_symbol_5prime
  
  DBI::dbDisconnect(con_cell_lines)
  
  gene_symbol_5prime
    
})

cellLineSelectorGeneCNVList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$sanger_model_ID
  
  gene_symbol <- con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol) %>%
    .$gene_symbol
  
  DBI::dbDisconnect(con_cell_lines)
  
  gene_symbol
})

cellLineSelectorGeneCNVCategoryList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")

  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(sanger_model_ID) %>%
    collect() %>%
    .$sanger_model_ID
  
  cn_category <- con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(cn_category) %>%
    distinct() %>%
    collect() %>%
    arrange(cn_category) %>%
    .$cn_category
  
  DBI::dbDisconnect(con_cell_lines)
  
  cn_category
})

cellLineSelectorHLAalleleList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  hla <- NULL
  
  if(local(input$cellLineSelectorHLAtypeSelect) != "none"){
    con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
    
    cell_line_names <- con_cell_lines %>%
      tbl("cell_line_meta") %>%
      dplyr::filter(tissue_name %in% presel_tissue) %>%
      dplyr::select(cell_line_name) %>%
      distinct() %>%
      collect() %>%
      .$cell_line_name
    
    hla <- con_cell_lines %>%
      tbl("cell_line_hla_type") %>%
      dplyr::filter(cell_line_name %in% local(cell_line_names)) %>%
      dplyr::select(matches(local(input$cellLineSelectorHLAtypeSelect))) %>%
      distinct() %>%
      collect() %>%
      unlist() %>%
      as.character %>%
      sort
    
    DBI::dbDisconnect(con_cell_lines)
  }

  
  hla[hla!="no"]
})


#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
observeEvent(input$cellLineSelectorLoadButton, {
  
  dt_meta <- cellLineSelectorDataTableMeta()

  output$cellLineSelectorDataTableMeta <- renderDataTable({
    if(!is.null(dt_meta)){
      dt_meta
    }
  })
  
  dt_screens <- cellLineSelectorDataTableScreens()
  
  output$cellLineSelectorDataTableScreens <- renderDataTable({
    
    if(!is.null(dt_screens)){
      dt_screens
    }
  })
  
  dt_expression <- cellLineSelectorDataTableExpression()
  
  output$cellLineSelectorDataTableExpression <- renderDataTable({
    if(!is.null(dt_expression)){
      dt_expression
    }
  })
  
  dt_mutation <- cellLineSelectorDataTableMutations()
  
  output$cellLineSelectorDataTableMutations <- renderDataTable({
    if(!is.null(dt_mutation)){
      dt_mutation
    }
  })
  
  dt_fusion <- cellLineSelectorDataTableFusions()

  output$cellLineSelectorDataTableFusions <- renderDataTable({
    if(!is.null(dt_fusion)){
      dt_fusion
    }
  })
  
  dt_cnv <- cellLineSelectorDataTableCNVs()

  output$cellLineSelectorDataTableCNVs <- renderDataTable({
    if(!is.null(dt_cnv)){
      dt_cnv
    }
  })
  
  dt_hla <- cellLineSelectorDataTableHLAs()
  
  output$cellLineSelectorDataTableHLAs <- renderDataTable({
    if(!is.null(dt_hla)){
      dt_hla
    }
  })
  
  if(is.null(dt_meta) & is.null(dt_screens) & is.null(dt_expression) & is.null(dt_mutation) & is.null(dt_fusion) & is.null(dt_cnv) & is.null(dt_hla)){
    output$cellLineSelectorInfo <- renderText({
      "ATTENTION: No cell lines found that fulfill filtering criteria. Please adjust the query!"
    })
  }else{
    output$cellLineSelectorInfo <- renderText({
      "SUCCESS: Loading completed!"
    })
  }
  
})

observeEvent(input$cellLineSelectorTissueSelect, {
  if(!is.null(input$cellLineSelectorTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'cellLineSelectorCheckTissueAll', value = FALSE)
  }
  
  if((isTRUE(input$cellLineSelectorCheckTissueAll) | (!is.null(input$cellLineSelectorTissueSelect)))) {
    enable("cellLineSelectorLoadButton")
    enable("cellLineSelectorGeneMutationSelect")
    enable("cellLineSelectorGeneDependencySelect")
    enable("cellLineSelectorGeneExpressionSelect")
    enable("cellLineSelectorGeneFusion3primeSelect")
    enable("cellLineSelectorGeneFusion5primeSelect")
    enable("cellLineSelectorGeneCNVSelect")
    enable("cellLineSelectorGeneCNVCategorySelect")
    enable("cellLineSelectorHLAtypeSelect")
    enable("cellLineSelectorHLAalleleSelect")
    
  }else{
    disable("cellLineSelectorLoadButton")
    disable("cellLineSelectorGeneMutationSelect")
    disable("cellLineSelectorGeneDependencySelect")
    disable("cellLineSelectorGeneExpressionSelect")
    disable("cellLineSelectorGeneFusion3primeSelect")
    disable("cellLineSelectorGeneFusion5primeSelect")
    disable("cellLineSelectorGeneCNVSelect")
    disable("cellLineSelectorGeneCNVCategorySelect")
    disable("cellLineSelectorHLAtypeSelect")
    disable("cellLineSelectorHLAalleleSelect")
  }

  #update selectb
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationSelect', choices = gene_list_cellLine$gene_symbol, selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneDependencySelect', choices = cellLineSelectorGeneDependencyList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneExpressionSelect', choices = cellLineSelectorGeneExpressionList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion3primeSelect', choices = cellLineSelectorGeneFusion3primeList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion5primeSelect', choices = cellLineSelectorGeneFusion5primeList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVSelect', choices = cellLineSelectorGeneCNVList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVCategorySelect', choices = cellLineSelectorGeneCNVCategoryList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorHLAalleleSelect', choices = cellLineSelectorHLAalleleList(), selected=NULL, server = TRUE)

  cellLineSelectorUpdateText()

}, ignoreNULL = FALSE)

observeEvent(input$cellLineSelectorCheckTissueAll, {
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'cellLineSelectorTissueSelect', choices = tissue_list_cellLine, selected=NULL, server = TRUE)
  }
  
  if((isTRUE(input$cellLineSelectorCheckTissueAll) | (!is.null(input$cellLineSelectorTissueSelect)))) {
    enable("cellLineSelectorLoadButton")
    enable("cellLineSelectorGeneMutationSelect")
    enable("cellLineSelectorGeneDependencySelect")
    enable("cellLineSelectorGeneExpressionSelect")
    enable("cellLineSelectorGeneFusion3primeSelect")
    enable("cellLineSelectorGeneFusion5primeSelect")
    enable("cellLineSelectorGeneCNVSelect")
    enable("cellLineSelectorGeneCNVCategorySelect")
    enable("cellLineSelectorHLAtypeSelect")
    enable("cellLineSelectorHLAalleleSelect")
  }else{
    disable("cellLineSelectorLoadButton")
    disable("cellLineSelectorGeneMutationSelect")
    disable("cellLineSelectorGeneDependencySelect")
    disable("cellLineSelectorGeneExpressionSelect")
    disable("cellLineSelectorGeneFusion3primeSelect")
    disable("cellLineSelectorGeneFusion5primeSelect")
    disable("cellLineSelectorGeneCNVSelect")
    disable("cellLineSelectorGeneCNVCategorySelect")
    disable("cellLineSelectorHLAtypeSelect")
    disable("cellLineSelectorHLAalleleSelect")
  }

  #update selectb
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationSelect', choices = gene_list_cellLine$gene_symbol, selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneDependencySelect', choices = cellLineSelectorGeneDependencyList(),  selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneExpressionSelect', choices = cellLineSelectorGeneExpressionList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion3primeSelect', choices = cellLineSelectorGeneFusion3primeList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion5primeSelect', choices = cellLineSelectorGeneFusion5primeList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVSelect', choices = cellLineSelectorGeneCNVList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVCategorySelect', choices = cellLineSelectorGeneCNVCategoryList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorHLAalleleSelect', choices = cellLineSelectorHLAalleleList(), selected=NULL, server = TRUE)

  cellLineSelectorUpdateText()

}, ignoreNULL = FALSE)


observeEvent(input$cellLineSelectorGeneMutationSelect, {

  #update gene selectb
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationProteinSelect', choices = cellLineSelectorGeneMutationProteinList(), selected=NULL, server = TRUE)

  if((!is.null(input$cellLineSelectorGeneMutationSelect))) {
    enable("cellLineSelectorGeneMutationProteinSelect")
  }else{
    disable("cellLineSelectorGeneMutationProteinSelect")
  }
  cellLineSelectorUpdateText()

}, ignoreNULL = FALSE)

observeEvent(input$cellLineSelectorHLAtypeSelect, {
  
  #update hla selectb
  updateSelectizeInput(session, 'cellLineSelectorHLAalleleSelect', choices = cellLineSelectorHLAalleleList(), selected=NULL, server = TRUE)
  
  cellLineSelectorUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$cellLineSelectorButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_cellLineSelector_meta_data.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading data"),
      value = 0,
      {
        meta <- cellLineSelectorMetaDataFrame()
        screen <- cellLineSelectorScreenDataFrame()
        expression <- cellLineSelectorExpressionDataFrame()
        mutations <- cellLineSelectorMutationsDataFrame()
        fusions <- cellLineSelectorFusionsDataFrame()
        cnvs <- cellLineSelectorCNVsDataFrame()
        hlas <- cellLineSelectorHLAsDataFrame()

        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        if("Cellline meta data" %in% input$cellLineSelectorDownloadCheck & nrow(meta)>0){
          #write each sheet to a csv file, save the name
          table <- meta
          fileName <- "cellLine_meta_data.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(1/9)
        
        if("Screen results" %in% input$cellLineSelectorDownloadCheck & nrow(meta)>0){
          #write each sheet to a csv file, save the name
          table <- screen
          fileName <- "screen_data.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(2/9)
        
        if("Expression levels" %in% input$cellLineSelectorDownloadCheck & nrow(meta)>0){
          #write each sheet to a csv file, save the name
          table <- expression
          fileName <- "expression_data.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(3/9)

        if("Gene mutations" %in% input$cellLineSelectorDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          table <- mutations
          fileName <- "gene_mutations.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(4/9)

        if("Gene fusions" %in% input$cellLineSelectorDownloadCheck & nrow(fusions)>0){
          #write each sheet to a csv file, save the name
          table <- fusions
          fileName <- "gene_fusions.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(5/9)

        if("Gene CNVs" %in% input$cellLineSelectorDownloadCheck & nrow(cnvs)>0){
          #write each sheet to a csv file, save the name
          table <- cnvs
          fileName <- "gene_cnvs.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(6/9)
        
        if("HLA types" %in% input$cellLineSelectorDownloadCheck & nrow(cnvs)>0){
          #write each sheet to a csv file, save the name
          table <- hlas
          fileName <- "hla_types.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(7/9)

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