# ----------------------------------------------------------------------------
# genome-wide sgRNA predictions
# ----------------------------------------------------------------------------
geneOfInterestGeneInputFile <- reactiveValues(data = NULL)


geneOfInterestUpdateText <- function(){
  output$geneOfInterestInfo <- renderText({
    if(is.null(input$geneOfInterestGeneSelect) & is.null(geneOfInterestGeneInputFile$data)){
      invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                       " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                       " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
    }else{
      "INFO: Click Load data!"
    }
  })
}
output$geneOfInterestInfo <- renderText({
  invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
})

#upon load display nothing
output$geneOfInterestScreenHitsTable <- renderDataTable({
})

output$geneOfInterestDEAHitsTable <- renderDataTable({
})

output$geneOfInterestScreenMetaTable <- renderDataTable({
})

output$geneOfInterestDEAMetaTable <- renderDataTable({
})

#Screen hits dataframe
geneOfInterestScreenHitsDataFrame <- reactive({
  
  if(!is.null(geneOfInterestGeneInputFile$data)){
    presel_genes_buff <- gwsGeneGeneList()
    genes_fileUpload <- c(paste0("\\(", (geneOfInterestGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (geneOfInterestGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes_both <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes_both<- input$geneOfInterestGeneSelect
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  #specify which table to select
  if(input$geneOfInterestSearchRadio == "guide_stats"){
    tableSelect <- "guide_stats"
  }else{
    tableSelect <- "gene_stats"
  }
  
  effect_threshold_depleting <- local(input$geneOfInterestScreenDepletingLFC)
  effect_threshold_enriching <- local(input$geneOfInterestScreenEnrichingLFC)
  
  con_geneOfInterest <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
  geneOfInterestDepleting <- data.frame()
  if("depletion" %in% input$geneOfInterestScreenInclude){
    geneOfInterestDepleting <- con_geneOfInterest %>%
      tbl(tableSelect) %>%
      dplyr::filter(gene_id %in% gene_select) %>%
      dplyr::filter_at(vars("lfc"), all_vars(. <= effect_threshold_depleting)) %>%
      dplyr::select(gene_id, matches("guide_id"), contrast_id, lfc, p, fdr) %>%
      collect()
  }
  geneOfInterestEnriching <- data.frame()
  if("enrichment" %in% input$geneOfInterestScreenInclude){
    geneOfInterestEnriching <- con_geneOfInterest %>%
      tbl(tableSelect) %>%
      dplyr::filter(gene_id %in% gene_select) %>%
      dplyr::filter_at(vars("lfc"), all_vars(. >= effect_threshold_enriching)) %>%
      dplyr::select(gene_id, matches("guide_id"), contrast_id, lfc, p, fdr) %>%
      collect()
  }
  DBI::dbDisconnect(con_geneOfInterest)
  
  out <- geneOfInterestDepleting %>%
    bind_rows(geneOfInterestEnriching)
  
  if(nrow(out) > 0){
    out %>%
      left_join(contrasts %>% dplyr::select(contrast_id, type, reference_type, library_id)) %>%
      left_join(features %>% dplyr::select(library_id, gene_id, symbol, entrez_id) %>% distinct) %>%
      dplyr::select(entrez_id, symbol, matches("guide_id"), contrast_id, library_id, type, reference_type, lfc, p, fdr)
  }
})

#DEA hits dataframe
geneOfInterestDEAHitsDataFrame <- reactive({
  
  if(!is.null(geneOfInterestGeneInputFile$data)){
    presel_genes_buff <- gwsGeneGeneList()
    genes_fileUpload <- c(paste0("\\(", (geneOfInterestGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (geneOfInterestGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes_both <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes_both<- input$geneOfInterestGeneSelect
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  effect_threshold_depleting <- local(input$geneOfInterestDEADepletingLFC)
  effect_threshold_enriching <- local(input$geneOfInterestDEAEnrichingLFC)
  
  #quantseq
  con_geneOfInterest <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/quantseq.db")
  
  geneOfInterestDepleting <- data.frame()
  if("depletion" %in% input$geneOfInterestDEAInclude){
    geneOfInterestDepleting <- con_geneOfInterest %>%
      tbl("quantseq_deseq2_results") %>%
      dplyr::filter(entrez_id %in% presel_entrez | symbol %in% presel_genes, padj <= 0.1) %>%
      dplyr::filter_at(vars("log2FoldChange"), all_vars(. <= effect_threshold_depleting)) %>%
      collect()
  }
  geneOfInterestEnriching <- data.frame()
  if("enrichment" %in% input$geneOfInterestDEAInclude){
    geneOfInterestEnriching <- con_geneOfInterest %>%
      tbl("quantseq_deseq2_results") %>%
      dplyr::filter(entrez_id %in% presel_entrez | symbol %in% presel_genes, padj <= 0.1) %>%
      dplyr::filter_at(vars("log2FoldChange"), all_vars(. >= effect_threshold_enriching)) %>%
      collect()
  }
  DBI::dbDisconnect(con_geneOfInterest)
  
  out_quantseq <- geneOfInterestDepleting %>%
    bind_rows(geneOfInterestEnriching) %>%
    mutate(dataset = "Quant-seq")
  
  #slamseq
  con_geneOfInterest <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")
  
  geneOfInterestDepleting <- data.frame()
  if("depletion" %in% input$geneOfInterestDEAInclude){
    geneOfInterestDepleting <- con_geneOfInterest %>%
      tbl("slamseq_deseq2_results") %>%
      dplyr::filter(entrez_id %in% presel_entrez | symbol %in% presel_genes, padj <= 0.1) %>%
      dplyr::filter_at(vars("log2FoldChange"), all_vars(. <= effect_threshold_depleting)) %>%
      collect()
  }
  geneOfInterestEnriching <- data.frame()
  if("enrichment" %in% input$geneOfInterestDEAInclude){
    geneOfInterestEnriching <- con_geneOfInterest %>%
      tbl("slamseq_deseq2_results") %>%
      dplyr::filter(entrez_id %in% presel_entrez | symbol %in% presel_genes, padj <= 0.1) %>%
      dplyr::filter_at(vars("log2FoldChange"), all_vars(. >= effect_threshold_enriching)) %>%
      collect()
  }
  DBI::dbDisconnect(con_geneOfInterest)
  
  out_slamseq <- geneOfInterestDepleting %>%
    bind_rows(geneOfInterestEnriching) %>%
    mutate(dataset = "SLAM-seq")
  
  #combine slamseq and quantseq
  out <- out_quantseq %>%
    bind_rows(out_slamseq) %>%
    dplyr::select(dataset, dplyr::everything())
  
  if(nrow(out) > 0){
    out
  }
})

#Screen meta dataframe
geneOfInterestScreenMetaDataFrame <- reactive({
  df <- geneOfInterestScreenHitsDataFrame()
  
  meta <- contrasts %>%
    dplyr::filter(contrast_id %in% df$contrast_id)
  
  if(nrow(meta) > 0){
    meta
  }
})

#DEA meta dataframe
geneOfInterestDEAMetaDataFrame <- reactive({
  df <- geneOfInterestDEAHitsDataFrame()
  
  meta <- deseq2_list_quantseq %>%
    bind_rows(deseq2_list_slamseq) %>%
    dplyr::filter(contrast_id %in% df$contrast_id)
  
  if(nrow(meta) > 0){
    meta
  }
})

#Screen hits datatable
geneOfInterestScreenHitsDataTable <- eventReactive(input$geneOfInterestLoadButton,{
  df <- geneOfInterestScreenHitsDataFrame()
  if(!is.null(df)) {
    output$geneOfInterestInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  fixedColumns = list(leftColumns = 2),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    output$geneOfInterestInfo <- renderText({
      "WARNING: No data found!"
    })
    NULL
  }
})

#DEA hits datatable
geneOfInterestDEAHitsDataTable <- eventReactive(input$geneOfInterestLoadButton,{
  df <- geneOfInterestDEAHitsDataFrame()
  
  if(!is.null(df)) {
    output$geneOfInterestInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  fixedColumns = list(leftColumns = 2),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    output$geneOfInterestInfo <- renderText({
      "WARNING: No data found!"
    })
    NULL
  }
})

#Screen meta datatable
geneOfInterestScreenMetaDataTable <- eventReactive(input$geneOfInterestLoadButton,{
  df <- geneOfInterestScreenMetaDataFrame()
  
  if(!is.null(df)) {
    output$geneOfInterestInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  fixedColumns = list(leftColumns = 1),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    output$geneOfInterestInfo <- renderText({
      "WARNING: No data found!"
    })
    NULL
  }
})

#DEA meta datatable
geneOfInterestDEAMetaDataTable <- eventReactive(input$geneOfInterestLoadButton,{
  df <- geneOfInterestDEAMetaDataFrame()
  
  if(!is.null(df)) {
    output$geneOfInterestInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    df %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  fixedColumns = list(leftColumns = 1),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    output$geneOfInterestInfo <- renderText({
      "WARNING: No data found!"
    })
    NULL
  }
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

geneOfInterestGeneList <- reactive({
  
  gene_list_slamseq %>%
    bind_rows(gene_list_quantseq) %>%
    dplyr::filter(species == input$geneOfInterestSpeciesSelect) %>%
    dplyr::select(Symbol=symbol, Entrez_id=entrez_id)
  
  if(input$geneOfInterestSpeciesSelect == "human"){
    gene_list <- gene_list_human 
  }else{
    gene_list <- gene_list_mouse 
  }
  
  gene_list %>%
    bind_rows(gene_list_slamseq) %>%
    distinct() %>%
    dplyr::mutate(gene = ifelse(is.na(Symbol), paste0("No symbol found (", EntrezID, ")"), paste0(Symbol , " (", EntrezID, ")"))) %>%
    arrange(gene) %>%
    .$gene
  
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------
observeEvent(input$geneOfInterestLoadButton, {
  
  #Screen hits
  output$geneOfInterestScreenHitsTable <- renderDataTable({
    dt <- geneOfInterestScreenHitsDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  #DEA hits
  output$geneOfInterestDEAHitsTable <- renderDataTable({
    dt <- geneOfInterestDEAHitsDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  #Screen meta
  output$geneOfInterestScreenMetaTable <- renderDataTable({
    dt <- geneOfInterestScreenMetaDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  
  #DEA meta
  output$geneOfInterestDEAMetaTable <- renderDataTable({
    dt <- geneOfInterestDEAMetaDataTable()
    if(!is.null(dt)){
      dt
    }
  })
})

observeEvent(input$geneOfInterestSpeciesSelect, {
  #update other species selects
  updateSpecies(input$geneOfInterestSpeciesSelect)
  #update gene selectbox
  updateSelectizeInput(session, 'geneOfInterestGeneSelect', choices = geneOfInterestGeneList(), server = TRUE)
  
})

observeEvent(input$geneOfInterestGeneSelect, {
  if((!is.null(input$geneOfInterestGeneSelect)) | !is.null(geneOfInterestGeneInputFile$data)){
    enable("geneOfInterestLoadButton")
    if(!is.null(input$geneOfInterestGeneSelect)){
      geneOfInterestGeneInputFile$data <- NULL
      reset('geneOfInterestGeneInputFile')
    }
  }else{
    disable("geneOfInterestLoadButton")
  }
  geneOfInterestUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$geneOfInterestGeneInputFile, {
  if(!is.null(input$geneOfInterestGeneInputFile)){
    updateSelectizeInput(session, 'geneOfInterestGeneSelect', choices = geneOfInterestGeneList(), server = TRUE)
    req(input$geneOfInterestGeneInputFile)
    geneOfInterestGeneInputFile$data <- read_tsv(input$geneOfInterestGeneInputFile$datapath, col_names = F)
  }else{
    geneOfInterestGeneInputFile$data <- NULL
  }
  if((!is.null(input$geneOfInterestGeneSelect)) | !is.null(geneOfInterestGeneInputFile$data)){
    enable("geneOfInterestLoadButton")
  }else{
    disable("geneOfInterestLoadButton")
  }
  geneOfInterestUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$geneOfInterestScreenHitsButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_geneOfInterest_screenHits.txt"
  },
  content = function(file) {
    geneOfInterestScreenHitsDataFrame() %>% write_tsv(file)
  }
)

output$geneOfInterestDEAHitsButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_geneOfInterest_DEAHits.txt"
  },
  content = function(file) {
    geneOfInterestDEAHitsDataFrame() %>% write_tsv(file)
  }
)

output$geneOfInterestScreenMetaButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_geneOfInterest_screenMetadata.txt"
  },
  content = function(file) {
    geneOfInterestScreenMetaDataFrame() %>% write_tsv(file)
  }
)

output$geneOfInterestDEAMetaButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_geneOfInterest_DEAMetadata.txt"
  },
  content = function(file) {
    geneOfInterestDEAMetaDataFrame() %>% write_tsv(file)
  }
)