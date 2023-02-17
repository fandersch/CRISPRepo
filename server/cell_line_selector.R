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
output$cellLineSelectorDataTableMutations <- renderDataTable({
})
output$cellLineSelectorDataTableFusions <- renderDataTable({
})
output$cellLineSelectorDataTableCNVs <- renderDataTable({
})

#--------------------------------------------------------
# query database and create dataframe
#--------------------------------------------------------

#
cellLineSelectorMetaDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  cellLineSelector_meta_data <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(Sanger_model_ID %in% local(c(model_ids))) %>%
    collect()

  cellLineSelector_meta_data
})

#
cellLineSelectorScreenDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
  cellLineSelector_meta_data_cellline <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(Sanger_model_ID %in% local(c(model_ids))) %>%
    collect() %>%
    .$cell_line_name

  contrasts_filtered <- contrasts %>%
    dplyr::filter(type == "dropout", dynamic_range <= -1.5, auc > 0.9, cellline_name %in% local(cellLineSelector_meta_data_cellline)) %>%
    distinct() %>%
    collect
  
  genes <- features %>%
    filter(library_id %in% local(contrasts_filtered$library_id)) %>%
    dplyr::select(gene_id, symbol, entrez_id) %>%
    distinct() %>%
    collect()
  
  gene_stats <- con %>%
    tbl("gene_stats") %>%
    dplyr::filter(contrast_id %in% local(contrasts_filtered$contrast_id)) %>%
    dplyr::select(contrast_id, gene_id, adjusted_effect) %>%
    dplyr::distinct() %>%
    collect() %>%
    dplyr::mutate(adjusted_effect = round(adjusted_effect,2)) %>%
    inner_join(contrasts_filtered %>% dplyr::select(contrast_id, contrast_id_QC)) %>%
    inner_join(genes) %>%
    dplyr::select(-contrast_id, -gene_id) %>%
    pivot_wider(names_from = "contrast_id_QC", values_from = "adjusted_effect") %>%
    dplyr::select(symbol, entrez_id, matches(contrasts_filtered$contrast_id_QC)) %>%
    arrange(symbol)
  
})

#
cellLineSelectorMutationsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
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
  
  cellLineSelector_mutation_data %>%
    dplyr::distinct() %>%
    collect()
})

#
cellLineSelectorFusionsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
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
  
  cellLineSelector_fusion_data %>%
    dplyr::distinct() %>%
    collect()
})

#
cellLineSelectorCNVsDataFrame <- reactive({
  model_ids <- cellLineSelectorModelList()
  
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

  cellLineSelector_cnv_data %>%
    dplyr::distinct() %>%
    collect()
})

#---------------------------------------------
# create datatable out of dataframe
#---------------------------------------------
cellLineSelectorDataTableMeta <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorMetaDataFrame()

  if (nrow(df) > 0) {
    output$cellLineSelectorInfo <- renderText({
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

cellLineSelectorDataTableScreens <- eventReactive(input$cellLineSelectorLoadButton,{
  
  df <- cellLineSelectorScreenDataFrame()
  
  if (nrow(df) > 0) {
    output$cellLineSelectorInfo <- renderText({
      "INFO: Loading completed!"
    })
    
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
    
    dt <- df %>%
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

cellLineSelectorDataTableMutations <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorMutationsDataFrame()

  if (nrow(df) > 0) {
    output$cellLineSelectorInfo <- renderText({
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

cellLineSelectorDataTableFusions <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorFusionsDataFrame()

  if (nrow(df) > 0) {
    output$cellLineSelectorInfo <- renderText({
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

cellLineSelectorDataTableCNVs <- eventReactive(input$cellLineSelectorLoadButton,{

  df <- cellLineSelectorCNVsDataFrame()

  if (nrow(df) > 0) {
    output$cellLineSelectorInfo <- renderText({
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
# Reactive cell line filtering
#----------------------------------------------------------------------------

cellLineSelectorModelList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$Sanger_model_ID
  
  #get cell lines with requested mutations
  cell_lines_mutation<-c()
  if(!is.null(input$cellLineSelectorGeneMutationSelect)){
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
  if(!is.null(input$cellLineSelectorGeneFusion3primeSelect)){
    
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
    if(!is.null(input$cellLineSelectorGeneFusion5primeSelect)){
      
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
  if(!is.null(input$cellLineSelectorGeneCNVSelect)){
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
  
  #get cell lines with requested gene dependency
  cell_lines_dependency<-c()
  if(!is.null(input$cellLineSelectorGeneDependencySelect)){
    
    gene_ids <- features %>%
      dplyr::filter(symbol %in% local(input$cellLineSelectorGeneDependencySelect)) %>%
      select(gene_id) %>%
      distinct() %>%
      collect() %>%
      .$gene_id
    
    contrasts_filtered <- contrasts %>%
      dplyr::filter(type == "dropout", dynamic_range <= -1.5, auc > 0.9) %>%
      distinct() %>%
      collect()
      
    contrast_ids <- con %>%
      tbl("gene_stats") %>%
      dplyr::filter(contrast_id %in% local(contrasts_filtered$contrast_id), gene_id %in% local(gene_ids), adjusted_effect <= local(input$cellLineSelectorGeneDependencySlider)) %>%
      dplyr::select(contrast_id, adjusted_effect) %>%
      dplyr::distinct() %>%
      collect() %>%
      group_by(contrast_id) %>%
      summarize(n_check=sum(adjusted_effect <= local(input$cellLineSelectorGeneDependencySlider)), n_genes=length(local(input$cellLineSelectorGeneDependencySelect))) %>%
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
      dplyr::select(Sanger_model_ID) %>%
      dplyr::distinct() %>%
      collect() %>%
      .$Sanger_model_ID
    
    model_ids <- model_ids[model_ids %in% cell_lines_dependency]
  }
  
  model_ids %>% unique
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
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    collect() %>%
    .$Sanger_model_ID
  
  con_cell_lines %>%
    tbl("cell_line_gene_mutations") %>%
    dplyr::filter(model_id %in% local(model_ids), gene_symbol %in% local(input$cellLineSelectorGeneMutationSelect)) %>%
    dplyr::select(protein_mutation) %>%
    distinct() %>%
    collect() %>%
    arrange(protein_mutation) %>%
    .$protein_mutation
    
})

cellLineSelectorGeneDependencyList <- reactive({
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  model_names <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(cell_line_name) %>%
    collect() %>%
    .$cell_line_name
  
  library_ids <- libraries %>%
    dplyr::filter(tissue_name %in% presel_tissue, cellline_name %in% local(model_names),  species == "human") %>%
    dplyr::select(library_id) %>%
    distinct %>%
    collect() %>%
    .$library_id
  
  gene_list_screens %>%
    dplyr::filter(library_id %in% local(library_ids)) %>%
    dplyr::select(symbol) %>%
    distinct() %>%
    collect() %>%
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

  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$Sanger_model_ID
  
  con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol_3prime) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol_3prime) %>%
    .$gene_symbol_3prime

})

cellLineSelectorGeneFusion5primeList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$Sanger_model_ID
  
  con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol_5prime) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol_5prime) %>%
    .$gene_symbol_5prime
    
})

cellLineSelectorGeneCNVList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }
  
  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    distinct() %>%
    collect() %>%
    .$Sanger_model_ID
  
  con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(gene_symbol) %>%
    distinct() %>%
    collect() %>%
    arrange(gene_symbol) %>%
    .$gene_symbol
    
})

cellLineSelectorGeneCNVCategoryList <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    presel_tissue <- tissue_list_cellLine
  }else{
    presel_tissue <- local(input$cellLineSelectorTissueSelect)
  }

  model_ids <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue) %>%
    dplyr::select(Sanger_model_ID) %>%
    collect() %>%
    .$Sanger_model_ID
  
  con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% local(model_ids)) %>%
    dplyr::select(cn_category) %>%
    distinct() %>%
    collect() %>%
    arrange(cn_category) %>%
    .$cn_category
})


# #----------------------------------------------------------------------------
# #  Observers
# #----------------------------------------------------------------------------
observeEvent(input$cellLineSelectorLoadButton, {
  output$cellLineSelectorDataTableMeta <- renderDataTable({
    dt <- cellLineSelectorDataTableMeta()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$cellLineSelectorDataTableScreens <- renderDataTable({
    dt <- cellLineSelectorDataTableScreens()
    if(!is.null(dt)){
      dt
    }
  })

  output$cellLineSelectorDataTableMutations <- renderDataTable({
    dt <- cellLineSelectorDataTableMutations()
    if(!is.null(dt)){
      dt
    }
  })

  output$cellLineSelectorDataTableFusions <- renderDataTable({
    dt <- cellLineSelectorDataTableFusions()
    if(!is.null(dt)){
      dt
    }
  })

  output$cellLineSelectorDataTableCNVs <- renderDataTable({
    dt <- cellLineSelectorDataTableCNVs()
    if(!is.null(dt)){
      dt
    }
  })
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
    enable("cellLineSelectorGeneFusion3primeSelect")
    enable("cellLineSelectorGeneFusion5primeSelect")
    enable("cellLineSelectorGeneCNVSelect")
    enable("cellLineSelectorGeneCNVCategorySelect")
  }else{
    disable("cellLineSelectorLoadButton")
    disable("cellLineSelectorGeneMutationSelect")
    disable("cellLineSelectorGeneDependencySelect")
    disable("cellLineSelectorGeneFusion3primeSelect")
    disable("cellLineSelectorGeneFusion5primeSelect")
    disable("cellLineSelectorGeneCNVSelect")
    disable("cellLineSelectorGeneCNVCategorySelect")
  }

  #update selectb
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationSelect', choices = gene_list_cellLine$gene_symbol, selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneDependencySelect', choices = cellLineSelectorGeneDependencyList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion3primeSelect', choices = cellLineSelectorGeneFusion3primeList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion5primeSelect', choices = cellLineSelectorGeneFusion5primeList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVSelect', choices = cellLineSelectorGeneCNVList(), selected=NULL,server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVCategorySelect', choices = cellLineSelectorGeneCNVCategoryList(), selected=NULL, server = TRUE)

  cellLineSelectorUpdateText()

}, ignoreNULL = FALSE)

observeEvent(input$cellLineSelectorCheckTissueAll, {
  if(isTRUE(input$cellLineSelectorCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'cellLineSelectorTissueSelect', choices = tissue_list_cellLine, selected=NULL, server = TRUE)
  }
  
  if((isTRUE(input$cellLineSelectorCheckTissueAll) | (!is.null(input$cellLineSelectorTissueSelect)))) {
    enable("cellLineSelectorGeneMutationSelect")
    enable("cellLineSelectorGeneDependencySelect")
    enable("cellLineSelectorGeneFusion3primeSelect")
    enable("cellLineSelectorGeneFusion5primeSelect")
    enable("cellLineSelectorGeneCNVSelect")
    enable("cellLineSelectorGeneCNVCategorySelect")
  }else{
    disable("cellLineSelectorGeneMutationSelect")
    disable("cellLineSelectorGeneDependencySelect")
    disable("cellLineSelectorGeneFusion3primeSelect")
    disable("cellLineSelectorGeneFusion5primeSelect")
    disable("cellLineSelectorGeneCNVSelect")
    disable("cellLineSelectorGeneCNVCategorySelect")
  }

  #update selectb
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationSelect', choices = gene_list_cellLine$gene_symbol, selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneDependencySelect', choices = cellLineSelectorGeneDependencyList(),  selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion3primeSelect', choices = cellLineSelectorGeneFusion3primeList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneFusion5primeSelect', choices = cellLineSelectorGeneFusion5primeList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVSelect', choices = cellLineSelectorGeneCNVList(), selected=NULL, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneCNVCategorySelect', choices = cellLineSelectorGeneCNVCategoryList(), selected=NULL, server = TRUE)

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
        mutations <- cellLineSelectorMutationsDataFrame()
        fusions <- cellLineSelectorFusionsDataFrame()
        cnvs <- cellLineSelectorCNVsDataFrame()

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
        shiny::incProgress(1/8)
        
        if("Screen results" %in% input$cellLineSelectorDownloadCheck & nrow(meta)>0){
          #write each sheet to a csv file, save the name
          table <- screen
          fileName <- "screen_data.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(2/8)

        if("Gene mutations" %in% input$cellLineSelectorDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          table <- mutations
          fileName <- "gene_mutations.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(3/8)

        if("Gene fusions" %in% input$cellLineSelectorDownloadCheck & nrow(fusions)>0){
          #write each sheet to a csv file, save the name
          table <- fusions
          fileName <- "gene_fusions.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(4/8)

        if("Gene CNVs" %in% input$cellLineSelectorDownloadCheck & nrow(cnvs)>0){
          #write each sheet to a csv file, save the name
          table <- cnvs
          fileName <- "gene_cnvs.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(5/8)

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