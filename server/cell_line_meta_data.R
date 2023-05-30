# ----------------------------------------------------------------------------
# cellLine
# ----------------------------------------------------------------------------
cellLine_geneInputFile <- reactiveValues(data = NULL)


cellLineUpdateText <- function(){
  output$cellLineInfo <- renderText({
    if(is.null(input$cellLineTissueSelect) & !isTRUE(input$cellLineCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) in the right panel!")
    }else{
      if(is.null(input$cellLineCellLineSelect) & !isTRUE(input$cellLineCheckCellLineAll)){
        invisible("INFO: Please select the cell line(s) in the right panel!")
      }else{
        if(is.null(input$cellLineGeneSelect) & !isTRUE(input$cellLineCheckGeneAll) & is.null(cellLine_geneInputFile$data)){
          invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                    " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                    " The file must have one gene per row (one single column): gene-symbol (no header)!"))
        }else{
          invisible("INFO: Click Load data!")
        }
      }
    }
  })
}

#upon load display nothing
output$cellLineDataTableMeta <- renderDataTable({
})
output$cellLineDataTableMutations <- renderDataTable({
})
output$cellLineDataTableFusions <- renderDataTable({
})
output$cellLineDataTableCNVs <- renderDataTable({
})

#query database and create dataframe
cellLineMetaDataFrame <- reactive({

  #get selected tissue
  if(isTRUE(input$cellLineCheckTissueAll)){
    presel_tissue <- cellLineTissueList()
  }else{
    presel_tissue <- local(input$cellLineTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$cellLineCheckCellLineAll)){
    presel_cell_line <- cellLineCellLineList()
  }else{
    presel_cell_line <- local(input$cellLineCellLineSelect)
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLine_meta_data <- con_cell_lines %>%
    tbl("cell_line_meta") %>%
    dplyr::filter(tissue_name %in% presel_tissue, cell_line_name %in% presel_cell_line) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLine_meta_data
})

cellLineMutationsDataFrame <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineCheckTissueAll)){
    presel_tissue <- cellLineTissueList()
  }else{
    presel_tissue <- local(input$cellLineTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$cellLineCheckCellLineAll)){
    presel_cell_line <- cellLineCellLineList()
  }else{
    presel_cell_line <- local(input$cellLineCellLineSelect)
  }
  
  #get model_id from selected cell lines
  presel_model_ids <- cellline_list_cellLine %>%
    filter(cell_line_name %in% presel_cell_line) %>%
    .$model_id
  
  if(!is.null(cellLine_geneInputFile$data)){
    presel_genes_buff <- cellLineGeneList()
    genes_fileUpload <- c(cellLine_geneInputFile$data$X1 %>% as.character)
    presel_genes <- grep(paste(paste0("^", genes_fileUpload, "$"),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$cellLineCheckGeneAll)){
      presel_genes <- cellLineGeneList()
    }else{
      presel_genes <- local(input$cellLineGeneSelect)
    }
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLine_meta_data <- con_cell_lines %>%
    tbl("cell_line_gene_mutations") %>%
    dplyr::filter(model_id %in% presel_model_ids, gene_symbol %in% presel_genes) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLine_meta_data
})

cellLineFusionsDataFrame <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineCheckTissueAll)){
    presel_tissue <- cellLineTissueList()
  }else{
    presel_tissue <- local(input$cellLineTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$cellLineCheckCellLineAll)){
    presel_cell_line <- cellLineCellLineList()
  }else{
    presel_cell_line <- local(input$cellLineCellLineSelect)
  }
  
  #get model_id from selected cell lines
  presel_model_ids <- cellline_list_cellLine %>%
    filter(cell_line_name %in% presel_cell_line) %>%
    .$model_id
  
  if(!is.null(cellLine_geneInputFile$data)){
    presel_genes_buff <- cellLineGeneList()
    genes_fileUpload <- c(cellLine_geneInputFile$data$X1 %>% as.character)
    presel_genes <- grep(paste(paste0("^", genes_fileUpload, "$"),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$cellLineCheckGeneAll)){
      presel_genes <- cellLineGeneList()
    }else{
      presel_genes <- local(input$cellLineGeneSelect)
    }
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLine_meta_data <- con_cell_lines %>%
    tbl("cell_line_gene_fusions") %>%
    dplyr::filter(model_id %in% presel_model_ids, gene_symbol_3prime %in% presel_genes | gene_symbol_5prime %in% presel_genes) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLine_meta_data
})

cellLineCNVsDataFrame <- reactive({
  
  #get selected tissue
  if(isTRUE(input$cellLineCheckTissueAll)){
    presel_tissue <- cellLineTissueList()
  }else{
    presel_tissue <- local(input$cellLineTissueSelect)
  }
  
  #get selected cell line
  if(isTRUE(input$cellLineCheckCellLineAll)){
    presel_cell_line <- cellLineCellLineList()
  }else{
    presel_cell_line <- local(input$cellLineCellLineSelect)
  }
  
  #get model_id from selected cell lines
  presel_model_ids <- cellline_list_cellLine %>%
    filter(cell_line_name %in% presel_cell_line) %>%
    .$model_id
  
  if(!is.null(cellLine_geneInputFile$data)){
    presel_genes_buff <- cellLineGeneList()
    genes_fileUpload <- c(cellLine_geneInputFile$data$X1 %>% as.character)
    presel_genes <- grep(paste(paste0("^", genes_fileUpload, "$"),collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    if(isTRUE(input$cellLineCheckGeneAll)){
      presel_genes <- cellLineGeneList()
    }else{
      presel_genes <- local(input$cellLineGeneSelect)
    }
  }
  
  con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")
  
  cellLine_meta_data <- con_cell_lines %>%
    tbl("cell_line_gene_cnv") %>%
    dplyr::filter(model_id %in% presel_model_ids, gene_symbol %in% presel_genes) %>%
    collect()
  
  DBI::dbDisconnect(con_cell_lines)
  
  cellLine_meta_data
})

#create datatable out of dataframe
cellLineDataTableMeta <- eventReactive(input$cellLineLoadButton,{
  
  df <- cellLineMetaDataFrame()
  
  if (nrow(df) > 0) {
    output$cellLineInfo <- renderText({
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

cellLineDataTableMutations <- eventReactive(input$cellLineLoadButton,{
  
  df <- cellLineMutationsDataFrame()
  
  if (nrow(df) > 0) {
    output$cellLineInfo <- renderText({
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

cellLineDataTableFusions <- eventReactive(input$cellLineLoadButton,{
  
  df <- cellLineFusionsDataFrame()
  
  if (nrow(df) > 0) {
    output$cellLineInfo <- renderText({
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

cellLineDataTableCNVs <- eventReactive(input$cellLineLoadButton,{
  
  df <- cellLineCNVsDataFrame()
  
  if (nrow(df) > 0) {
    output$cellLineInfo <- renderText({
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

cellLineTissueList <- reactive({
  
  cellline_list_cellLine %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

cellLineCellLineList <- reactive({
  
  preselTissue = cellline_list_cellLine %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    .$tissue_name
  
  if(!isTRUE(input$cellLineCheckTissueAll) & !is.null(input$cellLineTissueSelect)){
    preselTissue = cellline_list_cellLine %>%
      dplyr::filter(tissue_name %in% input$cellLineTissueSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
  }
  
  cellline_list_cellLine %>%
    dplyr::filter(tissue_name %in% preselTissue) %>%
    dplyr::select(cell_line_name) %>%
    arrange(cell_line_name) %>%
    .$cell_line_name
})

cellLineGeneList <- reactive({
  if(!isTRUE(input$cellLineCheckCellLineAll) & (is.null(input$cellLineCellLineSelect))){
    NULL
  }else{
    
  gene_list_cellLine %>%
      dplyr::select(gene_symbol) %>%
      arrange(gene_symbol) %>%
      .$gene_symbol
  }
})


#----------------------------------------------------------------------------
#  Observers
#----------------------------------------------------------------------------
observeEvent(input$cellLineLoadButton, {
  output$cellLineDataTableMeta <- renderDataTable({
    dt <- cellLineDataTableMeta()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$cellLineDataTableMutations <- renderDataTable({
    dt <- cellLineDataTableMutations()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$cellLineDataTableFusions <- renderDataTable({
    dt <- cellLineDataTableFusions()
    if(!is.null(dt)){
      dt
    }
  })
  
  output$cellLineDataTableCNVs <- renderDataTable({
    dt <- cellLineDataTableCNVs()
    if(!is.null(dt)){
      dt
    }
  })
})

observeEvent(input$cellLineTissueSelect, {
  if(!is.null(input$cellLineTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'cellLineCheckTissueAll', value = FALSE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'cellLineCellLineSelect', choices = cellLineCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'cellLineCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$cellLineCheckTissueAll) | (!is.null(input$cellLineTissueSelect))){
    enable("cellLineCellLineSelect")
    enable("cellLineCheckCellLineAll")
  }else{
    disable("cellLineCellLineSelect")
    disable("cellLineCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'cellLineGeneSelect', choices = cellLineGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
  
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLineCheckTissueAll, {
  if(isTRUE(input$cellLineCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'cellLineTissueSelect', choices = cellLineTissueList(), server = TRUE)
  }
  
  #update cell line selectbox
  updateSelectizeInput(session, 'cellLineCellLineSelect', choices = cellLineCellLineList(), server = TRUE)
  #select cell line checkbox
  updateCheckboxInput(session, 'cellLineCheckCellLineAll', value = FALSE)
  #update cell line selectb
  if(isTRUE(input$cellLineCheckTissueAll) | (!is.null(input$cellLineTissueSelect))){
    enable("cellLineCellLineSelect")
    enable("cellLineCheckCellLineAll")
  }else{
    disable("cellLineCellLineSelect")
    disable("cellLineCheckCellLineAll")
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'cellLineGeneSelect', choices = cellLineGeneList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
  
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLineCellLineSelect, {
  if(!is.null(input$cellLineCellLineSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'cellLineCheckCellLineAll', value = FALSE)
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'cellLineGeneSelect', choices = cellLineGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$cellLineCheckCellLineAll) | (!is.null(input$cellLineCellLineSelect))) & (isTRUE(input$cellLineCheckTissueAll) | (!is.null(input$cellLineTissueSelect)))) {
    enable("cellLineGeneSelect")
    enable("cellLine_inputFile")
    enable("cellLineCheckGeneAll")
  }else{
    disable("cellLineGeneSelect")
    disable("cellLine_inputFile")
    disable("cellLineCheckGeneAll")
  }
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLineCheckCellLineAll, {
  if(isTRUE(input$cellLineCheckCellLineAll)){
    #update library selectbox
    updateSelectizeInput(session, 'cellLineCellLineSelect', choices = cellLineCellLineList(), server = TRUE)
    
    #get selected tissue
    if(isTRUE(input$cellLineCheckTissueAll)){
      presel_tissue <- cellLineTissueList()
    }else{
      presel_tissue <- local(input$cellLineTissueSelect)
    }
    
    #get selected cell line
    if(isTRUE(input$cellLineCheckCellLineAll)){
      presel_cell_line <- cellLineCellLineList()
    }else{
      presel_cell_line <- local(input$cellLineCellLineSelect)
    }
  }
  
  #update gene selectb
  updateSelectizeInput(session, 'cellLineGeneSelect', choices = cellLineGeneList(), server = TRUE)
  #unselect checkbox gene
  updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
  
  if((isTRUE(input$cellLineCheckCellLineAll) | (!is.null(input$cellLineCellLineSelect))) & (isTRUE(input$cellLineCheckTissueAll) | (!is.null(input$cellLineTissueSelect)))) {
    enable("cellLineGeneSelect")
    enable("cellLine_inputFile")
    enable("cellLineCheckGeneAll")
  }else{
    disable("cellLineGeneSelect")
    disable("cellLine_inputFile")
    disable("cellLineCheckGeneAll")
  }
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLineCancelModal, {
  updateCheckboxInput(session, 'cellLineCheckCellLineAll', value = FALSE)
  removeModal()
})

observeEvent(input$cellLineGeneSelect, {
  if(!is.null(input$cellLineGeneSelect)){
    updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
    reset('cellLine_inputFile')
    cellLine_geneInputFile$data <- NULL
  }
  if(isTRUE(input$cellLineCheckGeneAll) | (!is.null(input$cellLineGeneSelect)) | !is.null(cellLine_geneInputFile$data)){
    enable("cellLineLoadButton")
  }else{
    disable("cellLineLoadButton")
  }
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLine_inputFile, {
  if(!is.null(input$cellLine_inputFile)){
    updateSelectizeInput(session, 'cellLineGeneSelect', choices = cellLineGeneList(), server = TRUE)
    updateCheckboxInput(session, 'cellLineCheckGeneAll', value = FALSE)
    req(input$cellLine_inputFile)
    cellLine_geneInputFile$data <- read_tsv(input$cellLine_inputFile$datapath, col_names = F)
  }else{
    cellLine_geneInputFile$data <- NULL
  }
  if(isTRUE(input$cellLineCheckGeneAll) | (!is.null(input$cellLineGeneSelect)) | !is.null(cellLine_geneInputFile$data)){
    enable("cellLineLoadButton")
  }else{
    disable("cellLineLoadButton")
  }
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$cellLineCheckGeneAll, {
  if(isTRUE(input$cellLineCheckGeneAll)){
    geneList <- cellLineGeneList()
    updateSelectizeInput(session, 'cellLineGeneSelect', choices = geneList, server = TRUE)
    reset('cellLine_inputFile')
    cellLine_geneInputFile$data <- NULL
  }
  if(isTRUE(input$cellLineCheckGeneAll) | (!is.null(input$cellLineGeneSelect)) | !is.null(cellLine_geneInputFile$data)){
    enable("cellLineLoadButton")
  }else{
    disable("cellLineLoadButton")
  }
  cellLineUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$cellLineButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_cellline_meta_data.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading data"),
      value = 0,
      {
        meta <- cellLineMetaDataFrame()
        mutations <- cellLineMutationsDataFrame()
        fusions <- cellLineFusionsDataFrame()
        cnvs <- cellLineCNVsDataFrame()
        
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;
        
        if("Cellline meta data" %in% input$cellLineDownloadCheck & nrow(meta)>0){
          #write each sheet to a csv file, save the name
          table <- meta
          fileName <- "cellline_meta_data.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(1/6)
        
        if("Gene mutations" %in% input$cellLineDownloadCheck & nrow(mutations)>0){
          #write each sheet to a csv file, save the name
          table <- mutations
          fileName <- "gene_mutations.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(2/6)
        
        if("Gene fusions" %in% input$cellLineDownloadCheck & nrow(fusions)>0){
          #write each sheet to a csv file, save the name
          table <- fusions
          fileName <- "gene_fusions.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(3/6)
        
        if("Gene CNVs" %in% input$cellLineDownloadCheck & nrow(cnvs)>0){
          #write each sheet to a csv file, save the name
          table <- cnvs
          fileName <- "gene_cnvs.txt"
          write.table(table,fileName, row.names = F, col.names = T)
          files <- c(fileName,files)
        }
        shiny::incProgress(4/6)
        
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