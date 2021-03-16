# ----------------------------------------------------------------------------
# genome-wide sgRNA predictions
# ----------------------------------------------------------------------------
correlationsGeneInputFile <- reactiveValues(data = NULL)

correlationsUpdateText <- function(){
  output$correlationsInfo <- renderText({
    if(is.null(input$correlationsGeneSelect) & is.null(correlationsGeneInputFile$data)){
      invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                       " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                       " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
    }else{
      "INFO: Click Load data!"
    }
  })
}
output$correlationsInfo <- renderText({
  invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
})

#upon load display nothing
output$correlationsDependencyExpressionTableOutput <- renderDataTable({
})

#upon load display nothing
output$correlationsCoEssentialityTableOutput <- renderDataTable({
})

#upon load display nothing
output$correlationsCoExpressionTableOutput <- renderDataTable({
})

#Dependency <> Expression
correlationsDependencyExpressionTable <- reactive({
  
  if(!is.null(correlationsGeneInputFile$data)){
    presel_genes_buff <- gene_list_correlations
    genes_fileUpload <- c(paste0("\\(", (correlationsGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (correlationsGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes <- input$correlationsGeneSelect
  }
  
  presel_genes <- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  correlationsDependencyExpression <- con_correlations %>%
    tbl("dependency_to_expression") %>%
    dplyr::filter(entrez_id_x %in% presel_gene_entrez) %>%
    collect() %>%
    group_by(entrez_id_x, dataset) %>%
    arrange(desc(abs(cor_coeff))) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 20 | abs(cor_coeff) >= input$correlationsSliderCoeff) %>%
    select(-rank) %>%
    ungroup %>%
    mutate(hit = ifelse(hit==1, "TRUE", "FALSE")) %>%
    dplyr::rename(number_Hits=nHits) %>%
    arrange(desc(number_Hits), symbol_y)


  correlationsDependencyExpression
})

#Co-essentiality
correlationsCoEssentialityTable <- reactive({
  
  if(!is.null(correlationsGeneInputFile$data)){
    presel_genes_buff <- gene_list_correlations
    genes_fileUpload <- c(paste0("\\(", (correlationsGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (correlationsGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes <- input$correlationsGeneSelect
  }
  
  presel_genes <- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  correlationsCoEssentiality <- con_correlations %>%
    tbl("co_essentiality") %>%
    dplyr::filter(entrez_id_x %in% presel_gene_entrez) %>%
    collect() %>%
    group_by(entrez_id_x) %>%
    arrange(desc(abs(cor_coeff))) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 20 | abs(cor_coeff) >= input$correlationsSliderCoeff) %>%
    select(-rank) %>%
    ungroup
  
  
  correlationsCoEssentiality
})

#Co-expression
correlationsCoExpressionTable <- reactive({
  
  if(!is.null(correlationsGeneInputFile$data)){
    presel_genes_buff <- gene_list_correlations
    genes_fileUpload <- c(paste0("\\(", (correlationsGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (correlationsGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes <- input$correlationsGeneSelect
  }
  
  presel_genes <- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  correlationsCoExpression <- con_correlations %>%
    tbl("co_expression") %>%
    dplyr::filter(entrez_id_x %in% presel_gene_entrez) %>%
    collect() %>%
    group_by(entrez_id_x) %>%
    arrange(desc(abs(cor_coeff))) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= 20 | abs(cor_coeff) >= input$correlationsSliderCoeff) %>%
    select(-rank) %>%
    ungroup
  
  
  correlationsCoExpression
})

correlationsDependencyExpressionTableOutput <- eventReactive(input$correlationsLoadButton,{
  correlationsDependencyExpression <- correlationsDependencyExpressionTable()
  
    if (nrow(correlationsDependencyExpression) > 0) {
      output$correlationsInfo <- renderText({
        "INFO: Loading completed!"
      })
      
      correlationsDependencyExpression %>% 
        datatable(extensions = c('FixedColumns','FixedHeader'),
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

correlationsCoEssentialityTableOutput <- eventReactive(input$correlationsLoadButton,{
  correlationsCoEssentiality <- correlationsCoEssentialityTable()
  
  if (nrow(correlationsCoEssentiality) > 0) {
    output$correlationsInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    correlationsCoEssentiality %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
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

correlationsCoExpressionTableOutput <- eventReactive(input$correlationsLoadButton,{
  correlationsCoExpression <- correlationsCoExpressionTable()
  
  if (nrow(correlationsCoExpression) > 0) {
    output$correlationsInfo <- renderText({
      "INFO: Loading completed!"
    })
    
    correlationsCoExpression %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
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

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------
observeEvent(input$correlationsLoadButton, {
  output$correlationsDependencyExpressionTableOutput <- renderDataTable({
    correlationsDependencyExpressionTableOutput()
  })
  
  output$correlationsCoEssentialityTableOutput <- renderDataTable({
    correlationsCoEssentialityTableOutput()
  })
  
  output$correlationsCoExpressionTableOutput <- renderDataTable({
    correlationsCoExpressionTableOutput()
  })
})

observeEvent(input$correlationsGeneSelect, {
  if((!is.null(input$correlationsGeneSelect)) | !is.null(correlationsGeneInputFile$data)){
    enable("correlationsLoadButton")
    if(!is.null(input$correlationsGeneSelect)){
      reset('correlationsGeneInputFile')
      correlationsGeneInputFile$data <- NULL
    }
  }else{
    disable("correlationsLoadButton")
  }
  correlationsUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$correlationsGeneInputFile, {
  if(!is.null(input$correlationsGeneInputFile)){
    updateSelectizeInput(session, 'correlationsGeneSelect', choices = gene_list_correlations, server = TRUE)
    req(input$correlationsGeneInputFile)
    correlationsGeneInputFile$data <- read_tsv(input$correlationsGeneInputFile$datapath, col_names = F)
  }else{
    correlationsGeneInputFile$data <- NULL
  }
  if((!is.null(input$correlationsGeneSelect)) | !is.null(correlationsGeneInputFile$data)){
    enable("correlationsLoadButton")
  }else{
    disable("correlationsLoadButton")
  }
  correlationsUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$correlationsDownloadCheck, {
  if(!is.null(input$correlationsDownloadCheck)){
    enable("correlationsButtonDownload")
  }else{
    disable("correlationsButtonDownload")
  }
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$correlationsButtonDownload <- downloadHandler(
  filename = function() {
    paste0("correlations_", paste(input$correlationsGeneSelect,collapse="_"), ".zip")
  },
  content = function(file) {
    
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    files <- NULL;
    
    if("Dependency <> Expression" %in% input$correlationsDownloadCheck & nrow(correlationsDependencyExpressionTable())>0){
      #write each sheet to a csv file, save the name
      table <- correlationsDependencyExpressionTable()
      fileName <- str_replace(str_replace(paste0(paste(input$correlationsGeneSelect,collapse="_"), "_dependency_expression.txt"), " (", "_"), ")", "")
      write.table(table,fileName, row.names = F, col.names = T)
      files <- c(fileName,files)
    }
    
    if("Co-Essentiality" %in% input$correlationsDownloadCheck & nrow(correlationsCoEssentialityTable())>0){
      #write each sheet to a csv file, save the name
      table <- correlationsCoEssentialityTable()
      fileName <- str_replace(str_replace(paste0(paste(input$correlationsGeneSelect,collapse="_"), "_co_essentiality.txt"), " (", "_"), ")", "")
      write.table(table,fileName, row.names = F, col.names = T)
      files <- c(fileName,files)
    }
    
    if("Co-Expression" %in% input$correlationsDownloadCheck & nrow(correlationsCoExpressionTable())>0){
      #write each sheet to a csv file, save the name
      table <- correlationsCoExpressionTable()
      fileName <- str_replace(str_replace(paste0(paste(input$correlationsGeneSelect,collapse="_"), "_co_expression.txt"), " (", "_"), ")", "")
      write.table(table,fileName, row.names = F, col.names = T)
      files <- c(fileName,files)
    }
    
    if(!is.null(files)){
      #create the zip file
      zip(file,files)
    }else{
      NULL
    }
  }
)
  
output$correlationsButtonDownloadPrimaryTables <- downloadHandler(
  filename = function() {
    "correlations_primary_data.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading", input$dataset, " Data"),
      value = 0,
      {
        files <- NULL;
        if("Dependency <> Expression" %in% input$correlationsDownloadPrimaryTablesCheck){
          fileName <- "correlations_primary_data/correlations_primary_data_dependency_expression.txt"
          files <- c(fileName,files)
        }
        if("Co-Essentiality" %in% input$correlationsDownloadPrimaryTablesCheck){
          fileName <- "correlations_primary_data/correlations_primary_data_co_essentiality.txt"
          files <- c(fileName,files)
        }
        if("Co-Expression" %in% input$correlationsDownloadPrimaryTablesCheck){
          fileName <- "correlations_primary_data/correlations_primary_data_co_expression.txt"
          files <- c(fileName,files)
        }
        shiny::incProgress(1/2)
        if(!is.null(files)){
          #create the zip file
          zip(file,files)
        }else{
          NULL
        }
      }
    )
  }
)