# ----------------------------------------------------------------------------
# genome-wide sgRNA predictions
# ----------------------------------------------------------------------------
dualSgRNAsTopCombinationsGeneInputFile <- reactiveValues(data = NULL)


dualSgRNAsTopCombinationsUpdateText <- function(){
  output$dualSgRNAsTopCombinationsInfo <- renderText({
    if(is.null(input$dualSgRNAsTopCombinationsGeneSelect) & is.null(dualSgRNAsTopCombinationsGeneInputFile$data)){
      invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                       " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                       " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
    }else{
      "INFO: Click Load data!"
    }
  })
}
output$dualSgRNAsTopCombinationsInfo <- renderText({
  invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
})

#upon load display nothing
output$dualSgRNAsTopCombinationsTableOutput <- renderDataTable({
})

dualSgRNAsTopCombinationsTable <- reactive({
  
  if(!is.null(dualSgRNAsTopCombinationsGeneInputFile$data)){
    presel_genes_buff <- dualSgRNAsTopCombinationsGeneList()
    genes_fileUpload <- c(paste0("\\(", (dualSgRNAsTopCombinationsGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (dualSgRNAsTopCombinationsGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes <- input$dualSgRNAsTopCombinationsGeneSelect
  }
  
  presel_genes <- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")
  
  #specify which sgRNAs table to join
  if(input$dualSgRNAsTopCombinationsSpeciesSelect == "human"){
    tabledualSgRNAsTopCombinations <- "top20_dual_guides_human"
  }else{
    tabledualSgRNAsTopCombinations <- "top20_dual_guides_mouse"
  }
  
  if(input$dualSgRNAsTopCombinationsSpeciesSelect == "all"){
    sgRNAs <- con_sgRNAs %>%
      tbl("top20_dual_guides_human") %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
      collect() %>%
      rbind(
        con_sgRNAs %>%
          tbl("top20_dual_guides_mouse") %>%
          dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
          collect()
      )
  }else{
    sgRNAs <- con_sgRNAs %>%
      tbl(tabledualSgRNAsTopCombinations) %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>% 
      collect()
  }
  
  DBI::dbDisconnect(con_sgRNAs)
  
  sgRNAs %>%
    dplyr::rename(first_sgRNA_23mer = first_sgRNA, matching_sgRNA_23mer=matching_sgRNA) %>%
    mutate_at(c("first_sgRNA_MHstrength", "matching_sgRNA_MHStrength", "score_penalty", "first_guide_rank_score", "second_guide_rank_score", "score_combined"),  round, 2)
})

dualSgRNAsTopCombinationsDataTableOutput <- eventReactive(input$dualSgRNAsTopCombinationsLoadButton,{
  sgRNAs <- dualSgRNAsTopCombinationsTable()
  
  if (nrow(sgRNAs) > 0) {
    if(!is.null(input$dualSgRNAsTopCombinationsGeneSelect)){
      output$dualSgRNAsTopCombinationsInfo <- renderText({
        "INFO: Loading completed!"
      })
    }
    
    sgRNAs %>% 
      datatable(extensions = c('FixedColumns','FixedHeader'),
                options = list(
                  autoWidth = FALSE,
                  headerCallback = JS(headerCallback),
                  scrollX=TRUE,
                  fixedColumns = list(leftColumns = 3),
                  columnDefs = list(list(className = 'dt-center', targets = "_all")),
                  pageLength = 25,
                  lengthMenu = c(25, 50, 100, 200),
                  searchHighlight = TRUE
                ),
                filter = list(position = 'top', clear = FALSE),
                rownames= FALSE)
  }else{
    if(!is.null(input$dualSgRNAsTopCombinationsGeneSelect)){
      output$dualSgRNAsTopCombinationsInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      dualSgRNAsTopCombinationsUpdateText()
    }
  }
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

dualSgRNAsTopCombinationsGeneList <- reactive({
  if(!is.null(gene_list_mouse) & !is.null(gene_list_mouse)){
    if(input$dualSgRNAsTopCombinationsSpeciesSelect == "all"){
      gene_list <- gene_list_mouse %>%
        rbind(gene_list_human)
    }else{
      if(input$dualSgRNAsTopCombinationsSpeciesSelect == "human"){
        gene_list <- gene_list_human 
      }else{
        gene_list <- gene_list_mouse 
      }
    }
    
    gene_list %>%
      dplyr::mutate(gene = ifelse(is.na(Symbol), paste0("No symbol found (", EntrezID, ")"), paste0(Symbol , " (", EntrezID, ")"))) %>%
      arrange(gene) %>%
      .$gene
  }
  
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------
observeEvent(input$dualSgRNAsTopCombinationsLoadButton, {
  output$dualSgRNAsTopCombinationsTableOutput <- renderDataTable({
    dualSgRNAsTopCombinationsDatatableOutput <- dualSgRNAsTopCombinationsDataTableOutput()
  })
})

observeEvent(input$dualSgRNAsTopCombinationsSpeciesSelect, {
  #update other species selects
  updateSpecies(input$dualSgRNAsTopCombinationsSpeciesSelect)
  #update gene selectbox
  updateSelectizeInput(session, 'dualSgRNAsTopCombinationsGeneSelect', choices = dualSgRNAsTopCombinationsGeneList(), server = TRUE)
})

observeEvent(input$dualSgRNAsTopCombinationsGeneSelect, {
  if((!is.null(input$dualSgRNAsTopCombinationsGeneSelect)) | !is.null(dualSgRNAsTopCombinationsGeneInputFile$data)){
    enable("dualSgRNAsTopCombinationsLoadButton")
    if(!is.null(input$dualSgRNAsTopCombinationsGeneSelect)){
      dualSgRNAsTopCombinationsGeneInputFile$data <- NULL
      reset('dualSgRNAsTopCombinationsGeneInputFile')
    }
  }else{
    disable("dualSgRNAsTopCombinationsLoadButton")
  }
  dualSgRNAsTopCombinationsUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$dualSgRNAsTopCombinationsGeneInputFile, {
  if(!is.null(input$dualSgRNAsTopCombinationsGeneInputFile)){
    updateSelectizeInput(session, 'dualSgRNAsTopCombinationsGeneSelect', choices = dualSgRNAsTopCombinationsGeneList(), server = TRUE)
    req(input$dualSgRNAsTopCombinationsGeneInputFile)
    dualSgRNAsTopCombinationsGeneInputFile$data <- read_tsv(input$dualSgRNAsTopCombinationsGeneInputFile$datapath, col_names = F)
  }else{
    dualSgRNAsTopCombinationsGeneInputFile$data <- NULL
  }
  if((!is.null(input$dualSgRNAsTopCombinationsGeneSelect)) | !is.null(dualSgRNAsTopCombinationsGeneInputFile$data)){
    enable("dualSgRNAsTopCombinationsLoadButton")
  }else{
    disable("dualSgRNAsTopCombinationsLoadButton")
  }
  dualSgRNAsTopCombinationsUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$dualSgRNAsTopCombinationsButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_dualSgRNAsTopCombinations.txt"
  },
  content = function(file) {
    dualSgRNAsTopCombinationsTable() %>% write_tsv(file)
  }
)