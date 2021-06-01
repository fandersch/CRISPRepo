# ----------------------------------------------------------------------------
# genome-wide sgRNA predictions
# ----------------------------------------------------------------------------
sgRNAsGeneInputFile <- reactiveValues(data = NULL)


sgRNAsUpdateText <- function(){
  output$sgRNAsInfo <- renderText({
    if(is.null(input$sgRNAsGeneSelect) & is.null(sgRNAsGeneInputFile$data)){
      invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                       " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                       " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
    }else{
      "INFO: Click Load data!"
    }
  })
}
output$sgRNAsInfo <- renderText({
  invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
})

#upon load display nothing
output$sgRNAsTableOutput <- renderDataTable({
})

sgRNAsTable <- reactive({
  
  if(!is.null(sgRNAsGeneInputFile$data)){
    presel_genes_buff <- sgRNAsGeneList()
    genes_fileUpload <- c(paste0("\\(", (sgRNAsGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (sgRNAsGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes <- input$sgRNAsGeneSelect
  }
  
  presel_genes <- presel_genes %>% strsplit(split="\\(|\\)")
  presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
  presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
  
  #specify which sgRNAs table to join
  if(input$sgRNAsSpeciesSelect == "human"){
    tableSgRNAs <- "sgRNAs_human"
  }else{
    tableSgRNAs <- "sgRNAs_mouse"
  }
  
  
  if(input$sgRNAsSpeciesSelect == "all"){
    sgRNAs <- con_sgRNAs %>%
      tbl("sgRNAs_human") %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
      # dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, `30-mer Position` = Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
      #        `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
      #        Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
      #        `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
      #        `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      rbind(
        con_sgRNAs %>%
          tbl("sgRNAs_mouse") %>%
          dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
          # dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, `30-mer Position` = Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
          #        `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
          #        Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
          #        `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
          #        `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
          collect()
        
      )
  }else{
    sgRNAs <- con_sgRNAs %>%
      tbl(tableSgRNAs) %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>% 
      # dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, `30-mer Position` = Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
      #        `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
      #        Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
      #        `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
      #        `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect()
  }
  
  sgRNAs <- sgRNAs %>% 
    rowwise() %>%
    dplyr::mutate(orientation=stringr::str_split(Position, pattern = "[()]")[[1]][2]) %>%
    dplyr::mutate(end=ifelse(orientation=="+", stringr::str_split(Position, pattern = "[-:(]")[[1]][3], stringr::str_split(Position, pattern = "[-:(]")[[1]][2])) %>%
    dplyr::mutate(cutting_position = ifelse(orientation=="+", as.numeric(end) - 3 - 3 - 3, as.numeric(end) + 3 + 3 + 3)) %>%
    dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, `30-mer Position` = Position, `Genomi cutting position`= cutting_position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
           `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target,
           Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
           `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation,
           `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check)
  
  
  if(sgRNAs$`Maps to genome` %>% as.character %>% unique %>% length == "1"){
    sgRNAs <- sgRNAs %>% dplyr::select(-`Maps to genome`)
  }
  
  sgRNAs
})

sgRNAsDataTableOutput <- eventReactive(input$sgRNAsLoadButton,{
  sgRNAs <- sgRNAsTable()
  
    if (nrow(sgRNAs) > 0) {
      if(!is.null(input$sgRNAsGeneSelect)){
        output$sgRNAsInfo <- renderText({
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
    if(!is.null(input$sgRNAsGeneSelect)){
      output$sgRNAsInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      sgRNAsUpdateText()
    }
  }
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

sgRNAsGeneList <- reactive({
  if(!is.null(gene_list_mouse) & !is.null(gene_list_mouse)){
    if(input$sgRNAsSpeciesSelect == "all"){
      gene_list <- gene_list_mouse %>%
        rbind(gene_list_human)
    }else{
      if(input$sgRNAsSpeciesSelect == "human"){
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
observeEvent(input$sgRNAsLoadButton, {
  output$sgRNAsTableOutput <- renderDataTable({
    sgRNAsDatatableOutput <- sgRNAsDataTableOutput()
  })
})

observeEvent(input$sgRNAsSpeciesSelect, {
  #update other species selects
  updateSpecies(input$sgRNAsSpeciesSelect)
  #update gene selectbox
  updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), server = TRUE)
  
})

observeEvent(input$sgRNAsGeneSelect, {
  if((!is.null(input$sgRNAsGeneSelect)) | !is.null(sgRNAsGeneInputFile$data)){
    enable("sgRNAsLoadButton")
    if(!is.null(input$sgRNAsGeneSelect)){
      sgRNAsGeneInputFile$data <- NULL
      reset('sgRNAsGeneInputFile')
    }
  }else{
    disable("sgRNAsLoadButton")
  }
  sgRNAsUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAsGeneInputFile, {
  if(!is.null(input$sgRNAsGeneInputFile)){
    updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), server = TRUE)
    req(input$sgRNAsGeneInputFile)
    sgRNAsGeneInputFile$data <- read_tsv(input$sgRNAsGeneInputFile$datapath, col_names = F)
  }else{
    sgRNAsGeneInputFile$data <- NULL
  }
  if((!is.null(input$sgRNAsGeneSelect)) | !is.null(sgRNAsGeneInputFile$data)){
    enable("sgRNAsLoadButton")
  }else{
    disable("sgRNAsLoadButton")
  }
  sgRNAsUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$sgRNAsButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_sgRNAs.txt"
  },
  content = function(file) {
    sgRNAsTable() %>% write_tsv(file)
  }
)