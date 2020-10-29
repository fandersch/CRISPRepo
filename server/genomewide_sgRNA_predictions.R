# ----------------------------------------------------------------------------
# genome-wide sgRNA predictions
# ----------------------------------------------------------------------------
sgRNAsUpdateText <- function(){
  output$sgRNAsInfo <- renderText({
    if(is.null(input$sgRNAsGeneSelect)){
      "INFO: Please select the gene(s) you want to browse in the right panel!"
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
  presel_genes <- input$sgRNAsGeneSelect %>% strsplit(split="\\(|\\)")
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
      filter(EntrezID %in% presel_gene_entrez) %>%
      select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      rbind(
        con_sgRNAs %>%
          tbl("sgRNAs_mouse") %>%
          filter(EntrezID %in% presel_gene_entrez) %>%
          select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
                 `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
                 Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
                 `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
                 `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
          collect()
        
      )
  }else{
    sgRNAs <- con_sgRNAs %>%
      tbl(tableSgRNAs) %>%
      filter(EntrezID %in% presel_gene_entrez) %>% 
      select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect()
  }
  if(sgRNAs$`Maps to genome` %>% as.character %>% unique %>% length == "1"){
    sgRNAs <- sgRNAs %>% select(-`Maps to genome`)
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
  #update gene selectbox
  #update other species selects
  updateSelectizeInput(session, 'gwsGeneSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse"), selected = input$sgRNAsSpeciesSelect, server = TRUE)
  updateSelectizeInput(session, 'gwsBrowseScreenSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse"), selected = input$sgRNAsSpeciesSelect, server = TRUE)
  updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), server = TRUE)
})

observeEvent(input$sgRNAsGeneSelect, {
  if((!is.null(input$sgRNAsGeneSelect))){
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
    paste0(paste(input$sgRNAsGeneSelect,collapse="_"), ".txt")
  },
  content = function(file) {
    sgRNAsTable() %>% write_tsv(file)
  }
)