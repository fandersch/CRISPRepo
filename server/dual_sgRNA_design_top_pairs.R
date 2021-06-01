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
  
  #specify which sgRNAs table to join
  if(input$dualSgRNAsTopCombinationsSpeciesSelect == "human"){
    tabledualSgRNAsTopCombinations <- "top10_dual_guides_human"
  }else{
    tabledualSgRNAsTopCombinations <- "top10_dual_guides_mouse"
  }
  
  
  if(input$dualSgRNAsTopCombinationsSpeciesSelect == "all"){
    sgRNAs <- con_sgRNAs %>%
      tbl("top10_dual_guides_human") %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
      dplyr::select(`Entrez ID` = EntrezID, 
                    Symbol,
                    first_sgRNA_sequence_23mer=original_sgRNA,
                    first_sgRNA_position_30mer=original_sgRNA_position,
                    first_sgRNA_exon_number=original_sgRNA_exon_number,
                    first_sgRNA_final_rank=original_final_rank,
                    first_sgRNA_validated_rank=original_final_validated_rank,
                    second_sgRNA_sequence_23mer=matching_sgRNA,
                    second_sgRNA_position_30mer=Position,
                    second_sgRNA_exon_number=exon,
                    second_sgRNA_final_rank=final_rank,
                    second_sgRNA_validated_rank=final_validated_rank,
                    second_sgRNA_VBCscore=`VBC.score`,
                    second_sgRNA_cleavage_activity=cleavage_activity,
                    second_sgRNA_Offtarget=Off_target,
                    targets_same_exon,
                    produces_frameshift,
                    proximity_1kb,
                    VBCscore_quartile=`VBC.score_top_group`,
                    cleavage_activity_quartile=`cleavage_act_top_group`,
                    produces_frameshift_in_most_transcripts,
                    produces_frameshift_in_all_transcripts,
                    most_common_dist,
                    transcripts_in_frame=TranscriptsInFrame,
                    n_transcripts_in_frame=nTranscriptsInFrame,
                    transcripts_out_of_frame=TranscriptsOutOfFrame,
                    n_transcripts_out_of_frame=nTranscriptsOutOfFrame,
                    combination_info=case,
                    penalty_score=score_penalty2,
                    combination_rank=score_rank) %>%
      collect() %>%
      rbind(
        con_sgRNAs %>%
          tbl("top10_dual_guides_mouse") %>%
          dplyr::filter(EntrezID %in% presel_gene_entrez) %>%
          dplyr::select(`Entrez ID` = EntrezID, 
                        Symbol,
                        first_sgRNA_sequence_23mer=original_sgRNA,
                        first_sgRNA_position_30mer=original_sgRNA_position,
                        first_sgRNA_exon_number=original_sgRNA_exon_number,
                        first_sgRNA_final_rank=original_final_rank,
                        first_sgRNA_validated_rank=original_final_validated_rank,
                        second_sgRNA_sequence_23mer=matching_sgRNA,
                        second_sgRNA_position_30mer=Position,
                        second_sgRNA_exon_number=exon,
                        second_sgRNA_final_rank=final_rank,
                        second_sgRNA_validated_rank=final_validated_rank,
                        second_sgRNA_VBCscore=`VBC.score`,
                        second_sgRNA_cleavage_activity=cleavage_activity,
                        second_sgRNA_Offtarget=Off_target,
                        targets_same_exon,
                        produces_frameshift,
                        proximity_1kb,
                        VBCscore_quartile=`VBC.score_top_group`,
                        cleavage_activity_quartile=`cleavage_act_top_group`,
                        produces_frameshift_in_most_transcripts,
                        produces_frameshift_in_all_transcripts,
                        most_common_dist,
                        transcripts_in_frame=TranscriptsInFrame,
                        n_transcripts_in_frame=nTranscriptsInFrame,
                        transcripts_out_of_frame=TranscriptsOutOfFrame,
                        n_transcripts_out_of_frame=nTranscriptsOutOfFrame,
                        combination_info=case,
                        penalty_score=score_penalty2,
                        combination_rank=score_rank) %>%
          collect()
      )
  }else{
    sgRNAs <- con_sgRNAs %>%
      tbl(tabledualSgRNAsTopCombinations) %>%
      dplyr::filter(EntrezID %in% presel_gene_entrez) %>% 
      dplyr::select(`Entrez ID` = EntrezID, 
                    Symbol,
                    first_sgRNA_sequence_23mer=original_sgRNA,
                    first_sgRNA_position_30mer=original_sgRNA_position,
                    first_sgRNA_exon_number=original_sgRNA_exon_number,
                    first_sgRNA_final_rank=original_final_rank,
                    first_sgRNA_validated_rank=original_final_validated_rank,
                    second_sgRNA_sequence_23mer=matching_sgRNA,
                    second_sgRNA_position_30mer=Position,
                    second_sgRNA_exon_number=exon,
                    second_sgRNA_final_rank=final_rank,
                    second_sgRNA_validated_rank=final_validated_rank,
                    second_sgRNA_VBCscore=`VBC.score`,
                    second_sgRNA_cleavage_activity=cleavage_activity,
                    second_sgRNA_Offtarget=Off_target,
                    targets_same_exon,
                    produces_frameshift,
                    proximity_1kb,
                    VBCscore_quartile=`VBC.score_top_group`,
                    cleavage_activity_quartile=`cleavage_act_top_group`,
                    produces_frameshift_in_most_transcripts,
                    produces_frameshift_in_all_transcripts,
                    most_common_dist,
                    transcripts_in_frame=TranscriptsInFrame,
                    n_transcripts_in_frame=nTranscriptsInFrame,
                    transcripts_out_of_frame=TranscriptsOutOfFrame,
                    n_transcripts_out_of_frame=nTranscriptsOutOfFrame,
                    combination_info=case,
                    penalty_score=score_penalty2,
                    combination_rank=score_rank) %>%
      collect()
  }
  
  sgRNAs
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