# ----------------------------------------------------------------------------
# sgRNA Info
# ----------------------------------------------------------------------------

sgRNAInfoUpdateText <- function(){
  output$sgRNAInfoInfo <- renderText({
    if(is.null(input$sgRNAInfoSelectGene)){
      "INFO: Please select the gene(s) of considered sgRNAs in the right panel!"
    }else{
      if(is.null(input$sgRNAInfoSelectGuide) & !isTRUE(input$sgRNAInfoCheckGuideAll)){
        "INFO: Please select the considered sgRNAs in the right panel!"
      }else{
        "INFO: Click Load data!"
      }
    }
  })
}

output$sgRNAInfoTableOutputScreens <- renderDataTable({
})

output$sgRNAInfoTableOutputPredictions <- renderDataTable({
})

output$sgRNAInfoTableOutputValidations <- renderDataTable({
})

#data load
sgRNAInfoTableScreens <- reactive({
  
  #specify which sgRNAs table to join
  if(input$sgRNAInfoSpeciesSelect == "human"){
    tableSgRNAs <- "sgRNAs_human"
  }else{
    tableSgRNAs <- "sgRNAs_mouse"
  }
  
  #retrieve selected genes
  presel_genes_both<- input$sgRNAInfoSelectGene %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  if(isTRUE(input$sgRNAInfoCheckGuideAll)){
    presel_guides <- sgRNAInfoGuideList()
  }else{
    presel_guides <- input$sgRNAInfoSelectGuide
  }
  
  #get gene_stats/guide_stats
  df_screens <- con %>%
    tbl("guide_stats") %>%
    filter(guide_id %in% local(presel_guides)) %>%
    left_join(features %>% select(-library_id) %>% distinct) %>%
    collect() %>%
    mutate_at(c("lfc","effect"), funs(round(., 3)))
  
  if (nrow(df_screens) > 0) {
    presel_contrasts <- df_screens$contrast_id %>% unique
    
    df_screens <- df_screens %>%
      select(contrast_id, guide_id, entrez_id, symbol,local(input$sgRNAInfoIndexRadio)) %>%
      spread(contrast_id, input$sgRNAInfoIndexRadio) %>%
      distinct()
  }
  df_screens
  
})

sgRNAInfoTablePredictions <- reactive({
  
  #specify which sgRNAs table to join
  if(input$sgRNAInfoSpeciesSelect == "human"){
    tableSgRNAs <- "sgRNAs_human"
  }else{
    tableSgRNAs <- "sgRNAs_mouse"
  }
  
  #retrieve selected genes
  presel_genes_both<- input$sgRNAInfoSelectGene %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  if(isTRUE(input$sgRNAInfoCheckGuideAll)){
    presel_guides <- sgRNAInfoGuideList()
  }else{
    presel_guides <- input$sgRNAInfoSelectGuide
  }
  
  sgRNAs_23mer <- features %>%
    filter(guide_id %in% presel_guides) %>%
    dplyr::mutate(sgRNA_23mer = substr(context, 5, nchar(context)-3)) %>%
    .$sgRNA_23mer
  
  if(input$sgRNAInfoSpeciesSelect == "all"){
    sgRNAs <- con_sgRNAs %>%
      tbl("sgRNAs_human") %>%
      filter(EntrezID %in% c(presel_entrez)) %>%
      select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      filter(str_detect(`20-mer + NGG`, paste(sgRNAs_23mer, collapse = "|"))) %>%
      rbind(
        con_sgRNAs %>%
          tbl("sgRNAs_mouse") %>%
          filter(EntrezID %in% c(presel_entrez)) %>%
          select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
                 `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
                 Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
                 `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
                 `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
          collect() %>%
          filter(str_detect(`20-mer + NGG`, paste(sgRNAs_23mer, collapse = "|")))
      )
  }else{
    #get guide ranks
    sgRNAs <- con_sgRNAs %>%
      tbl(tableSgRNAs) %>%
      filter(EntrezID %in% c(presel_entrez)) %>%
      select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      filter(str_detect(`20-mer + NGG`, paste(sgRNAs_23mer, collapse = "|")))
  }
  
  if(sgRNAs$`Maps to genome` %>% as.character %>% unique %>% length == "1"){
    sgRNAs <- sgRNAs %>% select(-`Maps to genome`)
  }
  
  sgRNAs
  
})

sgRNAInfoTableValidations <- reactive({
  
  #specify which sgRNAs table to join
  if(input$sgRNAInfoSpeciesSelect == "all"){
    tableSgRNAs <- c("validated_guides_dropoutHQ_human", "validated_guides_dropoutBroadAvana_human", "validated_guides_drugModifier_human", "validated_guides_facsBased_human", "validated_guides_dropoutHQ_mouse", "validated_guides_drugModifier_mouse", "validated_guides_facsBased_mouse")
  }else{
    if(input$sgRNAInfoSpeciesSelect == "human"){
      tableSgRNAs <- c("validated_guides_dropoutHQ_human", "validated_guides_dropoutBroadAvana_human", "validated_guides_drugModifier_human", "validated_guides_facsBased_human")
    }else{
      tableSgRNAs <- c("validated_guides_dropoutHQ_mouse", "validated_guides_drugModifier_mouse", "validated_guides_facsBased_mouse")
    }
  }
  
  #retrieve selected genes
  presel_genes_both<- input$sgRNAInfoSelectGene %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  if(isTRUE(input$sgRNAInfoCheckGuideAll)){
    presel_guides <- sgRNAInfoGuideList()
  }else{
    presel_guides <- input$sgRNAInfoSelectGuide
  }
  
  sgRNAs_23mer <- features %>%
    filter(guide_id %in% presel_guides) %>%
    dplyr::mutate(sgRNA_23mer = substr(context, 5, nchar(context)-3)) %>%
    .$sgRNA_23mer
  
  for(i in 1:length(tableSgRNAs)){
    sgRNAs_buff <- con_sgRNAs %>%
      tbl(tableSgRNAs[i]) %>%
      filter(EntrezID %in% c(presel_entrez)) %>%
      collect() %>%
      mutate(dataset = unlist(strsplit(tableSgRNAs[i], split="_"))[3:4] %>% paste(collapse = "_")) %>%
      separate(sgRNA_id, into =c("dummy", "seq_buff"), sep = "_", remove = FALSE) %>%
      filter(str_detect(seq_buff, paste(sgRNAs_23mer, collapse = "|"))) %>%
      select(-seq_buff, -dummy) %>%
      select(dataset, everything())
    
    if(!"outlier_high_confident_all_total" %in% colnames(sgRNAs_buff)){
      sgRNAs_buff <- sgRNAs_buff %>%
        mutate(outlier_high_confident_all_total = NA)
    }
    if(!exists("sgRNAs")){
      if(sgRNAs_buff %>% nrow > 0){
        sgRNAs <- sgRNAs_buff
      }
    }else{
      if(sgRNAs_buff %>% nrow > 0){
        sgRNAs <- sgRNAs %>% 
          rbind(sgRNAs_buff)
      }
    }
  }
  if(exists("sgRNAs")){
    sgRNAs
  }else{
    sgRNAs_buff
  }
  
})

sgRNAInfoDataTableScreens <- eventReactive(input$sgRNAInfoLoadButton,{
  df <- sgRNAInfoTableScreens()
  
  if (nrow(df) > 0) {
    if(!is.null(input$sgRNAInfoSelectGene)){
      output$sgRNAInfoInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    
    effect<-FALSE
    if(input$sgRNAInfoIndexRadio == "effect"){
      effect<-TRUE
    }
    
    sgRNAInfoDatatable <- df
    
    if (!is.null(sgRNAInfoDatatable) & nrow(sgRNAInfoDatatable) > 0) {
      #make color interval for heatmap
      if(effect){
        brks_smaller <- seq(-3, 0, length.out = 20)
        brks_bigger <- seq(0, 3, length.out = 20)
      }else{
        brks_smaller <- seq(-8, 0, length.out = 20)
        brks_bigger <- seq(0, 8, length.out = 20)
      }
      
      clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
      clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
      {paste0("rgb(", ., ",", ., ",255)")}
      
      brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
      clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
      
      #specify which columns should be frozen and which should have the heatmap
      nfreezeColumns <- 1
      nColorizeTableColumns <- 4
      
      colnames_sgRNAInfoDatatable <- colnames(sgRNAInfoDatatable)
      tooltip <- ''
      
      presel_contrasts <- colnames_sgRNAInfoDatatable[4:(colnames_sgRNAInfoDatatable %>% length)]
      
      for(i in 1:length(colnames_sgRNAInfoDatatable)){
        if(i < length(colnames_sgRNAInfoDatatable)){
          tooltip <- paste0(tooltip, "'", colnames_sgRNAInfoDatatable[i], "'",  ", " )
        }else{
          tooltip <- paste0(tooltip, "'", colnames_sgRNAInfoDatatable[i], "'")
        }
        
        if(colnames_sgRNAInfoDatatable[i] %in% presel_contrasts){
          colnames_sgRNAInfoDatatable[i] <- contrasts %>% select(contrast_id, contrast_id_QC) %>% filter(contrast_id == colnames_sgRNAInfoDatatable[i]) %>% .$contrast_id_QC
          
        }
      }
      
      colnames(sgRNAInfoDatatable) <- colnames_sgRNAInfoDatatable
      
      headerCallback <- c(
        "function(thead, data, start, end, display){",
        "  var $ths = $(thead).find('th');",
        "  $ths.css({'vertical-align': 'bottom', 'white-space': 'nowrap'});",
        "  var betterCells = [];",
        "  $ths.each(function(){",
        "    var cell = $(this);",
        "    var newDiv = $('<div>', {height: 'auto', width: 'auto'});",
        "    var newInnerDiv = $('<div>', {text: cell.text()});",
        "    newDiv.css({margin: 'auto'});",
        "    newInnerDiv.css({",
        "      transform: 'rotate(180deg)',",
        "      'writing-mode': 'tb-rl',",
        "      'white-space': 'nowrap'",
        "    });",
        "    newDiv.append(newInnerDiv);",
        "    betterCells.push(newDiv);",
        "  });",
        "  $ths.each(function(i){",
        "    $(this).html(betterCells[i]);",
        "  });",
        paste0("  var tooltips = [", tooltip, "];"),
        paste0("  for(var i=0; i<", length(colnames_sgRNAInfoDatatable), "; i++){"),
        "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
        "  }",
        "}"
      )
      sgRNAInfoDatatable %>% datatable(escape = FALSE,
                                       selection = 'none', 
                                       extensions = c('FixedColumns','FixedHeader'),
                                       options = list(autoWidth = FALSE, 
                                                      headerCallback = JS(headerCallback), 
                                                      scrollX=TRUE,
                                                      fixedColumns = list(leftColumns = nfreezeColumns), 
                                                      columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                      pageLength = 25, 
                                                      lengthMenu = c(25, 50, 100, 200)),
                                       filter = list(position = 'top', clear = FALSE),
                                       rownames= FALSE) %>%
        formatStyle(seq(nColorizeTableColumns, length(colnames_sgRNAInfoDatatable),1), backgroundColor = styleInterval(brks, clrs))
    }
    
  }else{
    if(!is.null(input$sgRNAInfoSelectGuide)){
      output$sgRNAInfoInfo <- renderText({
        "WARNING: No experimental data found!"
      })
    }else{
      sgRNAInfoUpdateText()
    }
    NULL
  }
})

sgRNAInfoDataTablePredictions <- eventReactive(input$sgRNAInfoLoadButton,{
  df <- sgRNAInfoTablePredictions()
  if (nrow(df) > 0) {
    if(!is.null(input$sgRNAInfoSelectGene)){
      output$sgRNAInfoInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    
    sgRNAInfoDatatable <- df
    
    
    sgRNAInfoDatatable %>% datatable(escape = FALSE,
                                     selection = 'none', 
                                     extensions = c('FixedColumns','FixedHeader'),
                                     options = list(autoWidth = FALSE, 
                                                    headerCallback = JS(headerCallback), 
                                                    scrollX=TRUE,
                                                    columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                    pageLength = 25, 
                                                    lengthMenu = c(25, 50, 100, 200)),
                                     filter = list(position = 'top', clear = FALSE),
                                     rownames= FALSE)
    
    
  }else{
    if(!is.null(input$sgRNAInfoSelectGuide)){
      output$sgRNAInfoInfo <- renderText({
        "WARNING: No sgRNA prediction data found!"
      })
    }else{
      sgRNAInfoUpdateText()
    }
    NULL
  }
})

sgRNAInfoDataTableValidations <- eventReactive(input$sgRNAInfoLoadButton,{
  df <- sgRNAInfoTableValidations()
  if (nrow(df) > 0) {
    if(!is.null(input$sgRNAInfoSelectGene)){
      output$sgRNAInfoInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    
    sgRNAInfoDatatable <- df
    
    
    sgRNAInfoDatatable %>% datatable(escape = FALSE,
                                     selection = 'none', 
                                     extensions = c('FixedColumns','FixedHeader'),
                                     options = list(autoWidth = FALSE, 
                                                    headerCallback = JS(headerCallback), 
                                                    scrollX=TRUE,
                                                    columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                    pageLength = 25, 
                                                    lengthMenu = c(25, 50, 100, 200)),
                                     filter = list(position = 'top', clear = FALSE),
                                     rownames= FALSE)
    
    
  }else{
    if(!is.null(input$sgRNAInfoSelectGuide)){
      output$sgRNAInfoInfo <- renderText({
        "WARNING: No sgRNA prediction data found!"
      })
    }else{
      gwsGeneUpdateText()
    }
    NULL
  }
})

# ----------------------------------------------------------------------------
# Selectbox lists
# ----------------------------------------------------------------------------
sgRNAInfoGeneList <- reactive({
  
  if(class(gene_list_screens)[1] == "tbl_SQLiteConnection" ){
    gene_list_screens <<- gene_list_screens %>%
      collect()
  }
  
  if(input$sgRNAInfoSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$sgRNAInfoSpeciesSelect
  }
  
  libraries_selected <- libraries %>%
      filter(species %in% speciesList) %>%
      .$library_id
    
  gene_list_screens %>%
    filter(library_id %in% libraries_selected) %>%
    select(symbol, entrez_id) %>%
    distinct() %>%
    dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
    arrange(gene) %>%
    .$gene
})

sgRNAInfoGuideList <- reactive({
  if(input$sgRNAInfoSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$sgRNAInfoSpeciesSelect
  }
  
  if((!is.null(input$sgRNAInfoSelectGene))){
    #retrieve selected genes
    presel_genes_both<- input$sgRNAInfoSelectGene %>% strsplit(split="\\(|\\)")
    presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
    presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
    
    features %>%
      filter(entrez_id %in% c(presel_entrez) | symbol %in% c(presel_genes)) %>%
      select(guide_id) %>%
      distinct() %>%
      collect() %>%
      arrange(guide_id) %>%
      .$guide_id
  }
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------

observeEvent(input$sgRNAInfoLoadButton, {
  output$sgRNAInfoTableOutputScreens <- renderDataTable({
    sgRNAInfoDataTableScreens()
  })
  
  output$sgRNAInfoTableOutputPredictions <- renderDataTable({
    sgRNAInfoDataTablePredictions()
  })
  
  output$sgRNAInfoTableOutputValidations <- renderDataTable({
    sgRNAInfoDataTableValidations()
  })
})

observeEvent(input$sgRNAInfoSpeciesSelect, {
  #update gene selectbox
  updateSelectizeInput(session, 'sgRNAInfoSelectGene', choices = sgRNAInfoGeneList(), server = TRUE)
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoSelectGene, {
  #update guide selectbox
  updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = input$sgRNAInfoSelectGuide, server = TRUE)
  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoSelectGuide, {
  if((!is.null(input$sgRNAInfoSelectGuide))){
    enable("sgRNAInfoLoadButton")
    updateCheckboxInput(session, 'sgRNAInfoCheckGuideAll', value = FALSE)
  }else{
    if(!isTRUE(input$sgRNAInfoCheckGuideAll)){
      disable("sgRNAInfoLoadButton")
    }
  }
  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoCheckGuideAll, {
  if((isTRUE(input$sgRNAInfoCheckGuideAll))){
    enable("sgRNAInfoLoadButton")
    updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = NULL, server = TRUE)
  }else{
    if((is.null(input$sgRNAInfoSelectGuide))){
      disable("sgRNAInfoLoadButton")
    }
  }
  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$sgRNAInfoButtonDownload <- downloadHandler(
  filename = function() {
    paste0(local(input$sgRNAInfoSelectGene), ".txt")
  },
  content = function(file) {
    sgRNAInfoTable() %>% write_tsv(file)
  }
)