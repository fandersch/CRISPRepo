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
  
  #get dropout screens
  if(input$sgRNAInfoSpeciesSelect == "all"){
    contrasts_dropout <- contrasts %>%
      filter(type=="dropout", reference_type=="reference")
  }else{
    contrasts_dropout <- contrasts %>%
      filter(species == local(input$sgRNAInfoSpeciesSelect), type=="dropout", reference_type=="reference")
  }
  
  #retrieve selected genes
  presel_genes_both<- input$sgRNAInfoSelectGene %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
  if(isTRUE(input$sgRNAInfoCheckGuideAll)){
    presel_guides <- sgRNAInfoGuideList()
    
    #get gene_stats/guide_stats
    df_screens <- con %>%
      tbl("guide_stats") %>%
      dplyr::select(contrast_id, guide_id, id_entrez_23mer, gene_id, local(input$sgRNAInfoIndexRadio)) %>%
      dplyr::filter(contrast_id %in% local(contrasts_dropout$contrast_id), gene_id %in% local(presel_entrez) | gene_id %in% local(presel_genes)) %>%
      collect() %>%
      dplyr::left_join(contrasts_dropout %>% dplyr::select(contrast_id, library_id)) %>%
      dplyr::left_join(features %>% select) %>%
      # left_join(gene_list_screens %>% dplyr::select(gene_id, symbol, entrez_id) %>% distinct) %>%
      mutate_at(local(input$sgRNAInfoIndexRadio), round, 3)
    
  }else{
    presel_guides <- input$sgRNAInfoSelectGuide
    
    #get gene_stats/guide_stats
    df_screens <- con %>%
      tbl("guide_stats") %>%
      dplyr::select(contrast_id, guide_id, id_entrez_23mer, gene_id, local(input$sgRNAInfoIndexRadio)) %>%
      dplyr::filter(contrast_id %in% local(contrasts_dropout$contrast_id), id_entrez_23mer %in% local(presel_guides)) %>%
      collect() %>%
      dplyr::left_join(contrasts_dropout %>% dplyr::select(contrast_id, library_id)) %>%
      left_join(features) %>%
      # left_join(gene_list_screens %>% dplyr::select(gene_id, symbol, entrez_id) %>% distinct) %>%
      mutate_at(local(input$sgRNAInfoIndexRadio), round, 3)
    
  }
  
  DBI::dbDisconnect(con)
  
  if (nrow(df_screens) > 0) {
    presel_contrasts <- df_screens$contrast_id %>% unique

    df_screens <- df_screens %>%
      dplyr::select(contrast_id, guide_id, entrez_id, symbol, sequence, sequence_matching, id_entrez_23mer, library_id, local(input$sgRNAInfoIndexRadio)) %>%
      arrange(contrast_id) %>%
      group_by(guide_id) %>%
      mutate(library_id = paste0(unique(library_id), collapse=", ")) %>%
      pivot_wider(names_from=contrast_id, values_from=input$sgRNAInfoIndexRadio) %>%
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
    dplyr::filter(id_entrez_23mer %in% presel_guides) %>%
    separate(id_entrez_23mer, into =c("dummy", "sgRNA_23mer"), sep = "_", remove = FALSE) %>%
    .$sgRNA_23mer %>%
    unique
  
  sgRNA_23mer_search_string <- paste(sgRNAs_23mer, collapse = "|")
  
  con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")
  
  if(input$sgRNAInfoSpeciesSelect == "all"){
    sgRNAs <- con_sgRNAs %>%
      tbl("sgRNAs_human") %>%
      dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
      dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      dplyr::filter(str_detect(`20-mer + NGG`, sgRNA_23mer_search_string)) %>%
      rbind(
        con_sgRNAs %>%
          tbl("sgRNAs_mouse") %>%
          dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
          dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
                 `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
                 Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
                 `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
                 `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
          collect() %>%
          dplyr::filter(str_detect(`20-mer + NGG`, sgRNA_23mer_search_string))
      )
  }else{
    #get guide ranks
    sgRNAs <- con_sgRNAs %>%
      tbl(tableSgRNAs) %>%
      dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
      dplyr::select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
             `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
             Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
             `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
             `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
      collect() %>%
      dplyr::filter(str_detect(`20-mer + NGG`, sgRNA_23mer_search_string))
  }
  
  DBI::dbDisconnect(con_sgRNAs)
  
  if(sgRNAs$`Maps to genome` %>% as.character %>% unique %>% length == "1"){
    sgRNAs <- sgRNAs %>% dplyr::select(-`Maps to genome`)
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
    dplyr::filter(id_entrez_23mer %in% presel_guides) %>%
    separate(id_entrez_23mer, into =c("dummy", "sgRNA_23mer"), sep = "_", remove = FALSE) %>%
    .$sgRNA_23mer %>%
    unique
  
  sgRNA_23mer_search_string <- paste(sgRNAs_23mer, collapse = "|")
  
  con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")
  
  for(i in 1:length(tableSgRNAs)){
    sgRNAs_buff <- con_sgRNAs %>%
      tbl(tableSgRNAs[i]) %>%
      dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
      collect() %>%
      mutate(Dataset = unlist(strsplit(tableSgRNAs[i], split="_"))[3:4] %>% paste(collapse = "_")) %>%
      separate(sgRNA_id, into =c("dummy", "seq_buff"), sep = "_", remove = FALSE) %>%
      dplyr::filter(str_detect(seq_buff, sgRNA_23mer_search_string)) %>%
      dplyr::select(-seq_buff, -dummy, -legacy_sequence, Library_origin = source, Mature_sgRNA = mature_sgRNA) %>%
      dplyr::select(Dataset, everything())
    
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
  
  DBI::dbDisconnect(con_sgRNAs)
  
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
    
    sgRNAInfoDatatable <- df
    
    if (!is.null(sgRNAInfoDatatable) & nrow(sgRNAInfoDatatable) > 0) {
      #make color interval for heatmap
      if(input$sgRNAInfoIndexRadio == "effect_essentialome"){
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
      nColorizeTableColumns <- 8
      
      colnames_sgRNAInfoDatatable <- colnames(sgRNAInfoDatatable)
      tooltip <- ''
      
      presel_contrasts <- colnames_sgRNAInfoDatatable[nColorizeTableColumns:(colnames_sgRNAInfoDatatable %>% length)]
      
      for(i in 1:length(colnames_sgRNAInfoDatatable)){
        if(i < length(colnames_sgRNAInfoDatatable)){
          tooltip <- paste0(tooltip, "'", colnames_sgRNAInfoDatatable[i], "'",  ", " )
        }else{
          tooltip <- paste0(tooltip, "'", colnames_sgRNAInfoDatatable[i], "'")
        }

        if(colnames_sgRNAInfoDatatable[i] %in% presel_contrasts){
          buff<-contrasts %>%
            dplyr::select(contrast_id, contrast_id_QC) %>%
            dplyr::filter(contrast_id == colnames_sgRNAInfoDatatable[i]) %>%
            .$contrast_id_QC
          
          colnames_sgRNAInfoDatatable[i] <- contrasts %>%
            dplyr::select(contrast_id, contrast_id_QC) %>%
            dplyr::filter(contrast_id == colnames_sgRNAInfoDatatable[i]) %>%
            .$contrast_id_QC
        }
      }
      
      colnames(sgRNAInfoDatatable) <- colnames_sgRNAInfoDatatable
      
      headerCallbackScreens <- c(
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
                                                      headerCallback = JS(headerCallbackScreens), 
                                                      scrollX=TRUE,
                                                      fixedColumns = list(leftColumns = nfreezeColumns), 
                                                      columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                      pageLength = 25, 
                                                      lengthMenu = c(25, 50, 100, 200),
                                                      searchHighlight = TRUE),
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
                                                    lengthMenu = c(25, 50, 100, 200),
                                                    searchHighlight = TRUE),
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
    
    colnames_sgRNAInfoDatatableValidation <- colnames(sgRNAInfoDatatable)
    tooltip <- ''
    
    for(i in 1:length(colnames_sgRNAInfoDatatableValidation)){
      if(colnames_sgRNAInfoDatatableValidation[i] == "gene_hits_total"){
        tooltip <- paste0(tooltip, "'The number of screens where the gene was considered as a hit (depletion/enrichment)'",  ", " )
      }else{
        if(colnames_sgRNAInfoDatatableValidation[i] == "good_guide_total"){
          tooltip <- paste0(tooltip, "'The number of screens where this sgRNA showed a stronger depletion/enrichment than the average of all sgRNAs that target this gene'",  ", " )
        }else{
          if(colnames_sgRNAInfoDatatableValidation[i] == "outlier_total"){
            tooltip <- paste0(tooltip, "'The number of screens where this sgRNA was classified as outlier'",  ", " )
          }else{
            if(colnames_sgRNAInfoDatatableValidation[i] == "good_guides_ratio"){
              tooltip <- paste0(tooltip, "'The number of good_guide_total divided by the number of gene_hits_total'",  ", " )
            }else{
              if(colnames_sgRNAInfoDatatableValidation[i] == "avgGuideLFCGeneHit"){
                tooltip <- paste0(tooltip, "'The average LFC of this sgRNA in screens where the gene was considered as a hit (depletion/enrichment)'",  ", " )
              }else{
                if(colnames_sgRNAInfoDatatableValidation[i] == "ratio_avgLFC_rank"){
                  tooltip <- paste0(tooltip, "'good_guides_ratio multiplied by avgGuideLFCGeneHit'",  ", " )
                }else{
                  if(colnames_sgRNAInfoDatatableValidation[i] == "outlier_high_confident_all_total"){
                    tooltip <- paste0(tooltip, "'The number of screens where this sgRNA was classified as outlier high confident'")
                  }else{
                    tooltip <- paste0(tooltip, "'", colnames_sgRNAInfoDatatableValidation[i], "'",  ", " )
                  }
                }
              }
            }
          }
        }
      }
    }
    
    headerCallbackValidation <- c(
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
      paste0("  for(var i=0; i<", length(colnames_sgRNAInfoDatatableValidation), "; i++){"),
      "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
      "  }",
      "}"
    )
    
    sgRNAInfoDatatable %>% datatable(escape = FALSE,
                                     selection = 'none', 
                                     extensions = c('FixedColumns','FixedHeader'),
                                     options = list(autoWidth = FALSE, 
                                                    headerCallback = JS(headerCallbackValidation), 
                                                    scrollX=TRUE,
                                                    columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                    pageLength = 25, 
                                                    lengthMenu = c(25, 50, 100, 200),
                                                    searchHighlight = TRUE),
                                     filter = list(position = 'top', clear = FALSE),
                                     rownames= FALSE)
    
    
  }else{
    if(!is.null(input$sgRNAInfoSelectGuide)){
      output$sgRNAInfoInfo <- renderText({
        "WARNING: No sgRNA validation data found!"
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
  
  if(input$sgRNAInfoSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$sgRNAInfoSpeciesSelect
  }
  
  libraries_selected <- libraries %>%
      dplyr::filter(species %in% speciesList) %>%
      .$library_id
    
  gene_list_screens %>%
    dplyr::filter(library_id %in% libraries_selected) %>%
    dplyr::select(symbol, entrez_id) %>%
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
      dplyr::filter(entrez_id %in% c(presel_entrez) | symbol %in% c(presel_genes)) %>%
      arrange(id_entrez_23mer) %>%
      .$id_entrez_23mer %>%
      unique
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
  #update other species selects
  updateSpecies(input$sgRNAInfoSpeciesSelect)
  #update gene selectbox
  updateSelectizeInput(session, 'sgRNAInfoSelectGene', choices = sgRNAInfoGeneList(), server = TRUE)
  
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoSelectGene, {
  
  updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = input$sgRNAInfoSelectGuide, server = TRUE)
  updateCheckboxInput(session, 'sgRNAInfoCheckGuideAll', value = FALSE)
  
  if((!is.null(input$sgRNAInfoSelectGene))){
    enable("sgRNAInfoCheckGuideAll")
    enable("sgRNAInfoSelectGuide")
  }else{
    disable("sgRNAInfoCheckGuideAll")
    disable("sgRNAInfoSelectGuide")
  }
  
  disable("sgRNAInfoLoadButton")

  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoSelectGuide, {
  if((!is.null(input$sgRNAInfoSelectGuide)) | isTRUE(input$sgRNAInfoCheckGuideAll)){
    enable("sgRNAInfoLoadButton")
  }else{
    disable("sgRNAInfoLoadButton")
  }
  
  if(!is.null(input$sgRNAInfoSelectGuide)){
    updateCheckboxInput(session, 'sgRNAInfoCheckGuideAll', value = FALSE)
  }
  
  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

observeEvent(input$sgRNAInfoCheckGuideAll, {
  if(isTRUE(input$sgRNAInfoCheckGuideAll) | !is.null(input$sgRNAInfoSelectGuide)){
    enable("sgRNAInfoLoadButton")
  }else{
    disable("sgRNAInfoLoadButton")
  }
  
  if(isTRUE(input$sgRNAInfoCheckGuideAll)){
    updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = NULL, server = TRUE)
  }

  sgRNAInfoUpdateText()
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$sgRNAInfoButtonDownload <- downloadHandler(
  
  filename = function() {
    "crisprepo_sgRNA_info.zip"
  },
  content = function(file){
    #go to a temp dir to avoid permission issues
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    files <- NULL;
    #write each sheet to a csv file, save the name
    table <- sgRNAInfoTableScreens()
    fileName <- "crisprepo_sgRNA_data.txt"
    write.table(table,fileName, row.names = F, col.names = T)
    files <- c(fileName,files)
    
    table <- sgRNAInfoTablePredictions()
    fileName <- "crisprepo_sgRNA_prediction_scores.txt"
    write.table(table,fileName, row.names = F, col.names = T)
    files <- c(fileName,files)
    
    table <- sgRNAInfoTableValidations()
    fileName <- "crisprepo_sgRNA_validation_scores.txt"
    write.table(table,fileName, row.names = F, col.names = T)
    files <- c(fileName,files)
    #create the zip file
    zip(file,files, compression_level = 2)
  }
  
)