# ----------------------------------------------------------------------------
# Gene Search
# ----------------------------------------------------------------------------

#necessary for reactive buttons, numerates buttons
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

gwsGeneUpdateText <- function(){
  output$gwsGeneInfo <- renderText({
    if(is.null(input$gwsGeneTissueSelect) & !isTRUE(input$gwsGeneCheckTissueAll)){
      "INFO: Please select the tissue(s) of considered screens in the right panel!"
    }else{
      if(is.null(input$gwsGeneLibrarySelect) & !isTRUE(input$gwsGeneCheckLibraryAll)){
        "INFO: Please select the librarie(s) of considered screens in the right panel!"
      }else{
        if(is.null(input$gwsGeneContrastSelect) & !isTRUE(input$gwsGeneCheckContrastAll)){
          "INFO: Please select the contrasts(s) you want to browse in the right panel!"
        }else{
          if(is.null(input$gwsGeneGeneSelect)){
            "INFO: Please select the gene(s) you want to browse in the right panel!"
          }else{
            "INFO: Click Load data!"
          }
        }
      }
    }
  })
}

#at initial load display nothing
output$gwsGeneTable <- renderDataTable({
})


gwsGeneDataFrame <- reactive({
  #specify which table to select
  if(input$gwsGeneSearchRadio == "guide_id"){
    tableSelect <- "guide_stats"
  }else{
    tableSelect <- "gene_stats"
  }
  
  #specify which sgRNAs table to join
  if(input$gwsGeneSpeciesSelect == "human"){
    tableSgRNAs <- "sgRNAs_human"
  }else{
    tableSgRNAs <- "sgRNAs_mouse"
  }
  
  #retrieve selected genes
  presel_genes_both<- input$gwsGeneGeneSelect %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  
  #retrieve selected contrasts
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    presel_contrasts <- gwsGeneContrastList()
  }else{
    presel_contrasts <- local(input$gwsGeneContrastSelect)
  }
  
  #query database
  if(input$gwsGeneDatasetSelect %in% c("facs")){
    database_con <- con_facs
    if(input$gwsGeneSpeciesSelect == "human"){
      gene_select <- presel_genes
    }else{
      gene_select <- presel_entrez
    }
    features_buff <- features_facs
    contrasts_buff <- contrasts_facs
    
  }else{
    database_con <- con
    gene_select <- presel_entrez
    features_buff <- features
    contrasts_buff <- contrasts
  }
  
  #get gene_stats/guide_stats
  df <- database_con %>%
    tbl(tableSelect) %>%
    filter(gene_id %in% gene_select, contrast_id %in% presel_contrasts) %>%
    left_join(contrasts_buff, copy=TRUE, by = "contrast_id")
  
  if(input$gwsGeneSearchRadio == "guide_id"){
    if(input$gwsGeneSpeciesSelect == "all"){
      sgRNAs <- con_sgRNAs %>%
        tbl("sgRNAs_human") %>%
        filter(EntrezID %in% c(presel_entrez)) %>%
        select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
        rbind(
          con_sgRNAs %>%
            tbl("sgRNAs_mouse") %>%
            filter(EntrezID %in% c(presel_entrez)) %>%
            select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
            collect() %>%
            mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation))
        )
    }else{
      #get guide ranks
      sgRNAs <- con_sgRNAs %>%
        tbl(tableSgRNAs) %>%
        filter(EntrezID %in% c(presel_entrez)) %>%
        select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation))
    }
    
    #join guide ranks
    df <- df %>%
      left_join(features_buff, by = c("gene_id", "guide_id")) %>%
      select(contrast_id, contrast_id_QC, guide_id, gene_id, lfc, effect, symbol, entrez_id, sequence, context, species) %>%
      distinct() %>%
      collect() %>%
      dplyr::mutate(sgRNA_23mer = substr(context, 5, nchar(context)-3)) %>%
      left_join(sgRNAs, copy=TRUE, by=c("sgRNA_23mer" = "sgRNA_23mer", "entrez_id" = "entrez_id")) %>%
      distinct() %>%
      mutate_at(c("lfc","effect"), funs(round(., 3))) %>%
      arrange(symbol)
    
  }else{
    df <- df %>%
      left_join(features_buff %>% select(-guide_id, -sequence) %>% distinct, by ="gene_id") %>%
      select(contrast_id, contrast_id_QC, gene_id, lfc, effect, symbol, entrez_id, species) %>%
      distinct() %>%
      collect() %>%
      mutate_at(c("lfc","effect"), funs(round(., 3))) %>%
      arrange(symbol)
  }
  
  if(input$gwsGeneSpeciesSelect == "all"){
    df_human <- df %>% 
      filter(species == "human") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_human")) %>%
      mutate(EntrezID_human = entrez_id) %>%
      rename(Symbol_human = symbol) %>%
      select(-entrez_id)
    
    df_mouse <- df %>% 
      filter(species == "mouse") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_mouse")) %>%
      mutate(EntrezID_mouse = entrez_id) %>%
      rename(Symbol_mouse = symbol) %>%
      select(-entrez_id)
    
    df <- df_human %>% rbind(df_mouse)
  }
  
  #remove duplicate entries to enable spread
  remove_duplicates <- df[df %>% select(contrast_id, local(input$gwsGeneSearchRadio)) %>% duplicated,]
  df <- df %>%
    anti_join(remove_duplicates)
  
  if (nrow(df) > 0) {
    presel_contrasts <- df$contrast_id %>% unique
    
    df <- df %>%
      select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), local(input$gwsGeneIndexRadio), matches("sequence"), matches("Length"), matches("sgRNA_23mer"), matches("VBC-score"), matches("rank_overall"), matches("rank_validation")) %>%
      spread(contrast_id, input$gwsGeneIndexRadio) %>%
      distinct() %>%
      select(contains("EntrezID_human"), contains("Symbol_human"), contains("EntrezID_mouse"), contains("Symbol_mouse"), everything())
    
    #create actionButtons for each rows
    Action_gene <- shinyInput(actionButton, nrow(df), 'button_gene_', label = HTML("gene <br/> predictions"), onclick = 'Shiny.onInputChange(\"select_button_gene\",  this.id)' )
    df <- df %>% cbind(Action_gene)
    if(input$gwsGeneSearchRadio == "guide_id"){
      Action_sgRNA <- shinyInput(actionButton, nrow(df), 'button_sgRNA_', label = HTML("sgRNA <br/> info"), onclick = 'Shiny.onInputChange(\"select_button_sgRNA\",  this.id)' )
      df <- df %>%
        cbind(Action_sgRNA)
    }
    
    #remove sensitive data for external view
    if(view=="external"){
      df %>%
        select(-contains("rank_overall"), -contains("rank_validation"), -contains("Action_sgRNA"), -contains("Action_gene"))
    }else{
      df
    }
  }
})

#make 
gwsGeneDataTable <- eventReactive(input$gwsGeneLoadButton,{
  df <- gwsGeneDataFrame()
  
  if (nrow(df) > 0) {
    if(!is.null(input$gwsGeneGeneSelect)){
      output$gwsGeneInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    
    effect<-FALSE
    if(input$gwsGeneIndexRadio == "effect"){
      effect<-TRUE
    }
    
    gwsGeneDatatable <- df
    
    if (!is.null(gwsGeneDatatable) & nrow(gwsGeneDatatable) > 0) {
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
      if(!is.null(gwsGeneDatatable$sequence)){
        nfreezeColumns <- 3
        nColorizeTableColumns <- 9
        nActionButtons <- 2
      }else{
        nfreezeColumns <- 2
        nColorizeTableColumns <- 3
        nActionButtons <- 1
      }
      
      if(input$gwsGeneSpeciesSelect == "all"){
        nfreezeColumns <- nfreezeColumns + 2
        nColorizeTableColumns <- nColorizeTableColumns + 2
      }
      
      if(input$gwsGeneSearchRadio == "guide_id"){
        presel_contrasts <- colnames(gwsGeneDatatable)[nColorizeTableColumns:(length(colnames(gwsGeneDatatable))-2)]
      }else{
        presel_contrasts <- colnames(gwsGeneDatatable)[nColorizeTableColumns:(length(colnames(gwsGeneDatatable))-1)]
      }
      
      colnames_gwsGeneDatatable <- colnames(gwsGeneDatatable)
      tooltip <- ''
      if(input$gwsGeneDatasetSelect %in% c("dropout")){
        for(i in 1:length(colnames_gwsGeneDatatable)){
          if(i < length(colnames_gwsGeneDatatable)){
            tooltip <- paste0(tooltip, "'", colnames_gwsGeneDatatable[i], "'",  ", " )
          }else{
            tooltip <- paste0(tooltip, "'", colnames_gwsGeneDatatable[i], "'")
          }
          
          if(colnames_gwsGeneDatatable[i] %in% presel_contrasts){
            colnames_gwsGeneDatatable[i] <- contrasts %>% select(contrast_id, contrast_id_QC) %>% filter(contrast_id == colnames_gwsGeneDatatable[i]) %>% .$contrast_id_QC
            
          }
        }
        
        colnames(gwsGeneDatatable) <- colnames_gwsGeneDatatable
      }
      
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
        paste0("  for(var i=0; i<", length(colnames_gwsGeneDatatable), "; i++){"),
        "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
        "  }",
        "}"
      )
      gwsGeneDatatable %>% datatable(escape = FALSE,
                                     selection = 'none', 
                                     extensions = c('FixedColumns','FixedHeader'),
                                     options = list(autoWidth = FALSE, 
                                                    headerCallback = JS(headerCallback), 
                                                    scrollX=TRUE,
                                                    fixedColumns = list(leftColumns = nfreezeColumns), 
                                                    columnDefs = list(list(className = 'dt-center', targets = "_all")), 
                                                    pageLength = 25, 
                                                    lengthMenu = c(25, 50, 100, 200),
                                                    searchHighlight = TRUE),
                                     filter = list(position = 'top', clear = FALSE),
                                     rownames= FALSE) %>%
        formatStyle(seq(nColorizeTableColumns, length(colnames_gwsGeneDatatable),1), backgroundColor = styleInterval(brks, clrs))
    }
    
  }else{
    if(!is.null(input$gwsGeneGeneSelect)){
      output$gwsGeneInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      gwsGeneUpdateText()
    }
  }
})

# ----------------------------------------------------------------------------
# Selectbox lists
# ----------------------------------------------------------------------------

gwsGeneTissueList <- reactive({
  if(class(pheno)[1] == "tbl_SQLiteConnection"){
    pheno <<- pheno %>%
      collect()
    
    libraries <<- libraries %>%
      collect()
    
    contrasts <<- contrasts %>%
      collect()
    
    contrasts_facs <<- contrasts_facs %>%
      collect()
    
    gene_list_screens <<- gene_list_screens %>%
        collect()
  }
  
  if(input$gwsGeneSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGeneSpeciesSelect
  }
  
  pheno %>%
    filter(species %in% speciesList) %>%
    select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
})

gwsGeneLibraryList <- reactive({
  if(input$gwsGeneSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGeneSpeciesSelect
  }
  
  libraries %>%
    filter(species %in% speciesList) %>%
    select(library_id) %>%
    arrange(library_id) %>%
    .$library_id
})

gwsGeneContrastList <- reactive({
  
  if(!is.null(input$gwsGeneLibrarySelect) | isTRUE(input$gwsGeneCheckLibraryAll)){
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    preselLibrary <- libraries %>%
      filter(species %in% speciesList) %>%
      select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
      preselLibrary = libraries %>%
        filter(species %in% speciesList) %>%
        filter(library_id %in% input$gwsGeneLibrarySelect) %>%
        .$library_id
    }
    
    preselTissue <- pheno %>%
      filter(species %in%  speciesList) %>%
      select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue = pheno %>%
        filter(species %in%  speciesList) %>%
        filter(tissue_name %in% input$gwsGeneTissueSelect) %>%
        .$tissue_name
    }
    
    if(input$gwsGeneDatasetSelect %in% c("facs")){
      contrasts_buff <- contrasts_facs
      if("human" %in%  speciesList){
        preselLibrary <- "zuber_library_original"
      }else{
        preselLibrary <- "zuber_library_mouse_original"
      }
    }else{
      contrasts_buff <- contrasts
    }
    contrasts_buff %>%
      filter(species %in%  speciesList) %>%
      filter(library_id %in% preselLibrary) %>%
      filter(tissue_name %in%  preselTissue) %>%
      filter(type == input$gwsGeneDatasetSelect) %>%
      select(contrast_id) %>%
      distinct() %>%
      .$contrast_id
  }
  
})

gwsGeneGeneList <- reactive({
  
  if(!is.null(input$gwsGeneContrastSelect) | isTRUE(input$gwsGeneCheckContrastAll)){
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    if(input$gwsGeneDatasetSelect %in% c("facs")){
      if(isTRUE(input$gwsGeneCheckContrastAll)){
        presel_contrasts <- gwsGeneContrastList()
      }else{
        presel_contrasts <- local(input$gwsGeneContrastSelect)
      }
      
      con_facs %>%
        tbl("gene_stats") %>%
        filter(contrast_id %in% presel_contrasts) %>%
        select(gene_id) %>%
        distinct() %>%
        left_join(con_facs %>% tbl("features") %>% select(gene_id, symbol, entrez_id)) %>%
        collect() %>%
        dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
        arrange(gene) %>%
        .$gene
      
    }else{
      preselLibrary <-  libraries %>%
        filter(species %in% speciesList) %>%
        select(library_id) %>%
        .$library_id
      
      if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
        preselLibrary = libraries %>%
          filter(library_id %in% input$gwsGeneLibrarySelect) %>%
          collect() %>%
          .$library_id
      }
      gene_list_screens %>%
        filter(library_id %in% preselLibrary) %>% 
        dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
        arrange(gene) %>%
        .$gene
    }
  }
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------
observeEvent(input$gwsGeneLoadButton, {
  output$gwsGeneTable <- renderDataTable({
    gwsGeneDatatable <- gwsGeneDataTable()
  })
})

#actionButton handler for datatable buttons
observeEvent(input$select_button_gene, {
  row <- as.numeric(strsplit(input$select_button_gene, "_")[[1]][3])
  if(input$gwsGeneSpeciesSelect == "all"){
    symbol_human <- gwsGeneDataFrame()[row,2]
    entrez_id_human <- gwsGeneDataFrame()[row,1]
    symbol_mouse <- gwsGeneDataFrame()[row,4]
    entrez_id_mouse <- gwsGeneDataFrame()[row,3]
    gene = c(ifelse(is.na(symbol_human), paste0("No symbol found (", entrez_id_human, ")"), paste0(symbol_human , " (", entrez_id_human, ")")),
             ifelse(is.na(symbol_mouse), paste0("No symbol found (", entrez_id_mouse, ")"), paste0(symbol_mouse , " (", entrez_id_mouse, ")")))
  }else{
    symbol <- gwsGeneDataFrame()[row,2]
    entrez_id <- gwsGeneDataFrame()[row,1]
    gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))
    
  }

  delay(250, updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), selected = gene, server = TRUE))
  updateTabItems(session, "tabs", "sgRNAsSidebar")
  enable("sgRNAsLoadButton")
  delay(500, click("sgRNAsLoadButton"))
})

#actionButton handler for datatable buttons
observeEvent(input$select_button_sgRNA, {
  row <- as.numeric(strsplit(input$select_button_sgRNA, "_")[[1]][3])
  if(input$gwsGeneSpeciesSelect == "all"){
    symbol_human <- gwsGeneDataFrame()[row,2]
    entrez_id_human <- gwsGeneDataFrame()[row,1]
    symbol_mouse <- gwsGeneDataFrame()[row,4]
    entrez_id_mouse <- gwsGeneDataFrame()[row,3]
    gene <- c(ifelse(is.na(symbol_human), paste0("No symbol found (", entrez_id_human, ")"), paste0(symbol_human , " (", entrez_id_human, ")")),
              ifelse(is.na(symbol_mouse), paste0("No symbol found (", entrez_id_mouse, ")"), paste0(symbol_mouse , " (", entrez_id_mouse, ")")))
    sequence <- gwsGeneDataFrame()[row,5]
    guide_id<-c(paste0(entrez_id_human, "_", "sequence"), paste0(entrez_id_human, "_", sequence))
  }else{
    symbol <- gwsGeneDataFrame()[row,2]
    entrez_id <- gwsGeneDataFrame()[row,1]
    gene <- ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))
    sequence <- gwsGeneDataFrame()[row,3]
    guide_id<-paste0(entrez_id, "_", sequence)
  }

  updateSelectizeInput(session, 'sgRNAInfoSelectGene', choices = sgRNAInfoGeneList(), selected = gene, server = TRUE)
  delay(250, updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = guide_id, server = TRUE))
  updateTabItems(session, "tabs", "sgRNAInfoSidebar")
  enable("sgRNAInfoLoadButton")
  delay(1500, click("sgRNAInfoLoadButton"))
})

observeEvent(input$gwsGeneSpeciesSelect, {
  if(input$gwsGeneSpeciesSelect == "all"){
    updateSelectizeInput(session, 'gwsGeneDatasetSelect', choices = dataset_selection_dropout_drug, selected = "dropout", server = TRUE)
  }else{
    updateSelectizeInput(session, 'gwsGeneDatasetSelect', choices = dataset_selection_all, selected = "dropout", server = TRUE)
  }
  #select checkbox tissue
  updateCheckboxInput(session, 'gwsGeneCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'gwsGeneTissueSelect', choices = gwsGeneTissueList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  #update library selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  #update other species selects
  updateSelectizeInput(session, 'sgRNAsSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = input$gwsGeneSpeciesSelect, server = TRUE)
  updateSelectizeInput(session, 'gwsBrowseScreenSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = input$gwsGeneSpeciesSelect, server = TRUE)
  #disable load button
  disable("gwsGeneLoadButton")
  gwsGeneUpdateText()
})

observeEvent(input$gwsGeneDatasetSelect, {
  #select checkbox tissue
  updateCheckboxInput(session, 'gwsGeneCheckTissueAll', value = FALSE)
  #update tissue selectbox
  #update tissue selectbox
  updateSelectizeInput(session, 'gwsGeneTissueSelect', choices = gwsGeneTissueList(), server = TRUE)
  if(input$gwsGeneDatasetSelect %in% c("synthetic", "facs")){
    if("human" ==  input$gwsGeneSpeciesSelect){
      #update library selectbox
      updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), selected = "hs_gw_zuber_v2", server = TRUE)
    }else{
      #update library selectbox
      updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), selected = "mm_gw_zuber_v1", server = TRUE)
    }
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
    disable("gwsGeneLibrarySelect")
    disable("gwsGeneCheckLibraryAll")
    #disable index
    updateRadioButtons(session, 'gwsGeneIndexRadio',selected = "lfc")
    disable("gwsGeneIndexRadio")
  }else{
    #update library selectbox
    updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
    #select library checkbox
    updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
    enable("gwsGeneIndexRadio")
  }
  #update contrasts select
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneTissueSelect, {
  if(!is.null(input$gwsGeneTissueSelect)){
    #reset tissue selectbox
    updateCheckboxInput(session, 'gwsGeneCheckTissueAll', value = FALSE)
  }
  if(input$gwsGeneDatasetSelect %in% c("synthetic", "facs")){
    if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
      enable("gwsGeneContrastSelect")
      enable("gwsGeneCheckContrastAll")
    }else{
      disable("gwsGeneContrastSelect")
      disable("gwsGeneCheckContrastAll")
    }
  }else{
    #select library checkbox
    updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
    #update library selectbox
    updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
    if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
      enable("gwsGeneLibrarySelect")
      enable("gwsGeneCheckLibraryAll")
    }else{
      disable("gwsGeneLibrarySelect")
      disable("gwsGeneCheckLibraryAll")
    }
  }
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCheckTissueAll, {
  if(isTRUE(input$gwsGeneCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsGeneTissueSelect', choices = gwsGeneTissueList(), server = TRUE)
  }
  if(input$gwsGeneDatasetSelect %in% c("synthetic", "facs")){
    if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
      enable("gwsGeneContrastSelect")
      enable("gwsGeneCheckContrastAll")
    }else{
      disable("gwsGeneContrastSelect")
      disable("gwsGeneCheckContrastAll")
    }
  }else{
    #select library checkbox
    updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
    #update library selectbox
    updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
    if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
      enable("gwsGeneLibrarySelect")
      enable("gwsGeneCheckLibraryAll")
    }else{
      disable("gwsGeneLibrarySelect")
      disable("gwsGeneCheckLibraryAll")
    }
  }
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  gwsGeneUpdateText()
  
})

observeEvent(input$gwsGeneLibrarySelect, {
  if(!is.null(input$gwsGeneLibrarySelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  }
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  
  if((isTRUE(input$gwsGeneCheckLibraryAll) | (!is.null(input$gwsGeneLibrarySelect))) & (isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect)))){
    enable("gwsGeneContrastSelect")
    enable("gwsGeneCheckContrastAll")
  }else{
    disable("gwsGeneContrastSelect")
    disable("gwsGeneCheckContrastAll")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCheckLibraryAll, {
  if(isTRUE(input$gwsGeneCheckLibraryAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  }
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  
  if((isTRUE(input$gwsGeneCheckLibraryAll) | (!is.null(input$gwsGeneLibrarySelect))) & (isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect)))){
    enable("gwsGeneContrastSelect")
    enable("gwsGeneCheckContrastAll")
  }else{
    disable("gwsGeneContrastSelect")
    disable("gwsGeneCheckContrastAll")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)



observeEvent(input$gwsGeneContrastSelect, {
  if(!is.null(input$gwsGeneContrastSelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  }
  
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  
  if(isTRUE(input$gwsGeneCheckContrastAll) | (!is.null(input$gwsGeneContrastSelect))){
    enable("gwsGeneGeneSelect")
  }else{
    disable("gwsGeneGeneSelect")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCheckContrastAll, {
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    contrastList <- gwsGeneContrastList()
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = contrastList, server = TRUE)
    
    showModal(modalDialog(
      title = "WARNING!", 
      paste0("WARNING: You have selected ", 
             length(contrastList), 
             " contrasts. Are you sure you want to load all contrasts for this selection? This might take a long time to load!"),
      footer = tagList(
        modalButton("OK"),
        actionButton("gwsGeneCancelModal", "Cancel")
      )
    ))
  }
  
  #update gene selectbox
  updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
  
  if(isTRUE(input$gwsGeneCheckContrastAll) | (!is.null(input$gwsGeneContrastSelect))){
    enable("gwsGeneGeneSelect")
  }else{
    disable("gwsGeneGeneSelect")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCancelModal, {
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  removeModal()
})

observeEvent(input$gwsGeneGeneSelect, {
  if((!is.null(input$gwsGeneGeneSelect))){
    enable("gwsGeneLoadButton")
  }else{
    disable("gwsGeneLoadButton")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$gwsGeneButtonDownload <- downloadHandler(
  filename = function() {
    table <- gwsGeneDataFrame()
    paste0(paste(table$symbol %>% unique,collapse="_"), ".txt")
  },
  content = function(file) {
    table <- gwsGeneDataFrame()
    table %>% select(-contains("Action_sgRNA"), -contains("Action_gene")) %>%
      write_tsv(file)
  }
)