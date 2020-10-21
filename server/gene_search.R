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
  
  statistics_columns_negative <- paste0(local(input$gwsGeneInclude), "_negative")
  statistics_columns_positive <- paste0(local(input$gwsGeneInclude), "_positive")
  
  database_con <- con
  gene_select <- presel_entrez
  
  #get gene_stats/guide_stats
  df <- con %>%
    tbl(tableSelect) %>%
    filter(gene_id %in% gene_select, contrast_id %in% presel_contrasts) %>%
    left_join(contrasts, copy=TRUE, by = "contrast_id")
  
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
      left_join(features, by = c("gene_id", "guide_id")) %>%
      select(contrast_id, contrast_id_QC, guide_id, gene_id, lfc, effect, symbol, entrez_id, sequence, guide_id, context, species) %>%
      distinct() %>%
      collect() %>%
      # dplyr::mutate(sgRNA_23mer = substr(context, 5, nchar(context)-3)) %>%
      left_join(sgRNAs, copy=TRUE, by=c("guide_id" = "sgRNA_23mer", "entrez_id" = "entrez_id")) %>%
      distinct() %>%
      mutate_at(c("lfc","effect"), funs(round(., 3))) %>%
      arrange(symbol)
    
  }else{
    
    df <- df %>%
      left_join(features %>% select(-guide_id, -sequence) %>% distinct, by ="gene_id") %>%
      select(contrast_id, contrast_id_QC, gene_id, lfc, effect, symbol, entrez_id, species, dplyr::one_of(statistics_columns_negative), dplyr::one_of(statistics_columns_positive)) %>%
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
    
    dt <- df %>%
      select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), 
             contains("Symbol_mouse"), contains("EntrezID_mouse"), local(input$gwsGeneIndexRadio), matches("sequence"), 
             matches("Length"), matches("sgRNA_23mer"), matches("VBC-score"), matches("rank_overall"), matches("rank_validation")) %>%
      spread(contrast_id, input$gwsGeneIndexRadio) %>%
      distinct() %>%
      select(contains("EntrezID_human"), contains("Symbol_human"), contains("EntrezID_mouse"), contains("Symbol_mouse"), everything())
    
    if(input$gwsGeneSearchRadio == "gene_id" & !is.null(local(input$gwsGeneInclude)) & length(local(input$gwsGeneInclude))!=0){
      column <- local(input$gwsGeneInclude)
      if("p"  %in% local(input$gwsGeneInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(p = ifelse(p_positive < p_negative, p_positive, p_negative)) %>%
              select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), p) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="p") %>%
              mutate_if(is.numeric, round, 3) %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_P'))
          )
      }
      if("fdr" %in% local(input$gwsGeneInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(fdr = ifelse(fdr_positive < fdr_negative, fdr_positive, fdr_negative)) %>%
              select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), fdr) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="fdr") %>%
              mutate_if(is.numeric, round, 3) %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_FDR'))
          )
      }
      if("guides" %in% local(input$gwsGeneInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(guides = ifelse(guides_positive > guides_negative, guides_positive, guides_negative)) %>%
              select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), guides) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="guides") %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_GUIDES'))
          )
      }
      if("guides_good" %in% local(input$gwsGeneInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(guides_good = ifelse(guides_good_positive > guides_good_negative, guides_good_positive, guides_good_negative)) %>%
              select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), guides_good) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="guides_good") %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_GUIDES_GOOD'))
          )
      }
    }
    
    dt <- dt %>%
      select(sort(current_vars())) %>%
      select(contains(local(input$gwsGeneSearchRadio)), contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), everything())
    
    #create actionButtons for each rows
    Action_gene <- shinyInput(actionButton, nrow(dt), 'button_gene_', label = HTML("gene <br/> predictions"), onclick = 'Shiny.onInputChange(\"select_button_gene\",  this.id)' )
    dt <- dt %>% cbind(Action_gene)
    if(input$gwsGeneSearchRadio == "guide_id"){
      Action_sgRNA <- shinyInput(actionButton, nrow(dt), 'button_sgRNA_', label = HTML("sgRNA <br/> info"), onclick = 'Shiny.onInputChange(\"select_button_sgRNA\",  this.id)' )
      dt <- dt %>%
        cbind(Action_sgRNA)
    }
    
    #remove sensitive data for external view
    if(view=="external"){
      dt %>%
        select(-contains("rank_overall"), -contains("rank_validation"), -contains("Action_sgRNA"), -contains("Action_gene"))
    }else{
      dt
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
      #set colorized column, depending on shown gene-level statistics
      colorInterval <- length(local(input$gwsGeneInclude))
      #set frozen/colourized columns
      if(input$gwsGeneSpeciesSelect == "all"){
        nfreezeColumns <- nfreezeColumns + 2
        nColorizeTableColumns <- nColorizeTableColumns + 2
      }
      
      #change to contrast_id_QC for selected columns when type is dropout
      presel_contrasts <- gwsGeneContrastList()
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
        formatStyle(seq(nColorizeTableColumns, length(colnames_gwsGeneDatatable),1+colorInterval), backgroundColor = styleInterval(brks, clrs))
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
    
    gene_list_screens <<- gene_list_screens %>%
        collect()
  }
  
  if(input$gwsGeneSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGeneSpeciesSelect
  }
  
  pheno %>%
    filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect) %>%
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
  
  if(input$gwsGeneCheckTissueAll == T){
    libraries %>%
      filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      select(library_id) %>%
      arrange(library_id) %>%
      .$library_id
  }else{
    libraries %>%
      filter(species %in% speciesList, tissue_name %in% input$gwsGeneTissueSelect, type %in% input$gwsGeneDatasetSelect) %>%
      select(library_id) %>%
      arrange(library_id) %>%
      .$library_id
  }
})

gwsGeneContrastList <- reactive({
  
  if(!is.null(input$gwsGeneLibrarySelect) | isTRUE(input$gwsGeneCheckLibraryAll)){
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    #get selected tissue
    if(isTRUE(input$gwsGeneCheckTissueAll)){
      preselTissue <- pheno %>%
        filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
        select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }else{
      if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
        preselTissue <- pheno %>%
          filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
          select(tissue_name) %>%
          distinct() %>%
          .$tissue_name
      }
    }
    
    #get selected library
    if(isTRUE(input$gwsGeneCheckLibraryAll)){
      preselLibrary <- libraries %>%
        filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
        select(library_id) %>%
        .$library_id
    }else{
      if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
        preselLibrary <- libraries %>%
          filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, library_id %in% input$gwsGeneLibrarySelect) %>%
          select(library_id) %>%
          .$library_id
      }
    }

    contrasts %>%
      filter(species %in%  speciesList, library_id %in% preselLibrary, tissue_name %in%  preselTissue, type == input$gwsGeneDatasetSelect) %>%
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
    
    if(isTRUE(input$gwsGeneCheckContrastAll)){
      presel_contrasts <- gwsGeneContrastList()
    }else{
      presel_contrasts <- local(input$gwsGeneContrastSelect)
    }
    
    #get selected tissue
    if(isTRUE(input$gwsGeneCheckTissueAll)){
      preselTissue <- pheno %>%
        filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
        select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }else{
      if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
        preselTissue <- pheno %>%
          filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
          select(tissue_name) %>%
          distinct() %>%
          .$tissue_name
      }
    }
   
    #get selected library
    if(isTRUE(input$gwsGeneCheckLibraryAll)){
      preselLibrary <- libraries %>%
        filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
        select(library_id) %>%
        .$library_id
    }else{
      if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
        preselLibrary <- libraries %>%
          filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, library_id %in% input$gwsGeneLibrarySelect) %>%
          select(library_id) %>%
          .$library_id
      }
    }
    
    gene_list_screens %>%
      filter(library_id %in% preselLibrary) %>% 
      dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
      arrange(gene) %>%
      .$gene
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

observeEvent(input$gwsGeneSearchRadio, {
  if(input$gwsGeneSearchRadio == "guide_id"){
    disable("gwsGeneInclude")
  }else{
    enable("gwsGeneInclude")
  }
})

observeEvent(input$gwsGeneSpeciesSelect, {
  updateSelectizeInput(session, 'gwsGeneDatasetSelect', choices = dataset_selection_all, selected = "dropout", server = TRUE)
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
  updateSelectizeInput(session, 'gwsGeneTissueSelect', choices = gwsGeneTissueList(), server = TRUE)
  if(input$gwsGeneDatasetSelect %in% c("drug_modifier", "facs_based")){
    updateRadioButtons(session, 'gwsGeneIndexRadio',selected = "lfc")
    disable("gwsGeneIndexRadio")
  }else{
    enable("gwsGeneIndexRadio")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
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
    if(input$gwsGeneSpeciesSelect == "all"){
      paste0(paste(c(table$Symbol_human %>% unique, table$Symbol_mouse %>% unique),collapse="_"), ".txt")
    }else{
      paste0(paste(table$symbol %>% unique,collapse="_"), ".txt")
    }
  },
  content = function(file) {
    table <- gwsGeneDataFrame()
    table %>% select(-contains("Action_sgRNA"), -contains("Action_gene")) %>%
      write_tsv(file)
  }
)