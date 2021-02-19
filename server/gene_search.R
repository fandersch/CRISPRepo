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
      if(is.null(input$gwsGeneCellLineSelect) & !isTRUE(input$gwsGeneCheckCellLineAll)){
        "INFO: Please select the cell line(s) of considered screens in the right panel!"
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
    }
  })
}

#at initial load display nothing
output$gwsGeneTable <- renderDataTable({
})

#upon load display nothing
output$gwsGeneContrastTable <- renderDataTable({
})

#upon load display nothing
output$gwsGeneSampleTable <- renderDataTable({
})

#query database and create dataframe
gwsGeneContrastDataFrame <- reactive({
  
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    presel_contrasts <- gwsGeneContrastList()
  }else{
    presel_contrasts <- local(input$gwsGeneContrastSelect)
  }
  
  df <- con %>%
    tbl("contrasts") %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct() %>%
    collect()
  
  df <- df %>%
    dplyr::select(contrast_id, contrast_id_QC, treatment, control, type, cellline_name, tissue_name, tissue_context, library_id, species, 
           auc, ssmd, dynamic_range, average_lfc_essentials, 
           norm_method_mageck = norm_method, fdr_method_mageck = fdr_method, lfc_method_mageck=lfc_method, cnv_correction_mageck=cnv_correction, 
           filter_mageck=filter, variance_estimation_mageck=variance_estimation)
  
  if(input$gwsGeneDisplayName == "short"){
    df <- df %>%
      dplyr::select(contrast_id_QC, contrast_id, everything())
  }
  
  if(!input$gwsGeneDatasetSelect %in% c("dropout")){
    df <- df %>%
      dplyr::select(-contrast_id_QC)
  } 
  
  df
})

#query database and create dataframe
gwsGeneSampleDataFrame <- reactive({
  
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    presel_contrasts <- gwsGeneContrastList()
  }else{
    presel_contrasts <- local(input$gwsGeneContrastSelect)
  }
  
  contrasts <- con %>%
    tbl("contrasts") %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct() %>%
    collect()
  
  sample_names <- c(contrasts$treatment %>% as.character %>% str_split(",") %>% unlist, contrasts$control %>% as.character %>% str_split(",") %>% unlist)

  df <- pheno %>%
    dplyr::filter(sample_id %in% sample_names) %>%
    distinct %>%
    dplyr::select(sample_id, cellline_name, tissue_name, tissue_context, library_id, species, type,
           condition, mutation, clone, time, treatment, dose, pcr_strategy, facs_gate, replicate, scientist, publication, legacy_sample_name,
           barcode, vbcf_id, raw_file) %>%
    collect
  
})

#create datatable out of dataframe
gwsGeneContrastDataTable <- eventReactive(input$gwsGeneLoadButton,{
  
  df <- gwsGeneContrastDataFrame()
  if (nrow(df) > 0) {
    nfreezeColumns <- 1
    
    dt <- df %>%
      DT::datatable(extensions = c('FixedColumns','FixedHeader'), 
                    options = list(
                      autoWidth = FALSE,
                      headerCallback = JS(headerCallback),
                      scrollX=TRUE,
                      fixedColumns = list(leftColumns = nfreezeColumns),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
                      searchHighlight = TRUE
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
      formatStyle(names(df),"white-space"="nowrap")
    
    output$gwsGeneInfo <- renderText({
      "Info: Loading completed!"
    })
    
    #display datatable
    dt 
  }else{
    output$gwsGeneInfo <- renderText({
      "WARNING: No data found!"
    })
  }
})

#create datatable out of dataframe
gwsGeneSampleDataTable <- eventReactive(input$gwsGeneLoadButton,{
  
  df <- gwsGeneSampleDataFrame()
  if (nrow(df) > 0) {
    nfreezeColumns <- 1
    
    dt <- df %>%
      DT::datatable(extensions = c('FixedColumns','FixedHeader'), 
                    options = list(
                      autoWidth = FALSE,
                      headerCallback = JS(headerCallback),
                      scrollX=TRUE,
                      fixedColumns = list(leftColumns = nfreezeColumns),
                      columnDefs = list(list(className = 'dt-center', targets = "_all")),
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
                      searchHighlight = TRUE
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
      formatStyle(names(df),"white-space"="nowrap")
    
    output$gwsGeneInfo <- renderText({
      "Info: Loading completed!"
    })
    
    #display datatable
    dt 
  }else{
    output$gwsGeneInfo <- renderText({
      "WARNING: No data found!"
    })
  }
})


gwsGeneDataFrame <- reactive({
  #specify which table to select
  if(input$gwsGeneSearchRadio == "guide_id"){
    tableSelect <- "guide_stats"
    search_radio <- "guide_id, gene_id"
  }else{
    tableSelect <- "gene_stats"
    search_radio <- "gene_id"
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
  
  gene_select <- presel_entrez
  
  if(statistics_columns_negative  == "_negative" | statistics_columns_positive == "_positive"){
    select_str <- paste("contrast_id", 
                        search_radio, 
                        local(input$gwsBrowseScreenIndexRadio), sep= ", ")
    
  }else{
    select_str <- paste("contrast_id", 
                        search_radio, 
                        local(input$gwsBrowseScreenIndexRadio), 
                        paste(statistics_columns_negative, collapse=", "), 
                        paste(statistics_columns_positive, collapse=", "), sep= ", ")
  }
  
  if(length(presel_contrasts)>=900){
    contrasts_filter_str <- c(paste(paste("contrast_id", paste0("'", presel_contrasts[1:899], "'"), sep="="), collapse=" OR "))
    i<-900
    while(i <= length(presel_contrasts)){
      end<-i+899
      if(end>length(presel_contrasts)){
        end<-length(presel_contrasts)
      }
      contrasts_filter_str <- c(contrasts_filter_str, paste(paste("contrast_id", paste0("'", presel_contrasts[i:end], "'"), sep="="), collapse=" OR "))
      i<-i+end
    }
  }else{
    contrasts_filter_str <- paste(paste("contrast_id", paste0("'", presel_contrasts, "'"), sep="="), collapse=" OR ")
  }
  
  gene_filter_str <- paste(paste("gene_id", paste0(presel_entrez), sep="="), collapse=" OR ")
  
  query <- paste0("SELECT ", select_str, " FROM ", tableSelect,
                  " WHERE (", gene_filter_str, ") ",
                  "AND (", contrasts_filter_str, ") ")
  
  df<- NULL
  for(z in 1:length(query)){
    # a chunk at a time
    res <- dbSendQuery(con, query[z])
    i<-1
    while(!dbHasCompleted(res)){
      chunk <- dbFetch(res, n = 5000000)
      if(is.null(df)){
        df <- chunk
      }else{
        df <- df %>% rbind(chunk)
      }
      if(i %% 10==0){
        gc()
      }
      i<-i+1
    }
    chunk<-NULL
    dbClearResult(res)
  }
  
  df <- df %>%
    left_join(contrasts, by="contrast_id")
  
  if(input$gwsGeneSearchRadio == "guide_id"){
    if(input$gwsGeneSpeciesSelect == "all"){
      sgRNAs <- con_sgRNAs %>%
        tbl("sgRNAs_human") %>%
        dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
        dplyr::select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
        mutate(guide_id = paste0(entrez_id, "_", sgRNA_23mer)) %>%
        dplyr::select(-entrez_id) %>%
        rbind(
          con_sgRNAs %>%
            tbl("sgRNAs_mouse") %>%
            dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
            dplyr::select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
            collect() %>%
            mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
            mutate(guide_id = paste0(entrez_id, "_", sgRNA_23mer)) %>%
            dplyr::select(-entrez_id)
        )
    }else{
      #get guide ranks
      sgRNAs <- con_sgRNAs %>%
        tbl(tableSgRNAs) %>%
        dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
        dplyr::select(entrez_id = EntrezID, sgRNA_23mer, `VBC-Score`=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
        mutate(guide_id = paste0(entrez_id, "_", sgRNA_23mer)) %>%
        dplyr::select(-entrez_id)
    }
    
    features_buff <- features %>% 
      dplyr::select(-context, -library_id) %>%
      dplyr::filter(gene_id %in% presel_entrez) %>%
      collect()
    
    #join guide ranks
    df <- df %>%
      left_join(features_buff, by = c("gene_id", "guide_id")) %>%
      dplyr::select(contrast_id, contrast_id_QC, guide_id, gene_id, local(input$gwsBrowseScreenIndexRadio), symbol, entrez_id, sequence, species) %>%
      left_join(sgRNAs, by="guide_id") %>%
      distinct() %>%
      mutate_at(c(local(input$gwsBrowseScreenIndexRadio)), round, 3) %>%
      arrange(symbol)
    
  }else{
    df <- df %>%
      left_join(gene_list_screens %>% dplyr::select(-library_id) %>% distinct, by ="gene_id") %>%
      dplyr::select(contrast_id, contrast_id_QC, gene_id, local(input$gwsBrowseScreenIndexRadio), symbol, entrez_id, species, dplyr::one_of(statistics_columns_negative), dplyr::one_of(statistics_columns_positive)) %>%
      mutate_at(c(local(input$gwsBrowseScreenIndexRadio)), round, 3) %>%
      arrange(symbol)
  }

  if(input$gwsGeneSpeciesSelect == "all"){
    df_human <- df %>% 
      dplyr::filter(species == "human") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_human")) %>%
      mutate(EntrezID_human = entrez_id) %>%
      dplyr::rename(Symbol_human = symbol) %>%
      dplyr::select(-entrez_id)
    
    df_mouse <- df %>% 
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_mouse")) %>%
      mutate(EntrezID_mouse = entrez_id) %>%
      dplyr::rename(Symbol_mouse = symbol) %>%
      dplyr::select(-entrez_id)
    
    df <- df_human %>% rbind(df_mouse)
    
    df_human<-NULL
    df_mouse<-NULL
    gc()
  }
  
  #remove duplicate entries to enable spread
  remove_duplicates <- df[df %>% dplyr::select(contrast_id, local(input$gwsGeneSearchRadio)) %>% duplicated,]
  df <- df %>%
    anti_join(remove_duplicates)
  
  if (nrow(df) > 0) {
    presel_contrasts <- df$contrast_id %>% unique
    
    dt <- df %>%
      dplyr::select(contrast_id, matches("entrez_id"), matches("symbol"), matches("Symbol_human"), matches("EntrezID_human"), 
                    matches("Symbol_mouse"), matches("EntrezID_mouse"), local(input$gwsGeneIndexRadio), matches("sequence"), 
             matches("Length"), matches("guide_id"), matches("VBC-score"), matches("rank_overall"), matches("rank_validation")) %>%
      pivot_wider(names_fro="contrast_id", values_from=local(input$gwsGeneIndexRadio))
      # distinct() %>%
    
    if(input$gwsGeneSearchRadio == "gene_id" & !is.null(local(input$gwsGeneInclude)) & length(local(input$gwsGeneInclude))!=0){
      column <- local(input$gwsGeneInclude)
      if("p"  %in% local(input$gwsGeneInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(p = ifelse(p_positive < p_negative, p_positive, p_negative)) %>%
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), p) %>%
              # distinct %>%
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
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), fdr) %>%
              # distinct %>%
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
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), guides) %>%
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
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), guides_good) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="guides_good") %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_GUIDES_GOOD'))
          )
      }
      dt <- dt %>%
        dplyr::select(sort(tidyselect::peek_vars()))
    }
    
    if(input$gwsGeneSpeciesSelect == "all"){
      dt <- dt %>%
        dplyr::select("Symbol_human", "EntrezID_human","Symbol_mouse","EntrezID_mouse", everything())
    }
      
    #create actionButtons for each rows
    Action_gene <- shinyInput(actionButton, nrow(dt), 'button_gene_', label = HTML("gene <br/> predictions"), onclick = 'Shiny.onInputChange(\"select_button_gene\",  this.id)' )
    dt <- dt %>% cbind(Action_gene)
    if(input$gwsGeneSearchRadio == "guide_id"){
      Action_sgRNA <- shinyInput(actionButton, nrow(dt), 'button_sgRNA_', label = HTML("sgRNA <br/> info"), onclick = 'Shiny.onInputChange(\"select_button_sgRNA\",  this.id)' )
      dt <- dt %>%
        cbind(Action_sgRNA)
    }
    
    df<-NULL
    gc()
    
    #remove sensitive data for external view
    if(view=="external"){
      dt %>%
        dplyr::select(-contains("rank_overall"), -contains("rank_validation"), -contains("Action_sgRNA"), -contains("Action_gene"))
    }else{
      dt
    }
  }
})

#make 
gwsGeneDataTable <- eventReactive(input$gwsGeneLoadButton,{
  gwsGeneDatatable <- gwsGeneDataFrame()
  if (nrow(gwsGeneDatatable) > 0) {
    if(!is.null(input$gwsGeneGeneSelect)){
      output$gwsGeneInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    
    if (!is.null(gwsGeneDatatable) & nrow(gwsGeneDatatable) > 0) {
      #make color interval for heatmap
      if(input$gwsGeneIndexRadio == "effect"){
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
      if(input$gwsGeneDatasetSelect %in% c("dropout") & input$gwsGeneDisplayName == "short"){
        for(i in 1:length(colnames_gwsGeneDatatable)){
          if(i < length(colnames_gwsGeneDatatable)){
            tooltip <- paste0(tooltip, "'", colnames_gwsGeneDatatable[i], "'",  ", " )
          }else{
            tooltip <- paste0(tooltip, "'", colnames_gwsGeneDatatable[i], "'")
          }
          colname_buff <- colnames_gwsGeneDatatable[i]
          if(str_ends(colname_buff, "_P")){
            colname_buff <- substr(colname_buff,1,nchar(colname_buff)-2)
          }
          if(str_ends(colname_buff, "_FDR")){
            colname_buff <- substr(colname_buff,1,nchar(colname_buff)-4)
          }
          if(str_ends(colname_buff, "_GUIDES")){
            colname_buff <- substr(colname_buff,1,nchar(colname_buff)-7)
          }
          if(str_ends(colname_buff, "_GUIDES_GOOD")){
            colname_buff <- substr(colname_buff,1,nchar(colname_buff)-12)
          }
          if(colname_buff %in% presel_contrasts){
            colnames_gwsGeneDatatable[i] <- str_replace(colnames_gwsGeneDatatable[i], colname_buff, (contrasts %>% dplyr::select(contrast_id, contrast_id_QC) %>% dplyr::filter(contrast_id == colname_buff) %>% .$contrast_id_QC))
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
    dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
})

gwsGeneCellLineList <- reactive({
  if(isTRUE(input$gwsGeneCheckTissueAll) | !is.null(input$gwsGeneTissueSelect)){
    
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      distinct %>%
      arrange(cellline_name) %>%
      .$cellline_name
  }
})


gwsGeneLibraryList <- reactive({
  if(isTRUE(input$gwsGeneCheckCellLineAll) | !is.null(input$gwsGeneCellLineSelect)){
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      distinct %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- pheno %>%
        dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% input$gwsGeneCellLineSelect) %>%
        dplyr::select(cellline_name) %>%
        distinct %>%
        arrange(cellline_name) %>%
        .$cellline_name
    }
    
    libraries %>%
      dplyr::filter(species %in% speciesList, 
             type %in% input$gwsGeneDatasetSelect, 
             tissue_name %in% preselTissue, 
             cellline_name %in% preselCellline) %>%
      dplyr::select(library_id) %>%
      distinct %>%
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
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- pheno %>%
        dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% input$gwsGeneCellLineSelect) %>%
        dplyr::select(cellline_name) %>%
        arrange(cellline_name) %>%
        .$cellline_name
    }
    
    #get selected library
    preselLibrary <- libraries %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% preselCellline) %>%
      dplyr::select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
      preselLibrary <- libraries %>%
        dplyr::filter(species %in% speciesList, 
               type %in% input$gwsGeneDatasetSelect, 
               tissue_name %in% preselTissue, 
               cellline_name %in% preselCellline,
               library_id %in% input$gwsGeneLibrarySelect) %>%
        dplyr::select(library_id) %>%
        .$library_id
    }
    
    presel_contrasts <- contrasts %>%
      dplyr::filter(species %in% speciesList,
             library_id %in% preselLibrary,
             tissue_name %in%  preselTissue,
             type == input$gwsGeneDatasetSelect,
             cellline_name %in% preselCellline
      ) 
    
    if(input$gwsGeneDatasetSelect == "dropout"){
      presel_contrasts <- presel_contrasts %>%
        dplyr::filter(abs(dynamic_range) >= input$gwsGeneQuality)
    }
    
    presel_contrasts %>%
      dplyr::select(contrast_id) %>%
      distinct() %>%
      .$contrast_id
  }
})

gwsGeneGeneList <- reactive({
  
  if((!is.null(input$gwsGeneContrastSelect) | isTRUE(input$gwsGeneCheckContrastAll)) & 
     (!is.null(input$gwsGeneLibrarySelect) | isTRUE(input$gwsGeneCheckLibraryAll)) & 
     (!is.null(input$gwsGeneTissueSelect) | isTRUE(input$gwsGeneCheckTissueAll))){
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
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- pheno %>%
        dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% input$gwsGeneCellLineSelect) %>%
        dplyr::select(cellline_name) %>%
        arrange(cellline_name) %>%
        .$cellline_name
    }
    
    #get selected library
    preselLibrary <- libraries %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% preselCellline) %>%
      dplyr::select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
      preselLibrary <- libraries %>%
        dplyr::filter(species %in% speciesList, 
               type %in% input$gwsGeneDatasetSelect, 
               tissue_name %in% preselTissue, 
               cellline_name %in% preselCellline,
               library_id %in% input$gwsGeneLibrarySelect) %>%
        dplyr::select(library_id) %>%
        .$library_id
    }
    
    gene_list_screens %>%
      dplyr::filter(library_id %in% preselLibrary) %>% 
      dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
      distinct %>%
      arrange(gene) %>%
      .$gene
  }
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------
observeEvent(input$gwsGeneLoadButton, {
  output$gwsGeneTable <- renderDataTable({
    dt <- gwsGeneDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  output$gwsGeneContrastTable <- renderDataTable({
    dt <- gwsGeneContrastDataTable()
    if(!is.null(dt)){
      dt
    }
  })
  output$gwsGeneSampleTable <- renderDataTable({
    dt <- gwsGeneSampleDataTable()
    if(!is.null(dt)){
      dt
    }
  })
})

#actionButton handler for datatable buttons
observeEvent(input$select_button_gene, {
  row <- as.numeric(strsplit(input$select_button_gene, "_")[[1]][3])
  if(input$gwsGeneSpeciesSelect == "all"){
    symbol_human <- gwsGeneDataFrame()[row,1]
    entrez_id_human <- gwsGeneDataFrame()[row,2]
    symbol_mouse <- gwsGeneDataFrame()[row,3]
    entrez_id_mouse <- gwsGeneDataFrame()[row,4]
    gene<-c()
    if(!is.na(entrez_id_human) & !is.na(symbol_human)){
      gene = c(ifelse(is.na(symbol_human), paste0("No symbol found (", entrez_id_human, ")"), paste0(symbol_human , " (", entrez_id_human, ")")))
    }
    if(!is.na(entrez_id_mouse) & !is.na(symbol_mouse)){
      gene = c(gene,
               ifelse(is.na(symbol_mouse), paste0("No symbol found (", entrez_id_mouse, ")"), paste0(symbol_mouse , " (", entrez_id_mouse, ")")))
    }
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
    symbol_human <- gwsGeneDataFrame()[row,1]
    entrez_id_human <- gwsGeneDataFrame()[row,2]
    symbol_mouse <- gwsGeneDataFrame()[row,3]
    entrez_id_mouse <- gwsGeneDataFrame()[row,4]
    gene<-c()
    if(!is.na(entrez_id_human) & !is.na(symbol_human)){
      gene = c(ifelse(is.na(symbol_human), paste0("No symbol found (", entrez_id_human, ")"), paste0(symbol_human , " (", entrez_id_human, ")")))
    }
    if(!is.na(entrez_id_mouse) & !is.na(symbol_mouse)){
      gene = c(gene,
               ifelse(is.na(symbol_mouse), paste0("No symbol found (", entrez_id_mouse, ")"), paste0(symbol_mouse , " (", entrez_id_mouse, ")")))
    }
    guide_id<-gwsGeneDataFrame()[row,7]
  }else{
    symbol <- gwsGeneDataFrame()[row,2]
    entrez_id <- gwsGeneDataFrame()[row,1]
    gene <- ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))
    guide_id <- gwsGeneDataFrame()[row,5]
  }
  
  updateSelectizeInput(session, 'sgRNAInfoSelectGene', choices = sgRNAInfoGeneList(), selected = gene, server = TRUE)
  delay(1000, updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = guide_id, server = TRUE))
  updateTabItems(session, "tabs", "sgRNAInfoSidebar")
  enable("sgRNAInfoLoadButton")
  delay(1500, click("sgRNAInfoLoadButton"))
})

observeEvent(input$gwsGeneSearchRadio, {
  if(input$gwsGeneSearchRadio == "guide_id"){
    disable("gwsGeneInclude")
    updateCheckboxInput(session, 'gwsGeneInclude', value="")
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
  updateSelectizeInput(session, 'sgRNAInfoSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = input$gwsGeneSpeciesSelect, server = TRUE)
  #disable load button
  disable("gwsGeneLoadButton")
  gwsGeneUpdateText()
})

observeEvent(input$gwsGeneQuality, {
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  #disable laod button
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
    disable("gwsGeneQuality")
  }else{
    enable("gwsGeneIndexRadio")
    enable("gwsGeneQuality")
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
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsGeneCheckTissueAll', value = FALSE)
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneCellLineSelect', choices = gwsGeneCellLineList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckCellLineAll', value = FALSE)
  
  #update contrasts selectb
  if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
    enable("gwsGeneCellLineSelect")
    enable("gwsGeneCheckCellLineAll")
  }else{
    disable("gwsGeneCellLineSelect")
    disable("gwsGeneCheckCellLineAll")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCheckTissueAll, {
  if(isTRUE(input$gwsGeneCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsGeneTissueSelect', choices = gwsGeneTissueList(), server = TRUE)
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneCellLineSelect', choices = gwsGeneCellLineList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckCellLineAll', value = FALSE)
  
  if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
    enable("gwsGeneCellLineSelect")
    enable("gwsGeneCheckCellLineAll")
  }else{
    disable("gwsGeneCellLineSelect")
    disable("gwsGeneCheckCellLineAll")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCellLineSelect, {
  if(!is.null(input$gwsGeneCellLineSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsGeneCheckCellLineAll', value = FALSE)
  }
  #update cellline selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  
  #select cellline checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  if(isTRUE(input$gwsGeneCheckCellLineAll) | (!is.null(input$gwsGeneCellLineSelect))){
    enable("gwsGeneLibrarySelect")
    enable("gwsGeneCheckLibraryAll")
  }else{
    disable("gwsGeneLibrarySelect")
    disable("gwsGeneCheckLibraryAll")
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCheckCellLineAll, {
  if(isTRUE(input$gwsGeneCheckCellLineAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsGeneCellLineSelect', choices = gwsGeneCellLineList(), server = TRUE)
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsGeneLibrarySelect', choices = gwsGeneLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsGeneCheckLibraryAll', value = FALSE)
  if(isTRUE(input$gwsGeneCheckTissueAll) | (!is.null(input$gwsGeneTissueSelect))){
    enable("gwsGeneLibrarySelect")
    enable("gwsGeneCheckLibraryAll")
  }else{
    disable("gwsGeneLibrarySelect")
    disable("gwsGeneCheckLibraryAll")
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsGeneContrastSelect', choices = gwsGeneContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

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
    
    presel_contrasts <- gwsGeneContrastList()
    colnames_gwsGeneDatatable <- colnames(table)
    if(input$gwsGeneDatasetSelect %in% c("dropout") & input$gwsGeneDisplayName == "short"){
      for(i in 1:length(colnames_gwsGeneDatatable)){
        colname_buff <- colnames_gwsGeneDatatable[i]
        if(str_ends(colname_buff, "_P")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-2)
        }
        if(str_ends(colname_buff, "_FDR")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-4)
        }
        if(str_ends(colname_buff, "_GUIDES")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-7)
        }
        if(str_ends(colname_buff, "_GUIDES_GOOD")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-12)
        }
        if(colname_buff %in% presel_contrasts){
          colnames_gwsGeneDatatable[i] <- str_replace(colnames_gwsGeneDatatable[i], colname_buff, (contrasts %>% dplyr::select(contrast_id, contrast_id_QC) %>% dplyr::filter(contrast_id == colname_buff) %>% .$contrast_id_QC))
        }
      }
      
      colnames(table) <- make.unique(colnames_gwsGeneDatatable)
    }
    
    
    table %>% dplyr::select(-contains("Action_sgRNA"), -contains("Action_gene")) %>%
      write_tsv(file)
  }
)