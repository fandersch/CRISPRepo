# ----------------------------------------------------------------------------
# Gene Search
# ----------------------------------------------------------------------------
gwsGeneGeneInputFile <- reactiveValues(data = NULL)

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
            if(is.null(input$gwsGeneGeneSelect)& is.null(gwsGeneGeneInputFile$data)){
              invisible(paste0(" INFO: Please select the gene(s) you want to browse in the right panel! ", HTML('<br/>'),
                               " You can also upload a list of genes with the provided file browser on the right side! ", HTML('<br/>'),
                               " The file must have one gene per row (one single column): either entrez-IDs or gene-symbol (no header, no mixed IDs)!"))
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
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
  df <- con %>%
    tbl("contrasts") %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct() %>%
    collect()
  
  DBI::dbDisconnect(con)
  
  df <- df %>%
    dplyr::select(contrast_id, contrast_id_QC, treatment, control, type, reference_type, cellline_name, tissue_name, tissue_context, library_id, species, 
                  auc, ssmd, dynamic_range, average_lfc_essentials, 
                  norm_method_mageck = norm_method, fdr_method_mageck = fdr_method, lfc_method_mageck=lfc_method, cnv_correction_mageck=cnv_correction, 
                  filter_mageck=filter, variance_estimation_mageck=variance_estimation)
  
  if(input$gwsGeneDisplayName == "short"){
    df <- df %>%
      dplyr::select(contrast_id_QC, contrast_id, everything())
  }
  
  if(!("dropout" %in% input$gwsGeneDatasetSelect)){
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
  
  contrasts <- contrasts %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct()
  
  sample_names <- c(contrasts$treatment %>% as.character %>% str_split(",") %>% unlist, contrasts$control %>% as.character %>% str_split(",") %>% unlist)
  
  df <- pheno %>%
    dplyr::filter(sample_id %in% sample_names) %>%
    distinct %>%
    dplyr::select(sample_id, cellline_name, tissue_name, tissue_context, library_id, species, type,
                  condition, mutation, clone, time, treatment, dose, pcr_strategy, facs_gate, replicate, scientist, publication, legacy_sample_name,
                  barcode, vbcf_id, raw_file)
  
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
                      pageLength = 10,
                      lengthMenu = c(10, 25, 50, 100, 200),
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
                      pageLength = 10,
                      lengthMenu = c(10, 25, 50, 100, 200),
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
    search_radio <- "guide_id, id_entrez_23mer, gene_id"
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
  
  if(!is.null(gwsGeneGeneInputFile$data)){
    presel_genes_buff <- gwsGeneGeneList()
    genes_fileUpload <- c(paste0("\\(", (gwsGeneGeneInputFile$data$X1 %>% as.character), "\\)"), paste0("^", (gwsGeneGeneInputFile$data$X1 %>% as.character), "\\s"))
    presel_genes_both <- grep(paste(genes_fileUpload,collapse="|"), presel_genes_buff %>% as.character, value=TRUE)
  }else{
    presel_genes_both<- input$gwsGeneGeneSelect
  }
  
  #retrieve selected genes
  presel_genes_both<- presel_genes_both %>% strsplit(split="\\(|\\)")
  presel_genes <- unlist(presel_genes_both)[c(TRUE, FALSE)] %>% trimws()
  presel_entrez <- unlist(presel_genes_both)[c(FALSE, TRUE)]
  gene_select <- c(presel_entrez, presel_genes)
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")
  
  #retrieve selected contrasts
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    presel_contrasts <- gwsGeneContrastList()
  }else{
    presel_contrasts <- local(input$gwsGeneContrastSelect)
  }
  
  include_columns <- local(input$gwsGeneInclude)
  selected_index <- local(input$gwsGeneIndexRadio)
  
  select_str <- paste(c("contrast_id", 
                        search_radio, 
                        selected_index, 
                        include_columns), collapse= ", ")
  
  
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
  
  gene_filter_str <- paste(paste("gene_id", paste0("'",gene_select,"'"), sep="="), collapse=" OR ")
  
  query <- paste0("SELECT ", select_str, " FROM ", tableSelect,
                  " WHERE (", gene_filter_str, ") ",
                  "AND (", contrasts_filter_str, ") ")
  
  df <- NULL
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
        dplyr::select(entrez_id = EntrezID, sgRNA_23mer, "VBC-Score"=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(entrez_id = as.character(entrez_id), rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
        mutate(guide_id = paste0(entrez_id, "_", sgRNA_23mer)) %>%
        dplyr::select(-entrez_id) %>%
        rbind(
          con_sgRNAs %>%
            tbl("sgRNAs_mouse") %>%
            dplyr::filter(EntrezID %in% c(presel_entrez)) %>%
            dplyr::select(entrez_id = EntrezID, sgRNA_23mer, "VBC-Score"=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
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
        dplyr::select(entrez_id = EntrezID, sgRNA_23mer, "VBC-Score"=VBC.score, Length = len_cloning_sgRNA, rank_overall = final_rank, rank_validation = final_validated_rank) %>%
        collect() %>%
        mutate(rank_overall = as.character(rank_overall), rank_validation = as.character(rank_validation)) %>%
        mutate(guide_id = paste0(entrez_id, "_", sgRNA_23mer)) %>%
        dplyr::select(-entrez_id)
    }
    
    features_buff <- features %>% 
      dplyr::filter(gene_id %in% presel_entrez | gene_id %in% presel_genes)
    
    #join guide ranks
    df <- df %>%
      left_join(features_buff, by = c("gene_id", "guide_id", "id_entrez_23mer", "library_id")) %>%
      dplyr::select(contrast_id, contrast_id_QC, guide_id, id_entrez_23mer, gene_id, selected_index, symbol, entrez_id, sequence, sequence_matching, species, dplyr::any_of(include_columns)) %>%
      left_join(sgRNAs, by=c("id_entrez_23mer" = "guide_id")) %>%
      mutate(guide_id = id_entrez_23mer) %>%
      dplyr::select(-id_entrez_23mer) %>%
      distinct() %>%
      arrange(symbol)
    
  }else{
    df <- df %>%
      left_join(gene_list_screens, by=c("gene_id", "library_id")) %>%
      dplyr::select(contrast_id, contrast_id_QC, gene_id, selected_index, symbol, entrez_id, species, dplyr::any_of(include_columns)) %>%
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
    
    #cleanup
    df_human<-NULL
    df_mouse<-NULL
    gc()
  }
  
  column_names_change <- c(effect = "effect_essentialome", adjusted_effect = "adjusted_effect_essentialome")
  selected_index <- str_remove(input$gwsGeneIndexRadio, pattern = "_essentialome")
  
  df <- df %>%
    dplyr::rename(any_of(column_names_change)) %>%
    discard(~all(is.na(.) | . =="")) %>%
    dplyr::select(contrast_id, dplyr::starts_with(c("symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse", "guide_id", 
                                                    "sequence", "sequence_matching", "Length", "VBC", "rank")), selected_index, dplyr::any_of(include_columns)) %>%
    pivot_wider(names_from=contrast_id, values_from=paste0(c(selected_index, include_columns)), names_glue = "{contrast_id}_{.value}") %>%
    arrange(dplyr::across(starts_with("symbol")))
  
  #order columns
  colnames(df) <- colnames(df) %>%
    str_replace(paste0(selected_index, "$"), paste0("1", selected_index)) %>%
    str_replace("p$", "2p") %>%
    str_replace("fdr$", "3fdr") %>%
    str_replace("guides_good$", "4guides_good") %>%
    str_replace("guides$", "5guides")
  
  df <- df %>%
    dplyr::select(order(colnames(df))) %>%
    dplyr::select(dplyr::starts_with(c("symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse", "guide_id", 
                                       "sequence", "sequence_matching", "Length", "VBC", "rank")), everything())
  
  colnames(df) <- colnames(df) %>%
    str_replace(paste0("1", selected_index, "$"), selected_index) %>%
    str_replace("2p$", "p") %>%
    str_replace("3fdr$", "fdr") %>%
    str_replace("4guides_good$", "guides_good") %>%
    str_replace("5guides$", "guides")
  
  DBI::dbDisconnect(con)
  DBI::dbDisconnect(con_sgRNAs)
  
  df
})

#make 
gwsGeneDataTable <- eventReactive(input$gwsGeneLoadButton,{
  dt <- gwsGeneDataFrame()
  
  if(isTRUE(input$gwsGeneCheckContrastAll)){
    presel_contrasts <- gwsGeneContrastList()
  }else{
    presel_contrasts <- local(input$gwsGeneContrastSelect)
  }
  
  selected_index <- str_remove(input$gwsGeneIndexRadio, pattern = "_essentialome")
  include_columns <- local(input$gwsGeneInclude)
  
  if(!is.null(dt) & nrow(dt) > 0) {
    
    #make color interval for heatmap
    if(input$gwsGeneIndexRadio %in% c("adjusted_effect_essentialome", "effect_essentialome")){
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

    displayed_table <<- dt

    #create actionButtons for each rows
    Action_gene <- shinyInput(actionButton, nrow(dt), 'button_gene_', label = HTML("gene <br/> predictions"), onclick = 'Shiny.onInputChange(\"select_button_gene\",  this.id)' )
    dt <- dt %>% cbind(Action_gene)
    if(input$gwsGeneSearchRadio == "guide_id"){
      Action_sgRNA <- shinyInput(actionButton, nrow(dt), 'button_sgRNA_', label = HTML("sgRNA <br/> info"), onclick = 'Shiny.onInputChange(\"select_button_sgRNA\",  this.id)' )
      dt <- dt %>%
        cbind(Action_sgRNA)
      
      #specify which columns should be frozen and which should have the heatmap
      nfreezeColumns <- 3
      
      nColorizeTableColumns <- (dt %>%
        dplyr::select(dplyr::starts_with(c("symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse", "guide_id", 
                                           "sequence", "sequence_matching", "Length", "VBC", "rank"))) %>% 
                                           ncol) + 1
                                           
      nActionButtons <- 2
    }else{
      nfreezeColumns <- 2
      nColorizeTableColumns <- 3
      nActionButtons <- 1
    }
    
    df<-NULL
    gc()
    
    #remove sensitive data for external view
    if(view=="external"){
      dt <- dt %>%
        dplyr::select(-contains("rank_overall"), -contains("rank_validation"), -contains("Action_sgRNA"), -contains("Action_gene"))
    }
    
    
    #set colorized column, depending on shown gene-level statistics
    colorInterval <- length(local(input$gwsGeneInclude))
    #set frozen/colourized columns
    if(input$gwsGeneSpeciesSelect == "all"){
      nfreezeColumns <- nfreezeColumns + 2
      nColorizeTableColumns <- nColorizeTableColumns
    }
    
    #change to contrast_id_QC for selected columns when type is dropout
    sample_info <- contrasts %>%
      filter(contrast_id %in% presel_contrasts) %>%
      dplyr::select(contrast_id, treatment, control) %>%
      separate_rows(c("treatment"), sep = ",") %>%
      left_join(pheno, by=c("treatment"="sample_id"))
    
    
    colnames_dt <- colnames(dt)
    tooltip <- ''
    if("dropout" %in% input$gwsGeneDatasetSelect & length(input$gwsGeneDatasetSelect) == 1){
      for(i in 1:length(colnames_dt)){
        colname_buff <- colnames_dt[i]
        
        if(str_ends(colname_buff, "_adjusted_effect")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-16)
        }
        if(str_ends(colname_buff, "_effect")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-7)
        }
        if(str_ends(colname_buff, "_lfc")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-4)
        }
        if(str_ends(colname_buff, "_p")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-2)
        }
        if(str_ends(colname_buff, "_fdr")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-4)
        }
        if(str_ends(colname_buff, "_guides")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-7)
        }
        if(str_ends(colname_buff, "_guides_good")){
          colname_buff <- substr(colname_buff,1,nchar(colname_buff)-12)
        }
        
        #change column name to contrast_id_QC
        if(colname_buff %in% presel_contrasts){
          colnames_dt[i] <- str_replace(colnames_dt[i], colname_buff, (contrasts %>% dplyr::select(contrast_id, contrast_id_QC) %>% dplyr::filter(contrast_id == colname_buff) %>% .$contrast_id_QC))
        }
        
        scientist <- sample_info %>%
          filter(contrast_id == colname_buff)  %>%
          .$scientist %>%
          unique %>%
          paste(collapse=", ")
        
        treatment_sample <- sample_info %>%
          filter(contrast_id == colname_buff)  %>%
          .$treatment %>%
          unique %>%
          paste(collapse=", ")
        
        control_sample <- sample_info %>%
          filter(contrast_id == colname_buff)  %>%
          .$control %>%
          unique %>%
          paste(collapse=", ")
        
        time <- sample_info %>%
          filter(contrast_id == colname_buff)  %>%
          .$time %>%
          unique %>%
          paste(collapse=", ")
        
        tooltip_check <- c(colname_buff,
                           sample_info %>%
                             filter(contrast_id == colname_buff)  %>%
                             .$scientist %>%
                             unique,
                           sample_info %>%
                             filter(contrast_id == colname_buff)  %>%
                             .$treatment %>%
                             unique,
                           sample_info %>%
                             filter(contrast_id == colname_buff)  %>%
                             .$control %>%
                             unique,
                           sample_info %>%
                             filter(contrast_id == colname_buff)  %>%
                             .$time %>%
                             unique)
        
        tooltip_check <- tooltip_check[!is.na(tooltip_check)]
        
        max_string_tooltip <- tooltip_check %>%
          nchar %>%
          max()
        
        if(colname_buff %in% presel_contrasts){
          tooltip_info <-  paste0("Full name: ", str_pad(colname_buff, max_string_tooltip, side ="right"), " ", 
                                  "Scientist:", scientist, " ", 
                                  "Treatment:", treatment_sample, " ", 
                                  "Control:",  str_pad(control_sample,max_string_tooltip, side ="right"), " ", 
                                  "Timepoint:", time)
          
        }else{
          tooltip_info <- colnames_dt[i]
        }
          
        
        if(i < length(colnames_dt)){
          tooltip <- paste0(tooltip, "'", tooltip_info, "'",  ", " )
        }else{
          tooltip <- paste0(tooltip, "'", tooltip_info, "'")
        }
        
      }
      if(input$gwsGeneDisplayName == "short"){
        colnames(dt) <-  make.unique(colnames_dt)
      }
    }
    
    headerCallback <- c(
      "function(thead, data, start, end, display){",
      "  var $ths = $(thead).find('th');",
      "  $ths.css({'vertical-align': 'bottom', 'white-space': 'pre-wrap'});",
      "  var betterCells = [];",
      "  $ths.each(function(){",
      "    var cell = $(this);",
      "    var newDiv = $('<div>', {height: 'auto', width: 'auto'});",
      "    var newInnerDiv = $('<div>', {text: cell.text()});",
      "    newDiv.css({margin: 'auto'});",
      "    newInnerDiv.css({",
      "      transform: 'rotate(180deg)',",
      "      'writing-mode': 'tb-rl',",
      "      'white-space': 'pre-wrap'",
      "    });",
      "    newDiv.append(newInnerDiv);",
      "    betterCells.push(newDiv);",
      "  });",
      "  $ths.each(function(i){",
      "    $(this).html(betterCells[i]);",
      "  });",
      paste0("  var tooltips = [", tooltip, "];"),
      paste0("  for(var i=0; i<", length(colnames_dt), "; i++){"),
      "    $('th:eq('+i+')',thead).attr('data-title', tooltips[i]);",
      "  }",
      "}"
    )
    dt <- dt %>% 
      DT::datatable(escape = FALSE,
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
      formatStyle(seq(nColorizeTableColumns, length(colnames_dt),1+colorInterval), 
                  backgroundColor = styleInterval(brks, clrs))
    
    if(!is.null(input$gwsGeneGeneSelect)){
      output$gwsGeneInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    gc()
    #display datatable
    dt
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
gwsGenenDatasetList <- reactive({
  if(input$gwsGenenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGenenSpeciesSelect
  }
  contrasts %>%
    dplyr::filter(species %in% speciesList, !is.na(type)) %>%
    .$type %>%
    as.character %>%
    unique
})

gwsGenenDataset <- reactive({
  if(isTRUE(input$gwsGenenCheckDatasetAll)){
    gwsGenenDatasetList()
  }else{
    input$gwsGenenDatasetSelect
  }
})

gwsGenenReferenceList <- reactive({
  if(input$gwsGenenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGenenSpeciesSelect
  }
  contrasts %>%
    dplyr::filter(species %in% speciesList, !is.na(reference_type)) %>%
    .$reference_type %>%
    as.character %>%
    unique
})

gwsGenenReference <- reactive({
  if(isTRUE(input$gwsGenenCheckReferenceAll)){
    gwsGenenReferenceList()
  }else{
    input$gwsGenenReferenceSelect
  }
})

gwsGeneTissueList <- reactive({

  if(input$gwsGeneSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGeneSpeciesSelect
  }
  
  contrasts %>%
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
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    contrasts %>%
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
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      distinct %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- contrasts %>%
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
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- contrasts %>%
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
                    type %in% input$gwsGeneDatasetSelect,
                    cellline_name %in% preselCellline
      ) 
    
    if("dropout" %in% input$gwsGeneDatasetSelect & length(input$gwsGeneDatasetSelect) == 1){
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
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsGeneCheckTissueAll) & !is.null(input$gwsGeneTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% input$gwsGeneTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsGeneDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsGeneCheckCellLineAll) & !is.null(input$gwsGeneCellLineSelect)){
      preselCellline  <- contrasts %>%
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
      dplyr::select(gene) %>%
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
    entrez_id_human <- gwsGeneDataFrame()[row,3]
    symbol_mouse <- gwsGeneDataFrame()[row,2]
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
    symbol <- gwsGeneDataFrame()[row,1]
    entrez_id <- gwsGeneDataFrame()[row,2]
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
    entrez_id_human <- gwsGeneDataFrame()[row,3]
    symbol_mouse <- gwsGeneDataFrame()[row,2]
    entrez_id_mouse <- gwsGeneDataFrame()[row,4]
    gene<-c()
    if(!is.na(entrez_id_human) & !is.na(symbol_human)){
      gene = c(ifelse(is.na(symbol_human), paste0("No symbol found (", entrez_id_human, ")"), paste0(symbol_human , " (", entrez_id_human, ")")))
    }
    if(!is.na(entrez_id_mouse) & !is.na(symbol_mouse)){
      gene = c(gene,
               ifelse(is.na(symbol_mouse), paste0("No symbol found (", entrez_id_mouse, ")"), paste0(symbol_mouse , " (", entrez_id_mouse, ")")))
    }
    guide_id<-gwsGeneDataFrame()[row,5]
  }else{
    symbol <- gwsGeneDataFrame()[row,1]
    entrez_id <- gwsGeneDataFrame()[row,2]
    gene <- ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))
    guide_id <- paste0(gwsGeneDataFrame()[row,3])
  }
  
  updateSelectizeInput(session, 'sgRNAInfoSelectGene', choices = sgRNAInfoGeneList(), selected = gene, server = TRUE)
  delay(1000, updateSelectizeInput(session, 'sgRNAInfoSelectGuide', choices = sgRNAInfoGuideList(), selected = guide_id, server = TRUE))
  updateTabItems(session, "tabs", "sgRNAInfoSidebar")
  enable("sgRNAInfoLoadButton")
  delay(1500, click("sgRNAInfoLoadButton"))
})

observeEvent(input$gwsGeneSearchRadio, {
  if(input$gwsGeneSearchRadio == "guide_id"){
    updateCheckboxInput(session, 'gwsGeneInclude', value="")
    disable("gwsGeneInclude")
    
    disable("gwsGeneIndexRadio")
    updateRadioButtons(session, 'gwsGeneIndexRadio', selected = "lfc")
    
    if("dropout" %in% input$gwsGeneDatasetSelect & length(input$gwsGeneDatasetSelect) == 1){
      enable("gwsGeneIndexRadio")
      updateRadioButtons(session,
                         "gwsGeneIndexRadio",
                         label = "Display data as:",
                         choices = list("Log-fold change" = "lfc", "Effect" = "effect_essentialome"),
                         selected = "lfc",
                         inline = F
      )
    }
  }else{
    enable("gwsGeneInclude")
    
    disable("gwsGeneIndexRadio")
    updateRadioButtons(session, 'gwsGeneIndexRadio', selected = "lfc")
    
    if("dropout" %in% input$gwsGeneDatasetSelect & length(input$gwsGeneDatasetSelect) == 1){
      updateRadioButtons(session, "gwsGeneIndexRadio", 
                         choices = list("Log-fold change" = "lfc", "Effect" = "effect_essentialome", "FDR-adjusted effect" = "adjusted_effect_essentialome"), 
                         selected = "adjusted_effect_essentialome"
      )
    }
  }
})

observeEvent(input$gwsGeneSpeciesSelect, {
  if(input$gwsGeneSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsGeneSpeciesSelect
  }
  
  all_types <- contrasts %>%
    dplyr::filter(species %in% speciesList, !is.na(type)) %>%
    .$type %>%
    unique
  
  dataset_selection_all<- setNames(all_types, all_types)
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
  if("dropout" %in% input$gwsGeneDatasetSelect & length(input$gwsGeneDatasetSelect) == 1){
    enable("gwsGeneIndexRadio")
    enable("gwsGeneQuality")
    updateRadioButtons(session, 'gwsGeneDisplayName', selected = "short")
    enable("gwsGeneDisplayName")
  }else{
    updateRadioButtons(session, 'gwsGeneIndexRadio', selected = "lfc")
    disable("gwsGeneIndexRadio")
    disable("gwsGeneQuality")
    updateRadioButtons(session, 'gwsGeneDisplayName', selected = "long")
    disable("gwsGeneDisplayName")
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
    enable("gwsGeneGeneInputFile")
  }else{
    disable("gwsGeneGeneSelect")
    disable("gwsGeneGeneInputFile")
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
    enable("gwsGeneGeneInputFile")
  }else{
    disable("gwsGeneGeneSelect")
    disable("gwsGeneGeneInputFile")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneCancelModal, {
  updateCheckboxInput(session, 'gwsGeneCheckContrastAll', value = FALSE)
  removeModal()
})

observeEvent(input$gwsGeneGeneSelect, {
  if((!is.null(input$gwsGeneGeneSelect)) | !is.null(gwsGeneGeneInputFile$data)){
    enable("gwsGeneLoadButton")
    if(!is.null(input$gwsGeneGeneSelect)){
      reset('gwsGeneGeneInputFile')
      gwsGeneGeneInputFile$data <- NULL
    }
  }else{
    disable("gwsGeneLoadButton")
  }
  gwsGeneUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsGeneGeneInputFile, {
  if(!is.null(input$gwsGeneGeneInputFile)){
    updateSelectizeInput(session, 'gwsGeneGeneSelect', choices = gwsGeneGeneList(), server = TRUE)
    req(input$gwsGeneGeneInputFile)
    gwsGeneGeneInputFile$data <- read_tsv(input$gwsGeneGeneInputFile$datapath, col_names = F)
  }else{
    gwsGeneGeneInputFile$data <- NULL
  }
  if((!is.null(input$gwsGeneGeneSelect)) | !is.null(gwsGeneGeneInputFile$data)){
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
    "crisprepo_gene_search.txt"
  },
  content = function(file) {
    
    displayed_table %>% dplyr::select(-contains("Action_sgRNA"), -contains("Action_gene")) %>%
      write_tsv(file)
  }
)