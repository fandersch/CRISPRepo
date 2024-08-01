# ----------------------------------------------------------------------------
# Browse Screen
# ----------------------------------------------------------------------------

gwsBrowseScreenUpdateText <- function(){
  output$gwsBrowseScreenInfo <- renderText({
    if(is.null(input$gwsBrowseScreenTissueSelect) & !isTRUE(input$gwsBrowseScreenCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) of considered screens in the right panel to start new query!")
    }else{
      if(is.null(input$gwsBrowseScreenCellLineSelect) & !isTRUE(input$gwsBrowseScreenCheckCellLineAll)){
        invisible("INFO: Please select the cell lines(s) of considered screens in the right panel!")
      }else{
        if(is.null(input$gwsBrowseScreenLibrarySelect) & !isTRUE(input$gwsBrowseScreenCheckLibraryAll)){
          invisible("INFO: Please select the librarie(s) of considered screens in the right panel!")
        }else{
          if(is.null(input$gwsBrowseScreenContrastSelect) & !isTRUE(input$gwsBrowseScreenCheckContrastAll)){
            invisible("INFO: Please select the contrasts(s) you want to browse in the right panel!")
          }else{
            invisible("INFO: Click Load data!")
          }
        }
      }
    }
  })
}

#upon load display nothing
output$gwsBrowseScreenTable <- renderDataTable({
  showModal(modalDialog(
    title = "ATTENTION!", 
    paste0("Do you want to load the default table of all FDR-adjusted scaled gene-level LFCs of HQ dropout screen? This will take around 30 seconds to process!"),
    footer = tagList(
      actionButton("gwsBrowseScreenLoadDefaultTableModal", "Yes"),
      modalButton("No")
    )
  ))
})

#upon load display nothing
output$gwsBrowseScreenContrastTable <- renderDataTable({
})

#upon load display nothing
output$gwsBrowseScreenSampleTable <- renderDataTable({
})

#query database and create dataframe
gwsBrowseScreenContrastDataFrame <- reactive({
  
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
    presel_contrasts <- gwsBrowseScreenContrastList()
  }else{
    presel_contrasts <- local(input$gwsBrowseScreenContrastSelect)
  }
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
  df <- con %>%
    tbl("contrasts") %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct() %>%
    collect()
  
  DBI::dbDisconnect(con)
  
  namekey <- c(contrast_id="contrast_id", contrast_id_QC="contrast_id_QC", treatment="treatment", control="control", type="type", reference_type="reference_type", cellline_name="cellline_name", 
               tissue_name= "tissue_name", tissue_context = "tissue_context", library_id="library_id", species="species", 
               auc= "auc", ssmd = "ssmd", dynamic_range = "dynamic_range", average_lfc_essentials = "average_lfc_essentials", 
               norm_method="norm_method_mageck", fdr_method="fdr_method_mageck", lfc_method="lfc_method_mageck", cnv_correction="cnv_correction_mageck",
               filter="filter_mageck", variance_estimation="variance_estimation_mageck")
  
  names(df) <- namekey[names(df)]
  df <- df %>%
    dplyr::select(contrast_id, contrast_id_QC, treatment, control, type, reference_type, cellline_name, tissue_name, tissue_context, library_id, species, 
           auc, ssmd, dynamic_range, average_lfc_essentials, 
           matches("norm_method_mageck"), matches("fdr_method_mageck"), matches("lfc_method_mageck"), matches("cnv_correction_mageck"), 
           matches("filter_mageck"), matches("variance_estimation_mageck"))
  
  if(input$gwsBrowseScreenDisplayName == "short"){
    df <- df %>%
      dplyr::select(contrast_id_QC, contrast_id, everything())
  }
  
  if(!("dropout" %in% input$gwsBrowseScreenDatasetSelect)){
    df <- df %>%
      dplyr::select(-contrast_id_QC)
  } 
  
  df
})

#query database and create dataframe
gwsBrowseScreenSampleDataFrame <- reactive({
  
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
    presel_contrasts <- gwsBrowseScreenContrastList()
  }else{
    presel_contrasts <- local(input$gwsBrowseScreenContrastSelect)
  }
  
  contrasts_buff <- contrasts %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct()
  
  sample_names <- c(contrasts_buff$treatment %>% as.character %>% str_split(",") %>% unlist, contrasts_buff$control %>% as.character %>% str_split(",") %>% unlist)

  df <- pheno %>%
    dplyr::filter(sample_id %in% sample_names) %>%
    distinct %>%
    dplyr::select(sample_id, cellline_name, tissue_name, tissue_context, library_id, species, type,
           condition, mutation, clone, time, treatment, dose, pcr_strategy, facs_gate, replicate, scientist, publication, legacy_sample_name,
           barcode, vbcf_id, raw_file)
})

#create datatable out of dataframe
gwsBrowseScreenContrastDataTable <- eventReactive(input$gwsBrowseScreenLoadButton,{
  
  df <- gwsBrowseScreenContrastDataFrame()
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
    
    
    output$gwsBrowseScreenInfo <- renderText({
      "Info: Loading completed!"
    })

    #display datatable
    dt 
  }else{
    output$gwsBrowseScreenInfo <- renderText({
      "WARNING: No data found!"
    })
  }
})

#create datatable out of dataframe
gwsBrowseScreenSampleDataTable <- eventReactive(input$gwsBrowseScreenLoadButton,{
  
  df <- gwsBrowseScreenSampleDataFrame()
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
    
    output$gwsBrowseScreenInfo <- renderText({
      "Info: Loading completed!"
    })
    
    #display datatable
    dt 
  }else{
    output$gwsBrowseScreenInfo <- renderText({
      "WARNING: No data found!"
    })
  }
})


#query database and create dataframe
gwsBrowseScreenDataFrame <- reactive({
  
  if(input$gwsBrowseScreenSearchRadio == "guide_id"){
    tableSelect <- "guide_stats"
    search_radio <- "guide_id, id_entrez_23mer, gene_id"
  }else{
    tableSelect <- "gene_stats"
    search_radio <- "gene_id"
  }
  
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
    presel_contrasts <- gwsBrowseScreenContrastList()
  }else{
    presel_contrasts <- local(input$gwsBrowseScreenContrastSelect)
  }
  
  include_columns <- local(input$gwsBrowseScreenInclude)
  selected_index <- local(input$gwsBrowseScreenIndexRadio)
  
  
 
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
  
  query <- paste0("SELECT ", select_str, " FROM ", tableSelect,
                  " WHERE NOT gene_id='AMBIGUOUS' AND NOT gene_id='UNMAPPED' AND NOT gene_id='NOFEATURE' AND NOT gene_id='SAFETARGETING' AND NOT gene_id='NONTARGETING' AND gene_id NOT NULL ",
                  "AND (", contrasts_filter_str, ") ")
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")
  
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
  
  DBI::dbDisconnect(con)
  
  df <- df %>%
    left_join(contrasts, by = "contrast_id") %>%
    left_join(gene_list_screens, by = c("gene_id", "library_id")) %>%
    dplyr::select(contrast_id, contrast_id_QC, contains("guide_id"), contains("id_entrez_23mer"), gene_id, selected_index, symbol, entrez_id, species, dplyr::any_of(include_columns)) %>%
    dplyr::mutate(guide_id = ifelse(!is.null(id_entrez_23mer) & !is.na(id_entrez_23mer), id_entrez_23mer, guide_id)) %>%
    dplyr::select(-matches("id_entrez_23mer")) %>%
    distinct()
  
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    
    df_human <- df %>% 
      dplyr::filter(species == "human") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_human")) %>%
      mutate(EntrezID_human = entrez_id) %>%
      dplyr::select(contrast_id, contrast_id_QC, contains("guide_id"), gene_id, selected_index, Symbol_human = symbol, EntrezID_human, Symbol_mouse,  EntrezID_mouse, dplyr::any_of(include_columns))
    
    df_mouse <- df %>% 
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_mouse")) %>%
      mutate(EntrezID_mouse = entrez_id) %>%
      dplyr::select(contrast_id, contrast_id_QC, contains("guide_id"),  gene_id, selected_index, Symbol_human, EntrezID_human, Symbol_mouse = symbol, EntrezID_mouse, dplyr::any_of(include_columns))
    
    df <- df_human %>% rbind(df_mouse)
    
    df_human<-NULL
    df_mouse<-NULL
  }
  
  column_names_change <- c(effect = "effect_essentialome", adjusted_effect = "adjusted_effect_essentialome")
  selected_index <- str_remove(input$gwsBrowseScreenIndexRadio, pattern = "_essentialome")
  
  df <- df %>%
    rename(any_of(column_names_change)) %>%
    discard(~all(is.na(.) | . =="")) %>%
    dplyr::select(contrast_id, any_of(c("guide_id", "symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse")), selected_index, dplyr::any_of(include_columns)) %>%
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
    dplyr::select(any_of(c("guide_id", "symbol", "entrez_id", "Symbol_human", "EntrezID_human", "Symbol_mouse", "EntrezID_mouse")), everything())
  
  colnames(df) <- colnames(df) %>%
    str_replace(paste0("1", selected_index, "$"), selected_index) %>%
    str_replace("2p$", "p") %>%
    str_replace("3fdr$", "fdr") %>%
    str_replace("4guides_good$", "guides_good") %>%
    str_replace("5guides$", "guides")
  
  df

})

#create datatable out of dataframe
gwsBrowseScreenDataTable <- eventReactive(input$gwsBrowseScreenLoadButton,{
  dt <- gwsBrowseScreenDataFrame()
  
  presel_contrasts <- gwsBrowseScreenContrastList()
  selected_index <- str_remove(input$gwsBrowseScreenIndexRadio, pattern = "_essentialome")
  include_columns <- local(input$gwsBrowseScreenInclude)
  
  if (nrow(dt) > 0) {
    # color codig for heatmap
    values <- dt %>% 
      dplyr::select(dplyr::ends_with(selected_index)) %>%
      as.data.frame() %>% as.matrix() %>% as.vector()
    
    brks_smaller <- seq(min(values,na.rm = TRUE), 0, .05)
    brks_bigger <- seq(0, max(values,na.rm = TRUE), .05)
    
    clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}
    clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
    {paste0("rgb(", ., ",", ., ",255)")}
    
    brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
    clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
    
    #decide mode for datatable scrolling
    nfreezeColumns <- 2
    if(input$gwsBrowseScreenSearchRadio == "guide_id"){
      nfreezeColumns <- 3
    }
    
    if(local(input$gwsBrowseScreenSpeciesSelect) == "all"){
          nfreezeColumns <- nfreezeColumns + 2
    }
      
    colorInterval <- length(local(input$gwsBrowseScreenInclude))
    
    #change to contrast_id_QC for selected columns when type is dropout
    sample_info <- contrasts %>%
      filter(contrast_id %in% presel_contrasts) %>%
      dplyr::select(contrast_id, treatment, control) %>%
      separate_rows(c("treatment"), sep = ",") %>%
      left_join(pheno, by=c("treatment"="sample_id"))
    

    colnames_dt <- colnames(dt)
    tooltip <- ''
    if("dropout" %in% input$gwsBrowseScreenDatasetSelect & length(input$gwsBrowseScreenDatasetSelect)==1){
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
          tooltip_info <-  paste0("Full name:", str_pad(colname_buff, max_string_tooltip, side ="right") , " ", 
                                  "Scientist:", scientist, " ", 
                                  "Treatment:", treatment_sample, " ", 
                                  "Control:", control_sample, " ", 
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
      if(input$gwsBrowseScreenDisplayName == "short"){
        colnames(dt) <-  make.unique(colnames_dt)
      }
    }
    
    displayed_table <<- dt
    
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
      "}")
    
    dt <- dt %>%
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
      formatStyle(seq(nfreezeColumns+1, length(colnames_dt),1+colorInterval),
                  backgroundColor = styleInterval(brks, clrs))
    
    if(!is.null(input$gwsBrowseScreenContrastSelect)){
      output$gwsBrowseScreenInfo <- renderText({
        "Info: Loading completed!"
      })
    }
    gc()
    #display datatable
    dt 
  }else{
    if(!is.null(input$gwsBrowseScreenContrastSelect)){
      output$gwsBrowseScreenInfo <- renderText({
        "WARNING: No data found!"
      })
    }else{
      gwsBrowseScreenUpdateText()
    }
  }
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

gwsBrowseScreenTissueList <- reactive({

  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsBrowseScreenSpeciesSelect
  }
  
  contrasts %>%
    dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
    dplyr::select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

gwsBrowseScreenCellLineList <- reactive({
  if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | !is.null(input$gwsBrowseScreenTissueSelect)){
    
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsBrowseScreenSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      distinct %>%
      arrange(cellline_name) %>%
      .$cellline_name
  }
})

gwsBrowseScreenLibraryList <- reactive({
  if(isTRUE(input$gwsBrowseScreenCheckCellLineAll) | !is.null(input$gwsBrowseScreenCellLineSelect)){
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsBrowseScreenSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      distinct %>%
      .$cellline_name
      
    if(!isTRUE(input$gwsBrowseScreenCheckCellLineAll) & !is.null(input$gwsBrowseScreenCellLineSelect)){
      preselCellline  <- contrasts %>%
        dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% input$gwsBrowseScreenCellLineSelect) %>%
        dplyr::select(cellline_name) %>%
        distinct %>%
        arrange(cellline_name) %>%
        .$cellline_name
    }
    
    libraries %>%
      dplyr::filter(species %in% speciesList, 
             type %in% input$gwsBrowseScreenDatasetSelect, 
             tissue_name %in% preselTissue, 
             cellline_name %in% preselCellline) %>%
      dplyr::select(library_id) %>%
      distinct %>%
      arrange(library_id) %>%
      .$library_id
  }
  
})

gwsBrowseScreenContrastList <- reactive({
  if(!is.null(input$gwsBrowseScreenLibrarySelect) | isTRUE(input$gwsBrowseScreenCheckLibraryAll)){
    if(local(input$gwsBrowseScreenSpeciesSelect) == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsBrowseScreenSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- contrasts %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- contrasts %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- contrasts %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckCellLineAll) & !is.null(input$gwsBrowseScreenCellLineSelect)){
      preselCellline  <- contrasts %>%
        dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% input$gwsBrowseScreenCellLineSelect) %>%
        dplyr::select(cellline_name) %>%
        arrange(cellline_name) %>%
        .$cellline_name
    }
    
    #get selected library
    preselLibrary <- libraries %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue, cellline_name %in% preselCellline) %>%
      dplyr::select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsBrowseScreenCheckLibraryAll) & !is.null(input$gwsBrowseScreenLibrarySelect)){
      preselLibrary <- libraries %>%
        dplyr::filter(species %in% speciesList, 
               type %in% input$gwsBrowseScreenDatasetSelect, 
               tissue_name %in% preselTissue, 
               cellline_name %in% preselCellline,
               library_id %in% input$gwsBrowseScreenLibrarySelect) %>%
        dplyr::select(library_id) %>%
        .$library_id
    }
    
    presel_contrasts <- contrasts %>%
      dplyr::filter(species %in% speciesList,
             library_id %in% preselLibrary,
             tissue_name %in%  preselTissue,
             type %in% input$gwsBrowseScreenDatasetSelect,
             cellline_name %in% preselCellline
             ) 
    
    if("dropout" %in% input$gwsBrowseScreenDatasetSelect & length(input$gwsBrowseScreenDatasetSelect)) {
      presel_contrasts <- presel_contrasts %>%
        dplyr::filter(abs(dynamic_range) >= input$gwsBrowseScreenQuality)
    }
    
    presel_contrasts %>%
      dplyr::select(contrast_id) %>%
      distinct() %>%
      .$contrast_id
  }
})


# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------

observeEvent(input$gwsBrowseScreenLoadButton, {
  output$gwsBrowseScreenTable <- renderDataTable({
    gwsBrowseScreenDataTable()
  })
  output$gwsBrowseScreenContrastTable <- renderDataTable({
    gwsBrowseScreenContrastDataTable()
  })
  output$gwsBrowseScreenSampleTable <- renderDataTable({
    gwsBrowseScreenSampleDataTable()
  })
})

observeEvent(input$gwsBrowseScreenSearchRadio, {
  if(input$gwsBrowseScreenSearchRadio == "guide_id"){
    updateCheckboxInput(session, 'gwsBrowseScreenInclude', value="")
    disable("gwsBrowseScreenInclude")
    
    disable("gwsBrowseScreenIndexRadio")
    updateRadioButtons(session, 'gwsBrowseScreenIndexRadio', selected = "lfc")
    
    if("dropout" %in% input$gwsBrowseScreenDatasetSelect & length(input$gwsBrowseScreenDatasetSelect) == 1){
      enable("gwsBrowseScreenIndexRadio")
      updateRadioButtons(session,
                         "gwsBrowseScreenIndexRadio",
                         label = "Display data as:",
                         choices = list("Log-fold change" = "lfc", "Effect" = "effect_essentialome"),
                         selected = "lfc",
                         inline = F
                         )
    }
  }else{
    enable("gwsBrowseScreenInclude")
    
    disable("gwsBrowseScreenIndexRadio")
    updateRadioButtons(session, 'gwsBrowseScreenIndexRadio', selected = "lfc")
    
    if("dropout" %in% input$gwsBrowseScreenDatasetSelect & length(input$gwsBrowseScreenDatasetSelect) == 1){
      updateRadioButtons(session, "gwsBrowseScreenIndexRadio", 
                         choices = list("Log-fold change" = "lfc", "Effect" = "effect_essentialome", "FDR-adjusted effect Here" = "adjusted_effect_essentialome"), 
                         selected = "adjusted_effect_essentialome"
      )
    }
  }
})

observeEvent(input$gwsBrowseScreenSpeciesSelect, {
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsBrowseScreenSpeciesSelect
  }
  
  all_types <- contrasts %>%
    dplyr::filter(species %in% speciesList, !is.na(type)) %>%
    .$type %>%
    unique
  
  dataset_selection_all<- setNames(all_types, all_types)
  
  updateSelectizeInput(session, 'gwsBrowseScreenDatasetSelect', choices = dataset_selection_all, selected = "dropout", server = TRUE)
  #select checkbox tissue
  updateCheckboxInput(session, 'gwsBrowseScreenCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenTissueSelect', choices = gwsBrowseScreenTissueList(), server = TRUE)
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  #update other species select
  updateSpecies(input$gwsBrowseScreenSpeciesSelect)
  
  #disable laod button
  disable("gwsBrowseScreenLoadButton")
  gwsBrowseScreenUpdateText()
})

observeEvent(input$gwsBrowseScreenQuality, {
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  #disable laod button
  disable("gwsBrowseScreenLoadButton")
  gwsBrowseScreenUpdateText()
})

observeEvent(input$gwsBrowseScreenDatasetSelect, {
  #select checkbox tissue
  updateCheckboxInput(session, 'gwsBrowseScreenCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenTissueSelect', choices = gwsBrowseScreenTissueList(), server = TRUE)
  
  if(!is.null(input$gwsBrowseScreenDatasetSelect)){
    if("dropout" %in% input$gwsBrowseScreenDatasetSelect & length(input$gwsBrowseScreenDatasetSelect) == 1){
      
      updateRadioButtons(session, 'gwsBrowseScreenIndexRadio', selected = "adjusted_effect_essentialome")
      enable("gwsBrowseScreenIndexRadio")
      enable("gwsBrowseScreenQuality")
      updateRadioButtons(session, 'gwsBrowseScreenDisplayName', selected = "short")
      enable("gwsBrowseScreenDisplayName")
    }else{
      updateRadioButtons(session, 'gwsBrowseScreenIndexRadio', selected = "lfc")
      disable("gwsBrowseScreenIndexRadio")
      disable("gwsBrowseScreenQuality")
      updateRadioButtons(session, 'gwsBrowseScreenDisplayName', selected = "long")
      disable("gwsBrowseScreenDisplayName")
    }
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #update contrasts select
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenTissueSelect, {
  if(!is.null(input$gwsBrowseScreenTissueSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsBrowseScreenCheckTissueAll', value = FALSE)
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenCellLineSelect', choices = gwsBrowseScreenCellLineList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckCellLineAll', value = FALSE)
  
  #update contrasts selectb
  if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
    enable("gwsBrowseScreenCellLineSelect")
    enable("gwsBrowseScreenCheckCellLineAll")
  }else{
    disable("gwsBrowseScreenCellLineSelect")
    disable("gwsBrowseScreenCheckCellLineAll")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCheckTissueAll, {
  if(isTRUE(input$gwsBrowseScreenCheckTissueAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsBrowseScreenTissueSelect', choices = gwsBrowseScreenTissueList(), server = TRUE)
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenCellLineSelect', choices = gwsBrowseScreenCellLineList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckCellLineAll', value = FALSE)
  
  if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
    enable("gwsBrowseScreenCellLineSelect")
    enable("gwsBrowseScreenCheckCellLineAll")
  }else{
    disable("gwsBrowseScreenCellLineSelect")
    disable("gwsBrowseScreenCheckCellLineAll")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCellLineSelect, {
  if(!is.null(input$gwsBrowseScreenCellLineSelect)){
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsBrowseScreenCheckCellLineAll', value = FALSE)
  }
  #update cellline selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  
  #select cellline checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #update contrasts selectb
  if(isTRUE(input$gwsBrowseScreenCheckCellLineAll) | (!is.null(input$gwsBrowseScreenCellLineSelect))){
    enable("gwsBrowseScreenLibrarySelect")
    enable("gwsBrowseScreenCheckLibraryAll")
  }else{
    disable("gwsBrowseScreenLibrarySelect")
    disable("gwsBrowseScreenCheckLibraryAll")
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCheckCellLineAll, {
  if(isTRUE(input$gwsBrowseScreenCheckCellLineAll)){
    #reset tissue selectbox
    updateSelectizeInput(session, 'gwsBrowseScreenCellLineSelect', choices = gwsBrowseScreenCellLineList(), server = TRUE)
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
    enable("gwsBrowseScreenLibrarySelect")
    enable("gwsBrowseScreenCheckLibraryAll")
  }else{
    disable("gwsBrowseScreenLibrarySelect")
    disable("gwsBrowseScreenCheckLibraryAll")
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenLibrarySelect, {
  if(!is.null(input$gwsBrowseScreenLibrarySelect)){
    #unselect library checkbox
    updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  }
  
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  
  if((isTRUE(input$gwsBrowseScreenCheckLibraryAll) | (!is.null(input$gwsBrowseScreenLibrarySelect))) & (isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect)))) {
    enable("gwsBrowseScreenContrastSelect")
    enable("gwsBrowseScreenCheckContrastAll")
  }else{
    disable("gwsBrowseScreenContrastSelect")
    disable("gwsBrowseScreenCheckContrastAll")
  }
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCheckLibraryAll, {
  if(isTRUE(input$gwsBrowseScreenCheckLibraryAll)){
    #update library selectbox
    updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  }
  #update contrasts selectb
  updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = gwsBrowseScreenContrastList(), server = TRUE)
  #unselect checkbox contrasts
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  
  if((isTRUE(input$gwsBrowseScreenCheckLibraryAll) | (!is.null(input$gwsBrowseScreenLibrarySelect))) & (isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect)))) {
    enable("gwsBrowseScreenContrastSelect")
    enable("gwsBrowseScreenCheckContrastAll")
  }else{
    disable("gwsBrowseScreenContrastSelect")
    disable("gwsBrowseScreenCheckContrastAll")
  }
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenContrastSelect, {
  #unselect library checkbox
  if(!is.null(input$gwsBrowseScreenContrastSelect)){
    updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  }
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll) | (!is.null(input$gwsBrowseScreenContrastSelect))){
    enable("gwsBrowseScreenLoadButton")
  }else{
    disable("gwsBrowseScreenLoadButton")
  }
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCheckContrastAll, {
  #unselect library selectbox
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
    contrastList <- gwsBrowseScreenContrastList()
    updateSelectizeInput(session, 'gwsBrowseScreenContrastSelect', choices = contrastList, server = TRUE)
    
    showModal(modalDialog(
      title = "WARNING!", 
      paste0("WARNING: You have selected ", 
             length(contrastList), 
             " contrasts. Are you sure you want to load all contrasts for this selection?"),
      footer = tagList(
        modalButton("OK"),
        actionButton("gwsBrowseScreenCancelModal", "Cancel")
      )
    ))
  }
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll) | (!is.null(input$gwsBrowseScreenContrastSelect))){
    enable("gwsBrowseScreenLoadButton")
  }else{
    disable("gwsBrowseScreenLoadButton")
  }
  gwsBrowseScreenUpdateText()
  
}, ignoreNULL = FALSE)

observeEvent(input$gwsBrowseScreenCancelModal, {
  updateCheckboxInput(session, 'gwsBrowseScreenCheckContrastAll', value = FALSE)
  removeModal()
})

observeEvent(input$gwsBrowseScreenLoadDefaultTableModal, {
  removeModal()
  output$gwsBrowseScreenTable <- renderDataTable({
    
    if(!is.null(default_df)){
      df<-default_df
      
      vals<-df[,3:ncol(df)]
      values_min <- vals %>% min(na.rm = TRUE)
      values_max <- vals %>% max(na.rm = TRUE)
      
      brks_smaller <- seq(values_min, 0, .05)
      brks_bigger <- seq(0, values_max, .05)
      
      clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
      {paste0("rgb(255,", ., ",", ., ")")}
      clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
      {paste0("rgb(", ., ",", ., ",255)")}
      
      brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
      clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
      
      dt <- df %>%
        DT::datatable(extensions = c('FixedColumns','FixedHeader'), 
                      options = list(
                        autoWidth = FALSE,
                        headerCallback = JS(headerCallback),
                        scrollX=TRUE,
                        fixedColumns = list(leftColumns = 2),
                        columnDefs = list(list(className = 'dt-center', targets = "_all")),
                        pageLength = 25,
                        lengthMenu = c(25, 50, 100, 200),
                        searchHighlight = TRUE
                        #fixedHeader = TRUE
                      ),
                      filter = list(position = 'top', clear = FALSE),
                      rownames= FALSE) %>%
        formatStyle(seq(3, ncol(df),1),
                    backgroundColor = styleInterval(brks, clrs))
      
      #display default datatable
      dt 
    }
  })
})


# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$gwsBrowseScreenButtonDownload <- downloadHandler(
  filename = function() {
    "crisprepo_browse_screen.txt"
  },
  content = function(file) {
    
    displayed_table %>% write_tsv(file)
    
  }
)

output$gwsBrowseScreenButtonDownloadPrimaryTables <- downloadHandler(
  filename = function() {
    "screen_data.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading", input$dataset, " Data"),
      value = 0,
      {
        files <- NULL;
        if("Human" %in% input$gwsBrowseScreenDownloadPrimaryTablesCheck){
          fileName <- "essentiality_data/human_scaledLFC_geneLevel_HQscreens.tsv"
          files <- c(fileName,files)
        }
        if("Mouse" %in% input$gwsBrowseScreenDownloadPrimaryTablesCheck){
          fileName <- "essentiality_data/mouse_scaledLFC_geneLevel_HQscreens.tsv"
          files <- c(fileName,files)
        }
        shiny::incProgress(1/2)
        if(!is.null(files)){
          #create the zip file
          zip(file,files, compression_level = 2)
        }else{
          NULL
        }
      }
    )
  }
)

output$gwsBrowseScreenButtonDownloadAdjustedPrimaryTables <- downloadHandler(
  filename = function() {
    "screen_data_adjusted.zip"
  },
  content = function(file) {
    shiny::withProgress(
      message = paste0("Downloading", input$dataset, " Data"),
      value = 0,
      {
        files <- NULL;
        if("Human" %in% input$gwsBrowseScreenDownloadPrimaryTablesCheck){
          fileName <- "essentiality_data/human_scaledLFC_fdr_adjusted_geneLevel_HQscreens.tsv"
          files <- c(fileName,files)
        }
        if("Mouse" %in% input$gwsBrowseScreenDownloadPrimaryTablesCheck){
          fileName <- "essentiality_data/mouse_scaledLFC_fdr_adjusted_geneLevel_HQscreens.tsv"
          files <- c(fileName,files)
        }
        shiny::incProgress(1/2)
        if(!is.null(files)){
          #create the zip file
          zip(file,files, compression_level = 2)
        }else{
          NULL
        }
      }
    )
  }
)