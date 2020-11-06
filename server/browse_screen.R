# ----------------------------------------------------------------------------
# Browse Screen
# ----------------------------------------------------------------------------

gwsBrowseScreenUpdateText <- function(){
  output$gwsBrowseScreenInfo <- renderText({
    if(is.null(input$gwsBrowseScreenTissueSelect) & !isTRUE(input$gwsBrowseScreenCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) of considered screens in the right panel!")
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
  
  df <- con %>%
    tbl("contrasts") %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    distinct() %>%
    collect()
  
  namekey <- c(contrast_id="contrast_id", contrast_id_QC="contrast_id_QC", treatment="treatment", control="control", type="type", cellline_name="cellline_name", 
               tissue_name= "tissue_name", tissue_context = "tissue_context", library_id="library_id", species="species", 
               auc= "auc", ssmd = "ssmd", dynamic_range = "dynamic_range", average_lfc_essentials = "average_lfc_essentials", 
               norm_method="norm_method_mageck", fdr_method="fdr_method_mageck", lfc_method="lfc_method_mageck", cnv_correction="cnv_correction_mageck",
               filter="filter_mageck", variance_estimation="variance_estimation_mageck")
  
  names(df) <- namekey[names(df)]
  df <- df %>%
    dplyr::select(contrast_id, contrast_id_QC, treatment, control, type, cellline_name, tissue_name, tissue_context, library_id, species, 
           auc, ssmd, dynamic_range, average_lfc_essentials, 
           matches("norm_method_mageck"), matches("fdr_method_mageck"), matches("lfc_method_mageck"), matches("cnv_correction_mageck"), 
           matches("filter_mageck"), matches("variance_estimation_mageck"))
  
  if(input$gwsBrowseScreenDisplayName == "short"){
    df <- df %>%
      dplyr::select(contrast_id_QC, contrast_id, everything())
  }
  
  if(!input$gwsBrowseScreenDatasetSelect %in% c("dropout")){
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
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
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
                      pageLength = 25,
                      lengthMenu = c(25, 50, 100, 200),
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
  }else{
    tableSelect <- "gene_stats"
  }
  
  if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
    presel_contrasts <- gwsBrowseScreenContrastList()
  }else{
    presel_contrasts <- local(input$gwsBrowseScreenContrastSelect)
  }
  
  statistics_columns_negative <- paste0(local(input$gwsBrowseScreenInclude), "_negative")
  statistics_columns_positive <- paste0(local(input$gwsBrowseScreenInclude), "_positive")
  
  df <- con %>%
    tbl(tableSelect) %>%
    dplyr::filter(contrast_id %in% presel_contrasts) %>%
    dplyr::filter(gene_id != "AMBIGUOUS") %>%
    dplyr::filter(gene_id != "UNMAPPED") %>%
    dplyr::filter(gene_id != "NOFEATURE") %>%
    dplyr::filter(gene_id != "SAFETARGETING") %>%
    dplyr::filter(gene_id != "NONTARGETING") %>%
    dplyr::filter(!is.na(gene_id)) %>%
    left_join(con %>% tbl("contrasts") %>% dplyr::select(contrast_id, contrast_id_QC, species, library_id), by = "contrast_id") %>%
    left_join(con %>% tbl("features") %>% dplyr::select(local(input$gwsBrowseScreenSearchRadio), hgnc_symbol, entrez_id, library_id) %>% distinct %>% rename(symbol=hgnc_symbol), by = c(local(input$gwsBrowseScreenSearchRadio), "library_id")) %>%
    dplyr::select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, symbol, entrez_id, species, dplyr::one_of(statistics_columns_negative), dplyr::one_of(statistics_columns_positive)) %>%
    distinct() %>%
    collect()
  
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    
    df_human <- df %>% 
      dplyr::filter(species == "human") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_human")) %>%
      mutate(EntrezID_human = entrez_id) %>%
      dplyr::select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, Symbol_human = symbol, EntrezID_human, Symbol_mouse,  EntrezID_mouse, dplyr::one_of(statistics_columns_negative), dplyr::one_of(statistics_columns_positive))
    
    df_mouse <- df %>% 
      dplyr::filter(species == "mouse") %>%
      left_join(dict_joined, by=c("symbol" = "Symbol_mouse")) %>%
      mutate(EntrezID_mouse = entrez_id) %>%
      dplyr::select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, Symbol_human, EntrezID_human, Symbol_mouse = symbol, EntrezID_mouse, dplyr::one_of(statistics_columns_negative), dplyr::one_of(statistics_columns_positive))
    
    df <- df_human %>% rbind(df_mouse)
  }
  df

})

#create datatable out of dataframe
gwsBrowseScreenDataTable <- eventReactive(input$gwsBrowseScreenLoadButton,{
  
  df <- gwsBrowseScreenDataFrame()
  if (nrow(df) > 0) {
    # color codig for heatmap
    values <- df %>% 
      dplyr::select(input$gwsBrowseScreenIndexRadio) %>%
      as.data.frame() %>% as.matrix() %>% as.vector()
    
    brks_smaller <- seq(min(values,na.rm = TRUE), 0, .05)
    brks_bigger <- seq(0, max(values,na.rm = TRUE), .05)
    
    clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}
    clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
    {paste0("rgb(", ., ",", ., ",255)")}
    
    brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
    clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
    
    if(input$gwsBrowseScreenSearchRadio == "guide_id"){
      nfreezeColumns <- 3
    }else{
      nfreezeColumns <- 2
    }
    
    df$lfc <- round(df$lfc, 3)
    df$effect <- round(df$effect, 3)

    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      
      dt <- df %>%
        dplyr::select(contrast_id, local(input$gwsBrowseScreenSearchRadio), Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, input$gwsBrowseScreenIndexRadio) %>%
        dplyr::select(-contains("gene_id")) %>%
        spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
        arrange(Symbol_human, Symbol_mouse) %>%
        dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse=EntrezID_mouse, everything())
      
      nfreezeColumns <- nfreezeColumns + 2
    }else{
      dt <- df %>%
        dplyr::select(contrast_id, local(input$gwsBrowseScreenSearchRadio), entrez_id, symbol, input$gwsBrowseScreenIndexRadio) %>%
        spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
        dplyr::select(-contains("gene_id")) %>%
        arrange(symbol)
    }
    if(input$gwsBrowseScreenSearchRadio == "gene_id" & !is.null(local(input$gwsBrowseScreenInclude)) & length(local(input$gwsBrowseScreenInclude))!=0){
      column <- local(input$gwsBrowseScreenInclude)
      if("p"  %in% local(input$gwsBrowseScreenInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(p = ifelse(p_positive < p_negative, p_positive, p_negative)) %>%
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), p) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="p") %>%
              mutate_if(is.numeric, round, 3) %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_P'))
            )
      }
      if("fdr" %in% local(input$gwsBrowseScreenInclude)){
        dt <- dt %>%
          left_join(
            df %>%
              mutate(fdr = ifelse(fdr_positive < fdr_negative, fdr_positive, fdr_negative)) %>%
              dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), fdr) %>%
              distinct %>%
              pivot_wider(names_from="contrast_id", values_from="fdr") %>%
              mutate_if(is.numeric, round, 3) %>%
              rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_FDR'))
          )
      }
      if("guides" %in% local(input$gwsBrowseScreenInclude)){
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
      if("guides_good" %in% local(input$gwsBrowseScreenInclude)){
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
    }
    
    dt <- dt %>%
      dplyr::select(sort(current_vars())) %>%
      dplyr::select(contains(local(input$gwsBrowseScreenSearchRadio)), contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), everything())
    
    colorInterval <- length(local(input$gwsBrowseScreenInclude))
    #change to contrast_id_QC for selected columns when type is dropout
    presel_contrasts <- gwsBrowseScreenContrastList()
    colnames_dt <- colnames(dt)
    tooltip <- ''
    if(input$gwsBrowseScreenDatasetSelect %in% c("dropout") & input$gwsBrowseScreenDisplayName == "short"){
      for(i in 1:length(colnames_dt)){
        if(i < length(colnames_dt)){
          tooltip <- paste0(tooltip, "'", colnames_dt[i], "'",  ", " )
        }else{
          tooltip <- paste0(tooltip, "'", colnames_dt[i], "'")
        }
        colname_buff <- colnames_dt[i]
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
          colnames_dt[i] <- str_replace(colnames_dt[i], colname_buff, (contrasts %>% dplyr::select(contrast_id, contrast_id_QC) %>% dplyr::filter(contrast_id == colname_buff) %>% .$contrast_id_QC))
        }
      }
      
      colnames(dt) <- colnames_dt
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
      paste0("  for(var i=0; i<", length(colnames_dt), "; i++){"),
      "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
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
  
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsBrowseScreenSpeciesSelect
  }
  
  pheno %>%
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
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    pheno %>%
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
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      distinct %>%
      .$cellline_name
      
    if(!isTRUE(input$gwsBrowseScreenCheckCellLineAll) & !is.null(input$gwsBrowseScreenCellLineSelect)){
      preselCellline  <- pheno %>%
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
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsBrowseScreenSpeciesSelect
    }
    
    #get selected tissue
    preselTissue <- pheno %>%
      dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect) %>%
      dplyr::select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue <- pheno %>%
        dplyr::filter(species %in%  speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        dplyr::select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    #get selected cell line
    preselCellline <- pheno %>%
      dplyr::filter(species %in% speciesList, type %in% input$gwsBrowseScreenDatasetSelect, tissue_name %in% preselTissue) %>%
      dplyr::select(cellline_name) %>%
      arrange(cellline_name) %>%
      .$cellline_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckCellLineAll) & !is.null(input$gwsBrowseScreenCellLineSelect)){
      preselCellline  <- pheno %>%
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
             type == input$gwsBrowseScreenDatasetSelect,
             cellline_name %in% preselCellline
             ) 
    
    if(input$gwsBrowseScreenDatasetSelect == "dropout"){
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
    disable("gwsBrowseScreenInclude")
  }else{
    enable("gwsBrowseScreenInclude")
  }
})

observeEvent(input$gwsBrowseScreenSpeciesSelect, {
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
  updateSelectizeInput(session, 'sgRNAsSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = input$gwsBrowseScreenSpeciesSelect, server = TRUE)
  updateSelectizeInput(session, 'gwsGeneSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = input$gwsBrowseScreenSpeciesSelect, server = TRUE)
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
  
  if(input$gwsBrowseScreenDatasetSelect %in% c("drug_modifier", "facs_based")){
    updateRadioButtons(session, 'gwsBrowseScreenIndexRadio',selected = "lfc")
    disable("gwsBrowseScreenIndexRadio")
    disable("gwsBrowseScreenQuality")
  }else{
    enable("gwsBrowseScreenIndexRadio")
    enable("gwsBrowseScreenQuality")
  }
  #update library selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
  #select library checkbox
  updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
  #enable index
  enable("gwsBrowseScreenIndexRadio")
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
             " contrasts. Are you sure you want to load all contrasts for this selection? This might take a long time to load!"),
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

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$gwsBrowseScreenButtonDownload <- downloadHandler(
  filename = function() {
    if(isTRUE(input$gwsBrowseScreenCheckContrastAll)){
      paste0(paste(c("all_contrasts", local(input$gwsBrowseScreenTissueSelect), local(input$gwsBrowseScreenCellLineSelect), local(input$gwsBrowseScreenLibrarySelect)),  collapse="_"), ".txt")
    }else{
      paste0(paste(local(input$gwsBrowseScreenContrastSelect),collapse="_"), ".txt")
    }
  },
  content = function(file) {
    df <- gwsBrowseScreenDataFrame()
    
    if (nrow(df) > 0) {
      if(input$gwsBrowseScreenSpeciesSelect == "all"){
        
        dt <- df %>%
          dplyr::select(contrast_id, local(input$gwsBrowseScreenSearchRadio), Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, input$gwsBrowseScreenIndexRadio) %>%
          dplyr::select(-contains("gene_id")) %>%
          spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
          arrange(Symbol_human, Symbol_mouse) %>%
          dplyr::select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse=EntrezID_mouse, everything())
        
        nfreezeColumns <- nfreezeColumns + 2
      }else{
        dt <- df %>%
          dplyr::select(contrast_id, local(input$gwsBrowseScreenSearchRadio), entrez_id, symbol, input$gwsBrowseScreenIndexRadio) %>%
          spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
          dplyr::select(-contains("gene_id")) %>%
          arrange(symbol)
      }
      if(input$gwsBrowseScreenSearchRadio == "gene_id" & !is.null(local(input$gwsBrowseScreenInclude)) & length(local(input$gwsBrowseScreenInclude))!=0){
        column <- local(input$gwsBrowseScreenInclude)
        if("p"  %in% local(input$gwsBrowseScreenInclude)){
          dt <- dt %>%
            left_join(
              df %>%
                mutate(p = ifelse(p_positive < p_negative, p_positive, p_negative)) %>%
                dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), p) %>%
                distinct %>%
                pivot_wider(names_from="contrast_id", values_from="p") %>%
                mutate_if(is.numeric, round, 3) %>%
                rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_P'))
            )
        }
        if("fdr" %in% local(input$gwsBrowseScreenInclude)){
          dt <- dt %>%
            left_join(
              df %>%
                mutate(fdr = ifelse(fdr_positive < fdr_negative, fdr_positive, fdr_negative)) %>%
                dplyr::select(contrast_id, contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), fdr) %>%
                distinct %>%
                pivot_wider(names_from="contrast_id", values_from="fdr") %>%
                mutate_if(is.numeric, round, 3) %>%
                rename_at(vars(-contains("entrez_id"), -contains("symbol"), -contains("Symbol_human"), -contains("EntrezID_human"), -contains("Symbol_mouse"), -contains("EntrezID_mouse")), ~ paste0(., '_FDR'))
            )
        }
        if("guides" %in% local(input$gwsBrowseScreenInclude)){
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
        if("guides_good" %in% local(input$gwsBrowseScreenInclude)){
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
      }
      
      dt <- dt %>%
        dplyr::select(sort(current_vars())) %>%
        dplyr::select(contains(local(input$gwsBrowseScreenSearchRadio)), contains("entrez_id"), contains("symbol"), contains("Symbol_human"), contains("EntrezID_human"), contains("Symbol_mouse"), contains("EntrezID_mouse"), everything())
      
      presel_contrasts <- gwsBrowseScreenContrastList()
      colnames_dt <- colnames(dt)
      tooltip <- ''
      if(input$gwsBrowseScreenDatasetSelect %in% c("dropout") & input$gwsBrowseScreenDisplayName == "short"){
        for(i in 1:length(colnames_dt)){
          colname_buff <- colnames_dt[i]
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
            colnames_dt[i] <- str_replace(colnames_dt[i], colname_buff, (contrasts %>% dplyr::select(contrast_id, contrast_id_QC) %>% dplyr::filter(contrast_id == colname_buff) %>% .$contrast_id_QC))
          }
        }
        
        colnames(dt) <- colnames_dt
      }
      
      dt %>% write_tsv(file)
    }
  }
)