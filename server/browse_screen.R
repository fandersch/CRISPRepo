# ----------------------------------------------------------------------------
# Browse Screen
# ----------------------------------------------------------------------------

gwsBrowseScreenUpdateText <- function(){
  output$gwsBrowseScreenInfo <- renderText({
    if(is.null(input$gwsBrowseScreenTissueSelect) & !isTRUE(input$gwsBrowseScreenCheckTissueAll)){
      invisible("INFO: Please select the tissue(s) of considered screens in the right panel!")
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
  })
}

#upon load display nothing
output$gwsBrowseScreenTable <- renderDataTable({
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
  
  #query database
  if(input$gwsBrowseScreenDatasetSelect %in% c("facs")){
    contrasts_buff <- contrasts_facs
    
    facs <- con_facs %>%
      tbl(tableSelect) %>%
      filter(contrast_id %in% presel_contrasts) %>%
      left_join(con_facs %>% tbl("contrasts"), by = "contrast_id") %>%
      left_join(con_facs %>% tbl("features") %>% select(symbol, entrez_id) %>% distinct, by = local(input$gwsBrowseScreenSearchRadio)) %>%
      select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, entrez_id, symbol) %>%
      distinct() %>%
      collect()
    
    remove_duplicates <- facs[facs %>% select(contrast_id, local(input$gwsBrowseScreenSearchRadio)) %>% duplicated,]
    
    df <- facs %>%
      anti_join(remove_duplicates)
    
  }else{
    df <- con %>%
      tbl(tableSelect) %>%
      filter(contrast_id %in% presel_contrasts) %>%
      filter(gene_id != "AMBIGUOUS") %>%
      filter(gene_id != "UNMAPPED") %>%
      filter(gene_id != "NOFEATURE") %>%
      filter(gene_id != "SAFETARGETING") %>%
      filter(gene_id != "NONTARGETING") %>%
      filter(!is.na(gene_id)) %>%
      left_join(con %>% tbl("contrasts"), by = "contrast_id") %>%
      left_join(con %>% tbl("features") %>% select(local(input$gwsBrowseScreenSearchRadio), hgnc_symbol, entrez_id, species) %>% distinct %>% rename(symbol=hgnc_symbol), by = local(input$gwsBrowseScreenSearchRadio)) %>%
      select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, symbol, entrez_id, species) %>%
      distinct() %>%
      collect()
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      df_human <- df %>% 
        filter(species == "human") %>%
        left_join(dict_joined, by=c("symbol" = "Symbol_human")) %>%
        mutate(EntrezID_human = entrez_id) %>%
        select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, Symbol_human = symbol, EntrezID_human, Symbol_mouse,  EntrezID_mouse)
      
      df_mouse <- df %>% 
        filter(species == "mouse") %>%
        left_join(dict_joined, by=c("symbol" = "Symbol_mouse")) %>%
        mutate(EntrezID_mouse = entrez_id) %>%
        select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, Symbol_human, EntrezID_human, Symbol_mouse = symbol, EntrezID_mouse)
      
      df <- df_human %>% rbind(df_mouse)
    }
    df
  }
  
})

#create datatable out of dataframe
gwsBrowseScreenDataTable <- eventReactive(input$gwsBrowseScreenLoadButton,{
  
  df <- gwsBrowseScreenDataFrame()
  if (nrow(df) > 0) {
    # color codig for heatmap
    values <- df %>% 
      select(input$gwsBrowseScreenIndexRadio) %>%
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
        select(contrast_id, local(input$gwsBrowseScreenSearchRadio), Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, input$gwsBrowseScreenIndexRadio) %>%
        select(-contains("gene_id")) %>%
        spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
        arrange(Symbol_human, Symbol_mouse) %>%
        select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse=EntrezID_mouse, everything())
      
      nfreezeColumns <- nfreezeColumns + 2
    }else{
      dt <- df %>%
        select(contrast_id, local(input$gwsBrowseScreenSearchRadio), entrez_id, symbol, input$gwsBrowseScreenIndexRadio) %>%
        spread(contrast_id, input$gwsBrowseScreenIndexRadio) %>%
        select(-contains("gene_id")) %>%
        arrange(symbol)
    }
    
    colnames_dt <- colnames(dt)
    contrast_ids <- df$contrast_id %>% unique
    tooltip <- ''
    
    if(input$gwsBrowseScreenDatasetSelect %in% c("dropout")){
      for(i in 1:length(colnames_dt)){
        if(i < length(colnames_dt)){
          tooltip <- paste0(tooltip, "'", colnames_dt[i], "'",  ", " )
        }else{
          tooltip <- paste0(tooltip, "'", colnames_dt[i], "'")
        }
        if(colnames_dt[i] %in% contrast_ids){
          colnames_dt[i] <- contrasts %>% select(contrast_id, contrast_id_QC) %>% filter(contrast_id == colnames_dt[i]) %>% .$contrast_id_QC
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
                      lengthMenu = c(25, 50, 100, 200)
                      #fixedHeader = TRUE
                    ),
                    filter = list(position = 'top', clear = FALSE),
                    rownames= FALSE) %>%
      formatStyle(seq(nfreezeColumns+1, length(colnames_dt),1),
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
    
    contrasts_facs <<- contrasts_facs %>%
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
    filter(species %in% speciesList) %>%
    select(tissue_name) %>%
    distinct() %>%
    arrange(tissue_name) %>%
    .$tissue_name
  
})

gwsBrowseScreenLibraryList <- reactive({
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$gwsBrowseScreenSpeciesSelect
  }
  libraries %>%
    filter(species %in% speciesList) %>%
    select(library_id) %>%
    arrange(library_id) %>%
    .$library_id
})

gwsBrowseScreenContrastList <- reactive({
  if(!is.null(input$gwsBrowseScreenLibrarySelect) | isTRUE(input$gwsBrowseScreenCheckLibraryAll)){
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsBrowseScreenSpeciesSelect
    }
    preselLibrary = libraries %>%
      filter(species %in% speciesList) %>%
      select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsBrowseScreenCheckLibraryAll) & !is.null(input$gwsBrowseScreenLibrarySelect)){
      preselLibrary = libraries %>%
        filter(species %in% speciesList) %>%
        filter(library_id %in% input$gwsBrowseScreenLibrarySelect) %>%
        .$library_id
    }
    
    preselTissue = pheno %>%
      filter(species %in%  speciesList) %>%
      select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$gwsBrowseScreenCheckTissueAll) & !is.null(input$gwsBrowseScreenTissueSelect)){
      preselTissue = pheno %>%
        filter(species %in% speciesList) %>%
        filter(tissue_name %in% input$gwsBrowseScreenTissueSelect) %>%
        .$tissue_name
    }
    if(input$gwsBrowseScreenDatasetSelect %in% c("facs")){
      contrasts_buff <- contrasts_facs
      if("human" ==  input$gwsBrowseScreenSpeciesSelect){
        preselLibrary <- "zuber_library_original"
      }else{
        if("human" ==  input$gwsBrowseScreenSpeciesSelect){
          preselLibrary <- "zuber_library_mouse_original"
        }else{
          preselLibrary <- c("zuber_library_original", "zuber_library_mouse_original")
        }
      }
    }else{
      contrasts_buff <- contrasts
    }
    contrasts_buff %>%
      filter(species %in% speciesList) %>%
      filter(library_id %in% preselLibrary) %>%
      filter(tissue_name %in%  preselTissue) %>%
      filter(type == input$gwsBrowseScreenDatasetSelect) %>%
      select(contrast_id) %>%
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
})

observeEvent(input$gwsBrowseScreenSpeciesSelect, {
  if(input$gwsBrowseScreenSpeciesSelect == "all"){
    updateSelectizeInput(session, 'gwsBrowseScreenDatasetSelect', choices = dataset_selection_dropout_drug, selected = "dropout", server = TRUE)
  }else{
    updateSelectizeInput(session, 'gwsBrowseScreenDatasetSelect', choices = dataset_selection_all, selected = "dropout", server = TRUE)
  }
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

observeEvent(input$gwsBrowseScreenDatasetSelect, {
  #select checkbox tissue
  updateCheckboxInput(session, 'gwsBrowseScreenCheckTissueAll', value = FALSE)
  #update tissue selectbox
  updateSelectizeInput(session, 'gwsBrowseScreenTissueSelect', choices = gwsBrowseScreenTissueList(), server = TRUE)
  
  if(input$gwsBrowseScreenDatasetSelect %in% c("synthetic", "facs")){
    if("human" ==  input$gwsBrowseScreenSpeciesSelect){
      #update library selectbox
      updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), selected = "hs_gw_zuber_v2", server = TRUE)
    }else{
      if("mouse" ==  input$gwsBrowseScreenSpeciesSelect){
        #update library selectbox
        updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), selected = "mm_gw_zuber_v1", server = TRUE)
      }else{
        #update library selectbox
        updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), selected = c("mm_gw_zuber_v1", "hs_gw_zuber_v2"), server = TRUE)
        
      }
    }
    #unselect checkbox tissue
    updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
    disable("gwsBrowseScreenLibrarySelect")
    disable("gwsBrowseScreenCheckLibraryAll")
    #disable index
    updateRadioButtons(session, 'gwsBrowseScreenIndexRadio',selected = "lfc")
    disable("gwsBrowseScreenIndexRadio")
  }else{
    #update library selectbox
    updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
    #select library checkbox
    updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
    #enable index
    enable("gwsBrowseScreenIndexRadio")
  }
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
  if(input$gwsBrowseScreenDatasetSelect %in% c("synthetic", "facs")){
    if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
      enable("gwsBrowseScreenContrastSelect")
      enable("gwsBrowseScreenCheckContrastAll")
    }else{
      disable("gwsBrowseScreenContrastSelect")
      disable("gwsBrowseScreenCheckContrastAll")
    }
  }else{
    #update library selectbox
    updateSelectizeInput(session, 'gwsBrowseScreenLibrarySelect', choices = gwsBrowseScreenLibraryList(), server = TRUE)
    #select library checkbox
    updateCheckboxInput(session, 'gwsBrowseScreenCheckLibraryAll', value = FALSE)
    #update contrasts selectb
    if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
      enable("gwsBrowseScreenLibrarySelect")
      enable("gwsBrowseScreenCheckLibraryAll")
    }else{
      disable("gwsBrowseScreenLibrarySelect")
      disable("gwsBrowseScreenCheckLibraryAll")
    }
  }
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
  if(input$gwsBrowseScreenDatasetSelect %in% c("synthetic", "facs")){
    if(isTRUE(input$gwsBrowseScreenCheckTissueAll) | (!is.null(input$gwsBrowseScreenTissueSelect))){
      enable("gwsBrowseScreenContrastSelect")
      enable("gwsBrowseScreenCheckContrastAll")
    }else{
      disable("gwsBrowseScreenContrastSelect")
      disable("gwsBrowseScreenCheckContrastAll")
    }
  }else{
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
      paste0(paste(c("all_contrasts", local(input$gwsBrowseScreenTissueSelect), local(input$gwsBrowseScreenLibrarySelect)),  collapse="_"), ".txt")
    }else{
      paste0(paste(local(input$gwsBrowseScreenContrastSelect),collapse="_"), ".txt")
    }
  },
  content = function(file) {
    df <- gwsBrowseScreenDataFrame()
    
    if (nrow(df) > 0) {
      df %>%
        select(contrast_id, local(input$gwsBrowseScreenSearchRadio), entrez_id, symbol, local(input$gwsBrowseScreenIndexRadio)) %>%
        spread(contrast_id, local(input$gwsBrowseScreenIndexRadio)) %>%
        #use function arrange_at() instead of arrange() because input$gwsBrowseScreenSearchRadio needs to be parsed toquosure, then unquoted it  !!
        arrange_at(local(input$gwsBrowseScreenSearchRadio)) %>% write_tsv(file)
    }
  }
)