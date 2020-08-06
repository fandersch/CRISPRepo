# crisprepo-shiny

# Copyright (c) 2018 Tobias Neumann, Jesse Lipp.
# 
# crisprepo-shiny is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# crisprepo-shiny is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

function(input, output, session) {
  ### Genome-wide Screens
  ####################
  # Browse Screen
  ####################

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
    
    presel = features %>%
      select(gene_id) %>% distinct %>% .$gene_id
    
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
        left_join(con_facs %>% tbl("features"), by = local(input$gwsBrowseScreenSearchRadio)) %>%
        select(contrast_id, contrast_id_QC, local(input$gwsBrowseScreenSearchRadio), lfc, effect, entrez_id, symbol) %>%
        distinct() %>%
        collect()
      
      remove_duplicates <- facs[facs %>% select(contrast_id, local(input$gwsBrowseScreenSearchRadio)) %>% duplicated,]
      
      df <- facs %>%
        anti_join(remove_duplicates)
      
    }else{
      df <- con %>%
        tbl(tableSelect) %>%
        filter(gene_id %in% presel, contrast_id %in% presel_contrasts) %>%
        left_join(con %>% tbl("contrasts"), by = "contrast_id") %>%
        left_join(con %>% tbl("features") %>% rename(symbol=hgnc_symbol), by = local(input$gwsBrowseScreenSearchRadio)) %>%
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
  
  #####################################
  # Browse Screen Selectbox Lists
  ####################################

  gwsBrowseScreenTissueList <- reactive({
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
  })
  
  
  #################################
  # Browse Screen Observers
  ################################
  
  observeEvent(input$gwsBrowseScreenLoadButton, {
    output$gwsBrowseScreenTable <- renderDataTable({
      gwsBrowseScreenDataTable()
    })
  })

  observeEvent(input$gwsBrowseScreenSpeciesSelect, {
    if(input$gwsBrowseScreenSpeciesSelect == "all"){
      updateSelectizeInput(session, 'gwsBrowseScreenDatasetSelect', choices = list("dropout" = "dropout", "drug_modifier" = "synthetic"), selected = "dropout", server = TRUE)
    }else{
      updateSelectizeInput(session, 'gwsBrowseScreenDatasetSelect', choices = list("dropout" = "dropout", "drug_modifier" = "synthetic", "facs_based" = "facs"), selected = "dropout", server = TRUE)
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
  
  ####################
  # Gene Search
  ####################
  
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
        left_join(features_buff, copy=TRUE, by = c("gene_id", "guide_id")) %>%
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
        left_join(features_buff %>% select(-guide_id, -sequence) %>% distinct, copy=TRUE, by ="gene_id") %>%
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
      Action <- shinyInput(actionButton, nrow(df), 'button_', label = "sgRNAs", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )
      df %>% cbind(Action)
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
          }else{
            nfreezeColumns <- 2
            nColorizeTableColumns <- 3
          }
          
          if(input$gwsGeneSpeciesSelect == "all"){
            nfreezeColumns <- nfreezeColumns + 2
            nColorizeTableColumns <- nColorizeTableColumns + 2
          }
          
          presel_contrasts <- colnames(gwsGeneDatatable)[nColorizeTableColumns:(length(colnames(gwsGeneDatatable))-1)]
          
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
                                                        lengthMenu = c(25, 50, 100, 200)),
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
  
  ####################################
  #Gene Search create selectbox lists
  ####################################
  
  gwsGeneTissueList <- reactive({
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
    if(input$gwsGeneSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$gwsGeneSpeciesSelect
    }
    
    preselLibrary = libraries %>%
      filter(species %in% speciesList) %>%
      select(library_id) %>%
      .$library_id
    
    if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
      preselLibrary = libraries %>%
        filter(species %in% speciesList) %>%
        filter(library_id %in% input$gwsGeneLibrarySelect) %>%
        .$library_id
    }
    
    preselTissue = pheno %>%
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
    
  })
  
  gwsGeneGeneList <- reactive({
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
      preselLibrary = libraries %>%
        filter(species %in% speciesList) %>%
        select(library_id) %>%
        .$library_id
      if(!isTRUE(input$gwsGeneCheckLibraryAll) & !is.null(input$gwsGeneLibrarySelect)){
        preselLibrary = libraries %>%
          filter(library_id %in% input$gwsGeneLibrarySelect) %>%
          collect() %>%
          .$library_id
      }
      features %>%
        filter(library_id %in% preselLibrary) %>% 
        select(symbol, entrez_id) %>%
        distinct() %>%
        dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (", entrez_id, ")"), paste0(symbol , " (", entrez_id, ")"))) %>%
        arrange(gene) %>%
        .$gene
    }
  })
  
  ###########################
  # Gene Search Observers
  ###########################
  
  observeEvent(input$gwsGeneLoadButton, {
    output$gwsGeneTable <- renderDataTable({
      gwsGeneDatatable <- gwsGeneDataTable()
    })
  })
  
  #actionButton handler for datatable buttons
  observeEvent(input$select_button, {
    row <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
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
    
    updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), selected = gene, server = TRUE)
    updateTabItems(session, "tabs", "sgRNAsSidebar")
  })
  
  observeEvent(input$gwsGeneSpeciesSelect, {
    if(input$gwsGeneSpeciesSelect == "all"){
      updateSelectizeInput(session, 'gwsGeneDatasetSelect', choices = list("dropout" = "dropout", "drug_modifier" = "synthetic"), selected = "dropout", server = TRUE)
    }else{
      updateSelectizeInput(session, 'gwsGeneDatasetSelect', choices = list("dropout" = "dropout", "drug_modifier" = "synthetic", "facs_based" = "facs"), selected = "dropout", server = TRUE)
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
  
  output$gwsGeneButtonDownload <- downloadHandler(
    filename = function() {
      table <- gwsGeneDataTable()
      paste0(paste(table$symbol %>% unique,collapse="_"), ".txt")
    },
    content = function(file) {
      table <- gwsGeneDataTable()
      table %>% select(-Action) %>%
        write_tsv(file)
    }
  )
  
  
  # ----------------------------------------------------------------------------
  # Libraries
  # ----------------------------------------------------------------------------
  
  libTable <- reactive({
    con %>% 
      tbl("features") %>% 
      filter(library_id %in% local(input$libSelectLibrary)) %>% 
      collect() %>%
      mutate(Length = nchar(sequence), order_entrezID = as.numeric(entrez_id), 
             sgRNA_23_mer = ifelse(!gene_id %in% c("AMBIGUOUS", "UNMAPPED", "NOFEATURE", "SAFETARGETING", "NONTARGETING"), paste0(gene_id, "_", substr(context, 5, nchar(context)-3)), NA)) %>%
      dplyr::arrange(ifelse(gene_id %in% c("AMBIGUOUS", "UNMAPPED", "NOFEATURE", "SAFETARGETING", "NONTARGETING"), 1, 0), order_entrezID, guide_id) %>%
      select("Guide-ID" = guide_id, "Entrez-ID" = entrez_id, "Gene-Symbol" = hgnc_symbol, "Ensembl-ID" = ensembl_id, Sequence = sequence, Length, 
             "sgRNA-ID-23-mer" = sgRNA_23_mer,
             "Library-ID" = library_id, Chromosome = chromosome, Strand = strand, "Genomic-Start-Position" = start, 
             "Genomic-End-Position" = end, "Perfect-matches-(PM)-total" = pm_total, "PM-with-PAM" = pm_pam, 
             "PM-with-PAM-in-CDS" = pm_pam_cds, "PM-with-PAM-in-CDS-unique" = pm_pam_cds_unique, 
             "30-mer-Genomic-Context" = context, "Legacy-ID" = legacy_id, "Legacy-Gene-Annotation" = legacy_group, Class = class)
      
      
  })
  
  output$libTableOutput <- renderDataTable({
    libTable()
  }, 
  filter = "bottom", 
  rownames= FALSE,
  extensions = c('FixedColumns','FixedHeader'),
  options = list(autoWidth = FALSE,
                 headerCallback = JS(headerCallback),
                 fixedColumns = list(leftColumns = 1),
                 scrollX=TRUE,
                 columnDefs = list(list(className = 'dt-center', targets = "_all")),
                 pageLength = 25,
                 lengthMenu = c(25, 50, 100, 200)
                 ))
  
  output$libBoxGuidesTotal <- renderInfoBox({
    infoBox(title = "Guides", value = libTable() %>% nrow())
  })
  
  output$libBoxGenesTotal <- renderInfoBox({
    infoBox(title = "Genes", value = libTable() %>% select("Entrez-ID") %>% distinct() %>% nrow())
  })
  
  output$libButtonDownload <- downloadHandler(
    filename = function() {
      paste0(local(input$libSelectLibrary), ".txt")
    },
    content = function(file) {
      libTable() %>% write_tsv(file)
    }
  )
  
  ####################
  # genome-wide sgRNA predictions
  ####################
  
  output$sgRNAsInfo <- renderText({
    invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
  })
  
  observeEvent(input$sgRNAsSpeciesSelect, {
    #update gene selectbox
    #update other species selects
    updateSelectizeInput(session, 'gwsGeneSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse"), selected = input$sgRNAsSpeciesSelect, server = TRUE)
    updateSelectizeInput(session, 'gwsBrowseScreenSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse"), selected = input$sgRNAsSpeciesSelect, server = TRUE)
    updateSelectizeInput(session, 'sgRNAsGeneSelect', choices = sgRNAsGeneList(), server = TRUE)
  })
  
  sgRNAsTable <- eventReactive(input$sgRNAsGeneSelect,{
    
    presel_genes<- input$sgRNAsGeneSelect %>% strsplit(split="\\(|\\)")
    presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
    presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
    
    #specify which sgRNAs table to join
    if(input$sgRNAsSpeciesSelect == "human"){
      tableSgRNAs <- "sgRNAs_human"
    }else{
      tableSgRNAs <- "sgRNAs_mouse"
    }
    
    
    if(input$sgRNAsSpeciesSelect == "all"){
      sgRNAs <- con_sgRNAs %>%
        tbl("sgRNAs_human") %>%
        filter(EntrezID %in% presel_gene_entrez) %>%
        select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
               `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
               Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
               `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
               `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
        collect() %>%
        rbind(
          con_sgRNAs %>%
            tbl("sgRNAs_mouse") %>%
            filter(EntrezID %in% presel_gene_entrez) %>%
            select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
                   `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
                   Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
                   `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
                   `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
            collect()
            
        )
    }else{
      sgRNAs <- con_sgRNAs %>%
        tbl(tableSgRNAs) %>%
        filter(EntrezID %in% presel_gene_entrez) %>% 
        select(`Entrez ID` = EntrezID, Symbol, `20-mer + NGG` = sgRNA_23mer, Position, Exon = exon, `Mature sgRNA` = mature_sgRNA,
               `VBC-Score`=VBC.score, `Frameshift ratio inDelphi` = inDelphi, `Cleavage activity` = cleavage_activity, Length = len_cloning_sgRNA, `Off-target predictions` = Off_target, 
               Length = len_cloning_sgRNA, `Prediction rank`= final_rank, `Validation rank` = final_validated_rank, `SNP targeting` = SNP_targeting,
               `Restriction-site- / Poly-A- / Multi-T-containing` = RS_PolyA_multiT_containing, `Failed validation / Single outlier validation` = failed_outlier_validation, 
               `Trans-species (human/mouse)` = transspecies, `Maps to genome` = check) %>%
        collect()
    }
    if(sgRNAs$`Maps to genome` %>% as.character %>% unique %>% length == "1"){
      sgRNAs <- sgRNAs %>% select(-`Maps to genome`)
    }
    sgRNAs
  })
  
  output$sgRNAsTableOutput <- renderDataTable({
    sgRNAs <- sgRNAsTable()
    if (nrow(sgRNAs) > 0) {
      sgRNAs %>% 
        datatable(extensions = c('FixedColumns','FixedHeader'),
          options = list(
          autoWidth = FALSE,
          headerCallback = JS(headerCallback),
          scrollX=TRUE,
          fixedColumns = list(leftColumns = 3),
          columnDefs = list(list(className = 'dt-center', targets = "_all")),
          pageLength = 25,
          lengthMenu = c(25, 50, 100, 200)
        ),
        filter = list(position = 'top', clear = FALSE),
        rownames= FALSE)
    }
  })
  
  sgRNAsGeneList <- reactive({
    if(input$sgRNAsSpeciesSelect == "human"){
      gene_list <- gene_list_human %>%
      dplyr::mutate(gene = ifelse(is.na(Symbol), paste0("No symbol found (", EntrezID, ")"), paste0(Symbol , " (", EntrezID, ")"))) %>%
        arrange(gene) %>%
        .$gene
    }else{
      if(input$sgRNAsSpeciesSelect == "mouse"){
        gene_list <- gene_list_mouse %>%
          dplyr::mutate(gene = ifelse(is.na(Symbol), paste0("No symbol found (", EntrezID, ")"), paste0(Symbol , " (", EntrezID, ")"))) %>%
          arrange(gene) %>%
          .$gene
      }else{
        gene_list <- gene_list_mouse %>% 
          rbind(gene_list_human) %>%
          dplyr::mutate(gene = ifelse(is.na(Symbol), paste0("No symbol found (", EntrezID, ")"), paste0(Symbol , " (", EntrezID, ")"))) %>%
          arrange(gene) %>%
          .$gene
      }
    }
    gene_list
  })
  
  output$sgRNAsButtonDownload <- downloadHandler(
    filename = function() {
      paste0(paste(input$sgRNAsGeneSelect,collapse="_"), ".txt")
    },
    content = function(file) {
      sgRNAsTable() %>% write_tsv(file)
    }
  )
  
  ####################
  # dual sgRNA designs
  ####################
  
  output$dualSgRNAsInfo <- renderText({
    invisible(paste("**************** TESTING MODE!!! THIS FEATURE IS STILL UNDER DEVELOPMENT, PLEASE DOUBLE CHECK YOUR RESULTS!!! **************", HTML('<br/>'), "INFO: Please upload a csv file with the provided file browser on the right side! ", HTML('<br/>'), " The file must have two columns: (1) the entrez-ID and (2) the 23 nucleotide long sgRNA sequence (no header)!"))
  })
  
  #upon loading do nothing
  output$dualSgRNAsTableOutput <- renderDataTable({
  })
  
  
  dualSgRNAsTable <- reactive({
    inFile <- input$dualSgRNAs_inputFile
    
    if(is.null(inFile)){
      return(NULL)
    }else{
      dualSgRNAs_input <- read_csv2(input$dualSgRNAs_inputFile$datapath, col_names =  FALSE, skip_empty_rows=TRUE)
      if(dualSgRNAs_input %>% ncol() != 2 | dualSgRNAs_input %>% nrow() == 0){
        dualSgRNAs_input <- read_csv(input$dualSgRNAs_inputFile$datapath, col_names =  FALSE, skip_empty_rows=TRUE)
        if(dualSgRNAs_input %>% ncol() != 2 | dualSgRNAs_input %>% nrow() == 0){
          return(NULL)
        }
      } 
    }

    entrez_list_human <- con_sgRNAs %>%
      tbl("sgRNAs_human") %>%
      select(EntrezID) %>%
      distinct %>%
      collect %>%
      .$EntrezID
    
    entrez_list_mouse <- con_sgRNAs %>%
      tbl("sgRNAs_mouse") %>%
      select(EntrezID) %>%
      distinct %>%
      collect() %>%
      .$EntrezID
  
    colnames(dualSgRNAs_input) <- c("entrezID", "sequence")
    
    dualSgRNAs_input <- dualSgRNAs_input %>%
      arrange(entrezID)
    
    dualSgRNAs_input$entrezID <- dualSgRNAs_input$entrezID %>% as.integer()
  
    entrez_old<-0
    dualSgRNAs_output <- NULL
    position_not_found_counter<-0
    entrez_not_found_counter<-0
    entrez_not_found<-c()
    position_not_found <- c()
    
    
    for(i in 1:nrow(dualSgRNAs_input)){
      
      input_sequence <- dualSgRNAs_input[i,2] %>% as.character
      input_entrez <- dualSgRNAs_input[i,1] %>% as.integer
      
      if(entrez_old!=input_entrez){
        if(input_entrez %in% entrez_list_human){
          sgRNA_candidates <- con_sgRNAs %>%
            tbl("sgRNAs_human") %>%
            filter(EntrezID == local(input_entrez)) %>%
            collect() %>%
            mutate_at(c("Symbol", "sgRNA_23mer", "sgRNA_ID", "Position", "mature_sgRNA", "Off_target", "guide_origin"), as.character) %>%
          
          entrez_old <- input_entrez
        }else{
          if(input_entrez %in% entrez_list_mouse){
            sgRNA_candidates <- con_sgRNAs %>%
              tbl("sgRNAs_mouse") %>%
              filter(EntrezID == local(input_entrez)) %>%
              collect() %>%
              mutate_at(c("Symbol", "sgRNA_23mer", "sgRNA_ID", "Position", "mature_sgRNA", "Off_target", "guide_origin"), as.character) %>%
              mutate(SNP_targeting = NA)
            
            entrez_old <- input_entrez
          }else{
            entrez_not_found <- c(entrez_not_found, input_entrez)
            entrez_not_found_counter<- entrez_not_found_counter+1
            position_not_found <- c(position_not_found, input_sequence)
            position_not_found_counter<- position_not_found_counter+1
            sgRNA_candidates <- NULL
          }
        }
      }
      
      if(!is.null(sgRNA_candidates)){
       
        #get position
        input_position <- sgRNA_candidates[sgRNA_candidates$sgRNA_23mer == input_sequence, "Position"]
        if(nrow(input_position)==1){
          input_position <- input_position %>% as.character
        }else{
          input_position <- NA
        }
        if(is.na(input_position) | input_position == ""){
          position_not_found_counter<-position_not_found_counter + 1
        }
        input_exon <- sgRNA_candidates[sgRNA_candidates$sgRNA_23mer == input_sequence, "exon"] %>% as.character()
        input_orientation <-  stringr::str_split(input_position, pattern = "[()]")[[1]][2]
        input_chr <-  stringr::str_split(input_position, pattern = "[-:(]")[[1]][1]
        input_start <-  ifelse(input_orientation=="+", stringr::str_split(input_position, pattern = "[-:(]")[[1]][2], stringr::str_split(input_position, pattern = "[-:(]")[[1]][3])
        input_end <-  ifelse(input_orientation=="+", stringr::str_split(input_position, pattern = "[-:(]")[[1]][3], stringr::str_split(input_position, pattern = "[-:(]")[[1]][2])
        
        
        # The 30nt include 4nt+23nt sgRNA + 3nt.
        # 3nt upstream of PAM
        input_genomic_cutting_position = ifelse(input_orientation=="+", as.numeric(input_end) - 3 - 3 - 3, as.numeric(input_end) + 3 + 3 + 3)

        sgRNAs_selected <- sgRNA_candidates %>%
          filter(!is.na(Position) & Position != "") %>%
          rowwise() %>%
          mutate(orientation =  stringr::str_split(`Position`, pattern = "[()]")[[1]][2],
                 chr = stringr::str_split(`Position`, pattern = "[-:(]")[[1]][1],
                 start =  ifelse(orientation=="+", stringr::str_split(`Position`, pattern = "[-:(]")[[1]][2], stringr::str_split(`Position`, pattern = "[-:(]")[[1]][3]),
                 end =  ifelse(orientation=="+", stringr::str_split(`Position`, pattern = "[-:(]")[[1]][3], stringr::str_split(`Position`, pattern = "[-:(]")[[1]][2])) %>%
          mutate(genomic_cutting_position = ifelse(orientation=="+", as.numeric(end) - 3 - 3 - 3, as.numeric(end) + 3 + 3 + 3)) %>%
          mutate(cutting_distance = input_genomic_cutting_position - genomic_cutting_position,
                 produces_frameshift = ifelse(cutting_distance %% 3 != 0, TRUE, FALSE),
                 proximity_100kb = ifelse(abs(cutting_distance)<=1000, TRUE, FALSE),
                 targets_same_exon = ifelse(input_exon==exon, TRUE, FALSE)) %>%
          # filter(proximity_100kb == TRUE, produces_frameshift == TRUE, check == TRUE) %>%
          mutate(original_sgRNA = input_sequence, original_sgRNA_position = input_position, original_sgRNA_exon_number = input_exon) %>%
          mutate(VBC.score = ifelse(VBC.score <= 0.001, NA, VBC.score)) %>%
          arrange(final_rank) %>%
          select(EntrezID, Symbol, original_sgRNA, original_sgRNA_position, original_sgRNA_exon_number,  matching_sgRNA=sgRNA_23mer, Position, VBC.score, Off_target, cutting_distance, exon, targets_same_exon, everything()) 
        
        #Limit number of reported dual-sgRNA-combinations per gene
        if(input$dualSgRNAs_LimitOutput == TRUE){
          limit <- ifelse(nrow(sgRNAs_selected) > input$dualSgRNAs_nOutput, input$dualSgRNAs_nOutput, nrow(sgRNAs_selected))
          sgRNAs_selected <- sgRNAs_selected[1:limit,]
        }
        if(i == 1){
          dualSgRNAs_output <- sgRNAs_selected
        }else{
          dualSgRNAs_output <- dualSgRNAs_output %>% rbind(sgRNAs_selected)
        }
      }
    }
    
    if(entrez_not_found_counter>0 | position_not_found_counter > 0){
      if(entrez_not_found_counter ==1){
        text_entrezID <- paste0(entrez_not_found_counter, " entrez-ID could not be found!<br>", "(", entrez_not_found, ")<br>")
      }else{
        if(entrez_not_found_counter >0){
          text_entrezID <- paste0(entrez_not_found_counter, " entrez-IDs could not be found!<br>", "(", paste(entrez_not_found, sep="", collapse = ", "), ")<br>")
        }else{
          text_entrezID <- ""
        }
      }
      if(position_not_found_counter ==1){
        text_position <- paste0(position_not_found_counter, " sgRNA genomic position could not be obtained by the provided 23-mer sgRNA sequence!", "(", position_not_found, ")<br>")
      }else{
        if(position_not_found_counter >0){
          text_position <- paste0(position_not_found_counter, " sgRNA genomic positions could not be obtained by the provided 23-mer sgRNA sequence!", "(", paste(position_not_found, sep="", collapse = ", "), ")<br>")
        }else{
          text_position <- ""
        }
      }
      showModal(modalDialog(
        title = "WARNING!", 
        HTML(paste0("WARNING:<br>", text_entrezID, text_position)),
        footer = tagList(
          modalButton("OK")
        )
      ))
    }
    dualSgRNAs_output 
  })
  
  dualSgRNAsDataTable <- eventReactive(input$dualSgRNALoadButton,{
   
    inFile <- input$dualSgRNAs_inputFile

    if (is.null(inFile)){
      return(NULL)
    }else{
        dualSgRNAs_input <- read_csv2(input$dualSgRNAs_inputFile$datapath, col_names =  FALSE, skip_empty_rows=TRUE)
      if(dualSgRNAs_input %>% ncol() != 2 | dualSgRNAs_input %>% nrow() == 0){
        dualSgRNAs_input <- read_csv(input$dualSgRNAs_inputFile$datapath, col_names =  FALSE, skip_empty_rows=TRUE)
        if(dualSgRNAs_input %>% ncol() != 2 | dualSgRNAs_input %>% nrow() == 0){
          showModal(modalDialog(
            title = "WARNING!", 
            HTML(paste0("WARNING:<br> Your uploaded file does not have two columns seperated by either ';' or ',' or has no lines. Please select another file!")),
            footer = tagList(
              modalButton("OK")
            )
          ))
          return(NULL)
        }
      }

      dualSgRNAs_output <- dualSgRNAsTable()
      if (nrow(dualSgRNAs_output) > 0) {
        return(
          dualSgRNAs_output %>%
            datatable(options = list(
              autoWidth = FALSE,
              headerCallback = JS(headerCallback),
              scrollX=TRUE,
              columnDefs = list(list(className = 'dt-center', targets = "_all")),
              pageLength = 25,
              lengthMenu = c(25, 50, 100, 200)
            ),
            filter = list(position = 'top', clear = FALSE),
            rownames= FALSE)
        )
      }
    }
  })
  
  observeEvent(input$dualSgRNALoadButton, {
    output$dualSgRNAsTableOutput <- renderDataTable({
      dualSgRNAsDataTable()
    })
  })
  
  observeEvent(input$dualSgRNAs_LimitOutput, {
    if(isTRUE(input$dualSgRNAs_LimitOutput)){
      enable("dualSgRNAs_nOutput")
    }else{
      disable("dualSgRNAs_nOutput")
    }
    
  }, ignoreNULL = FALSE)
  
  output$dualSgRNAsButtonDownload <- downloadHandler(
    filename = function() {
      paste0("dual_sgRNA_designs", ".txt")
    },
    content = function(file) {
      dualSgRNAsTable() %>% write_tsv(file)
    }
  )
  
  
  # ----------------------------------------------------------------------------
  # ExpressionData
  # ----------------------------------------------------------------------------
  
  
  expressionDataUpdateText <- function(){
    output$expressionDataInfo <- renderText({
      if(is.null(input$expressionDataTissueSelect) & !isTRUE(input$expressionDataCheckTissueAll)){
        invisible("INFO: Please select the tissue(s) in the right panel!")
      }else{
        if(is.null(input$expressionDataCellLineSelect) & !isTRUE(input$expressionDataCheckCellLineAll)){
          invisible("INFO: Please select the cell line(s) in the right panel!")
        }else{
          if(is.null(input$expressionDataGeneSelect) & !isTRUE(input$expressionDataCheckGeneAll)){
            invisible("INFO: Please select the gene(s) you want to browse in the right panel!")
          }else{
            invisible("INFO: Click Load data!")
          }
        }
      }
    })
  }
  
  #upon load display nothing
  output$expressionDataTable <- renderDataTable({
  })
  
  #query database and create dataframe
  expressionDataDataFrame <- reactive({
    #get selected species
    if(input$expressionDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$expressionDataSpeciesSelect
    }
    
    #get selected tissue
    if(isTRUE(input$expressionDataCheckTissueAll)){
      presel_tissue <- expressionDataTissueList()
    }else{
      presel_tissue <- local(input$expressionDataTissueSelect)
    }
    
    #get selected cell line
    if(isTRUE(input$expressionDataCheckCellLineAll)){
      presel_cell_line <- expressionDataCellLineList()
    }else{
      presel_cell_line <- local(input$expressionDataCellLineSelect)
    }
    
    if(isTRUE(input$expressionDataCheckGeneAll)){
      presel_genes <- expressionDataGeneList()
    }else{
      presel_genes <- local(input$expressionDataGeneSelect)
    }
    
    presel_genes<- presel_genes %>% strsplit(split="\\(|\\)")
    presel_gene_symbol <- unlist(presel_genes)[c(TRUE, FALSE)] %>% trimws()
    presel_gene_entrez <- unlist(presel_genes)[c(FALSE, TRUE)] %>% as.numeric
    
    sample_ids <- cellline_list_expressionData %>%
      filter(species %in% speciesList) %>%
      filter(tissue_name %in%  presel_tissue) %>%
      filter(cell_line_name %in% presel_cell_line) %>%
      select(sample_id) %>%
      .$sample_id
    
    df <- con_expression %>%
      tbl("expression_data_values") %>%
      select(sample_id, gene_symbol, entrez_id, expression_value) %>%
      filter(entrez_id %in% presel_gene_entrez, gene_symbol %in% presel_gene_symbol, sample_id %in% sample_ids) %>%
      distinct() %>%
      collect()
    
    df <- df %>%
      left_join(cellline_list_expressionData %>% select(sample_id, species, unit))

    if(input$expressionDataSpeciesSelect == "all"){
      df_human <- df %>%
        filter(species == "human") %>%
        left_join(dict_joined %>% select(Symbol_human, Symbol_mouse, EntrezID_mouse) %>% filter(!is.na(Symbol_human)), by=c("gene_symbol" = "Symbol_human")) %>%
        dplyr::rename(Symbol_human = gene_symbol, EntrezID_human = entrez_id)
      
      df_mouse <- df %>%
        filter(species == "mouse") %>%
        left_join(dict_joined %>% select(Symbol_human, Symbol_mouse, EntrezID_human) %>% filter(!is.na(Symbol_mouse)), by=c("gene_symbol" = "Symbol_mouse")) %>%
        dplyr::rename(Symbol_mouse = gene_symbol, EntrezID_mouse = entrez_id)
      
      df <- df_human %>% rbind(df_mouse)
    }

    df
    
  })
  
  #create datatable out of dataframe
  expressionDataDataTable <- eventReactive(input$expressionDataLoadButton,{
    
    df <- expressionDataDataFrame()
        if (nrow(df) > 0) {
          
      # # color codig for heatmap
      # values <- df %>%
      #   select(input$gwsBrowseScreenIndexRadio) %>%
      #   as.data.frame() %>% as.matrix() %>% as.vector()
      # 
      # brks_smaller <- seq(min(values), 0, .05)
      # brks_bigger <- seq(0, max(values), .05)
      # 
      # clrs_smaller <- round(seq(40, 255, length.out = (length(brks_smaller) + 1)), 0) %>%
      # {paste0("rgb(255,", ., ",", ., ")")}
      # clrs_bigger <- round(seq(255, 40, length.out = (length(brks_bigger))), 0) %>%
      # {paste0("rgb(", ., ",", ., ",255)")}
      # 
      # brks <- c(as.vector(brks_smaller), as.vector(brks_bigger))
      # clrs <- c(as.vector(clrs_smaller), as.vector(clrs_bigger))
      # 
      nfreezeColumns <- 2
      
      if(input$expressionDataSpeciesSelect == "all"){
        
        dt <- df %>%
          select(sample_id, expression_value, Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, unit) %>%
          spread(sample_id, expression_value) %>%
          arrange(Symbol_human, Symbol_mouse) %>%
          select(Symbol_human, EntrezID_human, Symbol_mouse, EntrezID_mouse, everything())

        nfreezeColumns <- nfreezeColumns + 2
      }else{
        dt <- df %>%
          select(sample_id, gene_symbol, entrez_id, expression_value, unit) %>%
          spread(sample_id, expression_value) %>%
          arrange(gene_symbol)
      }

      # #tooltips
      # colnames_dt <- colnames(dt)
      # contrast_ids <- df$contrast_id %>% unique
      # tooltip <- ''
      # 
      # if(input$gwsBrowseScreenDatasetSelect %in% c("dropout")){
      #   for(i in 1:length(colnames_dt)){
      #     if(i < length(colnames_dt)){
      #       tooltip <- paste0(tooltip, "'", colnames_dt[i], "'",  ", " )
      #     }else{
      #       tooltip <- paste0(tooltip, "'", colnames_dt[i], "'")
      #     }
      #     if(colnames_dt[i] %in% contrast_ids){
      #       colnames_dt[i] <- contrasts %>% select(contrast_id, contrast_id_QC) %>% filter(contrast_id == colnames_dt[i]) %>% .$contrast_id_QC
      #     }
      #   }
      #   colnames(dt) <- colnames_dt
      # }

      # headerCallback <- c(
      #   "function(thead, data, start, end, display){",
      #   "  var $ths = $(thead).find('th');",
      #   "  $ths.css({'vertical-align': 'bottom', 'white-space': 'nowrap'});",
      #   "  var betterCells = [];",
      #   "  $ths.each(function(){",
      #   "    var cell = $(this);",
      #   "    var newDiv = $('<div>', {height: 'auto', width: 'auto'});",
      #   "    var newInnerDiv = $('<div>', {text: cell.text()});",
      #   "    newDiv.css({margin: 'auto'});",
      #   "    newInnerDiv.css({",
      #   "      transform: 'rotate(180deg)',",
      #   "      'writing-mode': 'tb-rl',",
      #   "      'white-space': 'nowrap'",
      #   "    });",
      #   "    newDiv.append(newInnerDiv);",
      #   "    betterCells.push(newDiv);",
      #   "  });",
      #   "  $ths.each(function(i){",
      #   "    $(this).html(betterCells[i]);",
      #   "  });",
      #   paste0("  var tooltips = [", tooltip, "];"),
      #   paste0("  for(var i=0; i<", length(colnames_dt), "; i++){"),
      #   "    $('th:eq('+i+')',thead).attr('title', tooltips[i]);",
      #   "  }",
      #   "}")

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
                      rownames= FALSE)
        # formatStyle(seq(nfreezeColumns+1, length(colnames_dt),1),
        #             backgroundColor = styleInterval(brks, clrs))
        
      if(!is.null(input$expressionDataGeneSelect)){
        output$expressionDataInfo <- renderText({
          "Info: Loading completed!"
        })
      }
      #display datatable
      dt
    }else{
      if(!is.null(input$expressionDataGeneSelect)){
        output$expressionDataInfo <- renderText({
          "WARNING: No data found!"
        })
      }else{
        expressionDataUpdateText()
      }
    }
  })
  
  #####################################
  # Expression Data Selectbox Lists
  ####################################
  
  expressionDataTissueList <- reactive({
    if(input$expressionDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$expressionDataSpeciesSelect
    }
    
    cellline_list_expressionData %>%
      filter(species %in% speciesList) %>%
      select(tissue_name) %>%
      distinct() %>%
      arrange(tissue_name) %>%
      .$tissue_name
    
  })
  
  expressionDataCellLineList <- reactive({
    if(input$expressionDataSpeciesSelect == "all"){
      speciesList <- c("human", "mouse")
    }else{
      speciesList <- input$expressionDataSpeciesSelect
    }
    
    preselTissue = cellline_list_expressionData %>%
      filter(species %in%  speciesList) %>%
      select(tissue_name) %>%
      distinct() %>%
      .$tissue_name
    
    if(!isTRUE(input$expressionDataCheckTissueAll) & !is.null(input$expressionDataTissueSelect)){
      preselTissue = cellline_list_expressionData %>%
        filter(species %in% speciesList) %>%
        filter(tissue_name %in% input$expressionDataTissueSelect) %>%
        select(tissue_name) %>%
        distinct() %>%
        .$tissue_name
    }
    
    cellline_list_expressionData %>%
      filter(species %in% speciesList, tissue_name %in% preselTissue) %>%
      select(cell_line_name) %>%
      arrange(cell_line_name) %>%
      .$cell_line_name
  })
  
  expressionDataGeneList <- reactive({
    if(!isTRUE(input$expressionDataCheckCellLineAll) & (is.null(input$expressionDataCellLineSelect))){
      NULL
    }else{
      if(input$expressionDataSpeciesSelect == "all"){
        speciesList <- c("human", "mouse")
      }else{
        speciesList <- input$expressionDataSpeciesSelect
      }
      
      gene_list_expressionData %>%
        filter(species %in% speciesList) %>%
        select(gene_symbol, entrez_id) %>%
        collect() %>%
        dplyr::mutate(gene = ifelse(is.na(gene_symbol), paste0("No symbol found (", entrez_id, ")"), paste0(gene_symbol , " (", entrez_id, ")"))) %>%
        arrange(gene) %>%
        .$gene
    }
  })
  
  
  #################################
  # Browse Screen Observers
  ################################
  
  observeEvent(input$expressionDataLoadButton, {
    output$expressionDataTable <- renderDataTable({
      expressionDataDataTable()
    })
  })
  
  observeEvent(input$expressionDataSpeciesSelect, {
    #select checkbox tissue
    updateCheckboxInput(session, 'expressionDataCheckTissueAll', value = FALSE)
    #update tissue selectbox
    updateSelectizeInput(session, 'expressionDataTissueSelect', choices = expressionDataTissueList(), server = TRUE)
    #update cell linne selectbox
    updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
    #select cell line checkbox
    updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
    #update contrasts selectb
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
    #unselect checkbox contrasts
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    #disable laod button
    disable("expressionDataLoadButton")
    expressionDataUpdateText()
  })
  
  observeEvent(input$expressionDataTissueSelect, {
    if(!is.null(input$expressionDataTissueSelect)){
      #unselect checkbox tissue
      updateCheckboxInput(session, 'expressionDataCheckTissueAll', value = FALSE)
    }
    
    #update cell line selectbox
    updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
    #select cell line checkbox
    updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
    #update cell line selectb
    if(isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect))){
      enable("expressionDataCellLineSelect")
      enable("expressionDataCheckCellLineAll")
    }else{
      disable("expressionDataCellLineSelect")
      disable("expressionDataCheckCellLineAll")
    }
    
    #update gene selectb
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
    #unselect checkbox contrasts
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  observeEvent(input$expressionDataCheckTissueAll, {
    if(isTRUE(input$expressionDataCheckTissueAll)){
      #reset tissue selectbox
      updateSelectizeInput(session, 'expressionDataTissueSelect', choices = expressionDataTissueList(), server = TRUE)
    }
    
    #update cell line selectbox
    updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
    #select cell line checkbox
    updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
    #update cell line selectb
    if(isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect))){
      enable("expressionDataCellLineSelect")
      enable("expressionDataCheckCellLineAll")
    }else{
      disable("expressionDataCellLineSelect")
      disable("expressionDataCheckCellLineAll")
    }
    
    #update gene selectb
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
    #unselect checkbox contrasts
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  observeEvent(input$expressionDataCellLineSelect, {
    if(!is.null(input$expressionDataCellLineSelect)){
      #unselect library checkbox
      updateCheckboxInput(session, 'expressionDataCheckCellLineAll', value = FALSE)
    }
    
    #update gene selectb
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
    #unselect checkbox gene
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    
    if((isTRUE(input$expressionDataCheckCellLineAll) | (!is.null(input$expressionDataCellLineSelect))) & (isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect)))) {
      enable("expressionDataGeneSelect")
      enable("expressionDataCheckGeneAll")
    }else{
      disable("expressionDataGeneSelect")
      disable("expressionDataCheckGeneAll")
    }
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  observeEvent(input$expressionDataCheckCellLineAll, {
    if(isTRUE(input$expressionDataCheckCellLineAll)){
      #update library selectbox
      updateSelectizeInput(session, 'expressionDataCellLineSelect', choices = expressionDataCellLineList(), server = TRUE)
    }
    
    #update gene selectb
    updateSelectizeInput(session, 'expressionDataGeneSelect', choices = expressionDataGeneList(), server = TRUE)
    #unselect checkbox gene
    updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    
    if((isTRUE(input$expressionDataCheckCellLineAll) | (!is.null(input$expressionDataCellLineSelect))) & (isTRUE(input$expressionDataCheckTissueAll) | (!is.null(input$expressionDataTissueSelect)))) {
      enable("expressionDataGeneSelect")
      enable("expressionDataCheckGeneAll")
    }else{
      disable("expressionDataGeneSelect")
      disable("expressionDataCheckGeneAll")
    }
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  observeEvent(input$expressionDataGeneSelect, {
    #unselect library checkbox
    if(!is.null(input$expressionDataGeneSelect)){
      updateCheckboxInput(session, 'expressionDataCheckGeneAll', value = FALSE)
    }
    if(isTRUE(input$expressionDataCheckGeneAll) | (!is.null(input$expressionDataGeneSelect))){
      enable("expressionDataLoadButton")
    }else{
      disable("expressionDataLoadButton")
    }
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  observeEvent(input$expressionDataCheckGeneAll, {
    #unselect library selectbox
    if(isTRUE(input$expressionDataCheckGeneAll)){
      geneList <- expressionDataGeneList()
      updateSelectizeInput(session, 'expressionDataGeneSelect', choices = geneList, server = TRUE)
    }
    if(isTRUE(input$expressionDataCheckGeneAll) | (!is.null(input$expressionDataGeneSelect))){
      enable("expressionDataLoadButton")
    }else{
      disable("expressionDataLoadButton")
    }
    expressionDataUpdateText()
    
  }, ignoreNULL = FALSE)
  
  output$expressionDataButtonDownload <- downloadHandler(
    filename = function() {
      if(isTRUE(input$expressionDataCheckGeneAll)){
        paste0(paste(c("expresssion_all_genes_", local(input$expressionDataTissueSelect), local(input$expressionDataCellLineSelect)),  collapse="_"), ".txt")
      }else{
        str_replace_all(string = paste0("expression_", paste(local(input$expressionDataGeneSelect),collapse="_"), ".txt"), pattern = " ", replacement = "_")
      }
    },
    content = function(file) {
      df <- expressionDataDataFrame()
      
      if (nrow(df) > 0) {
        dt <- df %>%
          select(sample_id, gene_symbol, entrez_id, expression_value, unit) %>%
          spread(sample_id, expression_value) %>%
          arrange(gene_symbol) %>% write_tsv(file)
      }
    }
  )
  
  
##################
# Header callback
##################
  
  #rotate vertical
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
    "}"
  )
  
}