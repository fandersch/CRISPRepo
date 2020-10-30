# ----------------------------------------------------------------------------
# dual sgRNA designs
# ----------------------------------------------------------------------------

output$dualSgRNAsInfo <- renderText({
  invisible(paste("INFO: Please upload a csv file with the provided file browser on the right side! ", HTML('<br/>'), " The file must have two columns: (1) the entrez-ID and (2) the 23 nucleotide long sgRNA sequence (no header)!"))
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
          mutate_at(c("Symbol", "sgRNA_23mer", "sgRNA_ID", "Position", "mature_sgRNA", "Off_target", "guide_origin"), as.character)
          
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
               produces_frameshift = ifelse(cutting_distance %% 3 != 0 & input_exon==exon, TRUE, ifelse(cutting_distance %% 3 != 0 & input_exon!=exon), NA, FALSE),
               proximity_1kb = ifelse(abs(cutting_distance)<=1000, TRUE, FALSE),
               targets_same_exon = ifelse(input_exon==exon, TRUE, FALSE)) %>%
        filter(proximity_1kb == TRUE, is.na(produces_frameshift), produces_frameshift == TRUE) %>%
        mutate(original_sgRNA = input_sequence, original_sgRNA_position = input_position, original_sgRNA_exon_number = input_exon) %>%
        mutate(VBC.score = ifelse(VBC.score <= 0.001, NA, VBC.score)) %>%
        arrange(EntrezID, desc(proximity_1kb), desc(produces_frameshift), final_rank) %>%
        select(EntrezID, Symbol, original_sgRNA, original_sgRNA_position, original_sgRNA_exon_number,  matching_sgRNA=sgRNA_23mer, Position, VBC.score, Off_target, cutting_distance, exon, targets_same_exon, everything()) %>%
        rename(maps_to_genome = check)

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
  
  if(length(sgRNAs_selected$maps_to_genome %>% unique)==1){
    sgRNAs_selected <- sgRNAs_selected %>%
      select(-maps_to_genome)
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
            lengthMenu = c(25, 50, 100, 200),
            searchHighlight = TRUE
          ),
          filter = list(position = 'top', clear = FALSE),
          rownames= FALSE)
      )
    }
  }
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$dualSgRNAsButtonDownload <- downloadHandler(
  filename = function() {
    paste0("dual_sgRNA_designs", ".txt")
  },
  content = function(file) {
    dualSgRNAsTable() %>% write_tsv(file)
  }
)
