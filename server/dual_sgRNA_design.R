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
    dplyr::select(EntrezID) %>%
    distinct %>%
    collect %>%
    .$EntrezID
  
  entrez_list_mouse <- con_sgRNAs %>%
    tbl("sgRNAs_mouse") %>%
    dplyr::select(EntrezID) %>%
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
  species <- NULL
  
  
  for(i in 1:nrow(dualSgRNAs_input)){
    
    firstSgRNA_sequence <- dualSgRNAs_input[i,2] %>% as.character
    firstSgRNA_entrez <- dualSgRNAs_input[i,1] %>% as.integer
    
    #get all sgRNAs for this entrezID (only if entrezID changes)
    if(entrez_old!=firstSgRNA_entrez){
      if(firstSgRNA_entrez %in% entrez_list_human){
        sgRNA_candidates <- con_sgRNAs %>%
          tbl("sgRNAs_human") %>%
          dplyr::filter(EntrezID == local(firstSgRNA_entrez)) %>%
          collect() %>%
          mutate_at(c("Symbol", "sgRNA_23mer", "sgRNA_ID", "Position", "mature_sgRNA", "Off_target", "guide_origin"), as.character)
          
          entrez_old <- firstSgRNA_entrez
          species <- "human"
      }else{
        if(firstSgRNA_entrez %in% entrez_list_mouse){
          sgRNA_candidates <- con_sgRNAs %>%
            tbl("sgRNAs_mouse") %>%
            dplyr::filter(EntrezID == local(firstSgRNA_entrez)) %>%
            collect() %>%
            mutate_at(c("Symbol", "sgRNA_23mer", "sgRNA_ID", "Position", "mature_sgRNA", "Off_target", "guide_origin"), as.character) %>%
            mutate(SNP_targeting = NA)
          
          entrez_old <- firstSgRNA_entrez
          species <- "mouse"
        }else{
          entrez_not_found <- c(entrez_not_found, firstSgRNA_entrez)
          entrez_not_found_counter<- entrez_not_found_counter+1
          position_not_found <- c(position_not_found, firstSgRNA_sequence)
          position_not_found_counter<- position_not_found_counter+1
          sgRNA_candidates <- NULL
          species <- NULL
        }
      }
    }
    
    if(!is.null(sgRNA_candidates) & !is.null(species)){
      
      #get position from matching 23mer
      firstSgRNA_position <- sgRNA_candidates[sgRNA_candidates$sgRNA_23mer == firstSgRNA_sequence, "Position"]
      #must retrieve only one matching sgRNA position
      if(nrow(firstSgRNA_position)==1){
        firstSgRNA_position <- firstSgRNA_position %>% as.character
      }else{
        firstSgRNA_position <- NA
      }
      if(is.na(firstSgRNA_position) | firstSgRNA_position == ""){
        position_not_found_counter<-position_not_found_counter + 1
      }else{
        
        firstSgRNA_exon <- sgRNA_candidates[sgRNA_candidates$sgRNA_23mer == firstSgRNA_sequence, "exon"] %>% as.character()
        firstSgRNA_orientation <-  stringr::str_split(firstSgRNA_position, pattern = "[()]")[[1]][2]
        firstSgRNA_chr <-  stringr::str_split(firstSgRNA_position, pattern = "[-:(]")[[1]][1]
        firstSgRNA_start <-  ifelse(firstSgRNA_orientation=="+", stringr::str_split(firstSgRNA_position, pattern = "[-:(]")[[1]][2], stringr::str_split(firstSgRNA_position, pattern = "[-:(]")[[1]][3])
        firstSgRNA_end <-  ifelse(firstSgRNA_orientation=="+", stringr::str_split(firstSgRNA_position, pattern = "[-:(]")[[1]][3], stringr::str_split(firstSgRNA_position, pattern = "[-:(]")[[1]][2])
        # The 30nt include 4nt+23nt sgRNA + 3nt. Cutting site is described as 3nt upstream of PAM
        firstSgRNA_genomic_cutting_position = ifelse(firstSgRNA_orientation=="+", as.numeric(firstSgRNA_end) - 3 - 3 - 3, as.numeric(firstSgRNA_end) + 3 + 3 + 3)
        
        sgRNAs_selected <- sgRNA_candidates %>%
          dplyr::filter(!is.na(Position) & Position != "") %>%
          rowwise() %>%
          mutate(second_sgRNA_orientation =  stringr::str_split(`Position`, pattern = "[()]")[[1]][2],
                 second_sgRNA_chr = stringr::str_split(`Position`, pattern = "[-:(]")[[1]][1],
                 second_sgRNA_start =  ifelse(second_sgRNA_orientation=="+", stringr::str_split(`Position`, pattern = "[-:(]")[[1]][2], stringr::str_split(`Position`, pattern = "[-:(]")[[1]][3]),
                 second_sgRNA_end =  ifelse(second_sgRNA_orientation=="+", stringr::str_split(`Position`, pattern = "[-:(]")[[1]][3], stringr::str_split(`Position`, pattern = "[-:(]")[[1]][2]), 
                 second_sgRNA_genomic_cutting_position = ifelse(second_sgRNA_orientation=="+", as.numeric(second_sgRNA_end) - 3 - 3 - 3, as.numeric(second_sgRNA_end) + 3 + 3 + 3),
                 genomic_cutting_distance = firstSgRNA_genomic_cutting_position - second_sgRNA_genomic_cutting_position,
                 produces_frameshift = ifelse(genomic_cutting_distance %% 3 != 0 & firstSgRNA_exon==exon, TRUE, ifelse(genomic_cutting_distance %% 3 == 0 & firstSgRNA_exon==exon, FALSE, NA)),
                 proximity_1kb = ifelse(abs(genomic_cutting_distance)<=1000, TRUE, FALSE),
                 targets_same_exon = ifelse(firstSgRNA_exon==exon, TRUE, FALSE),
                 first_sgRNA_sequence_23mer = firstSgRNA_sequence, 
                 first_sgRNA_position_30mer = firstSgRNA_position, 
                 first_sgRNA_exon_number = firstSgRNA_exon,
                 first_sgRNA_genomic_cutting_position = firstSgRNA_genomic_cutting_position,
                 VBC.score = ifelse(VBC.score <= 0.001, NA, VBC.score)) %>%
          dplyr::filter(proximity_1kb == TRUE, check == TRUE) %>%
          arrange(EntrezID, desc(proximity_1kb), desc(produces_frameshift), final_rank) %>%
          dplyr::select(EntrezID, Symbol, 
                 first_sgRNA_sequence_23mer, first_sgRNA_position_30mer, first_sgRNA_exon_number, first_sgRNA_genomic_cutting_position,
                 second_sgRNA_sequence_23mer=sgRNA_23mer, second_sgRNA_mature_sgRNA=mature_sgRNA, second_sgRNA_position_30mer=Position, second_sgRNA_exon=exon, second_sgRNA_genomic_cutting_position, 
                 genomic_cutting_distance, targets_same_exon, proximity_1kb, produces_frameshift, 
                 VBC.score, Off_target, inDelphi, cleavage_activity, everything()) %>%
          dplyr::select(-sgRNA_ID, -guide_origin, -second_sgRNA_orientation, -second_sgRNA_chr, -second_sgRNA_start, -second_sgRNA_end, -check)
        
        chr_buff <-""
        x<-1
        for(x in 1:nrow(sgRNAs_selected)){
          
          sgRNAs_selected[x,"produces_frameshift_in_most_transcripts"] <- NA
          sgRNAs_selected[x,"produces_frameshift_in_all_transcripts"] <- NA
          sgRNAs_selected[x,"most_common_dist"] <- NaN
          sgRNAs_selected[x,"first_sgRNA_targeted_transcripts"] <- ""
          sgRNAs_selected[x,"nFirst_sgRNA_targeted_transcripts"] <- NaN
          sgRNAs_selected[x,"second_sgRNA_targeted_transcripts"] <- ""
          sgRNAs_selected[x,"nSecond_sgRNA_targeted_transcripts"] <- NaN
          sgRNAs_selected[x,"TranscriptsInFrame"] <- ""
          sgRNAs_selected[x,"nTranscriptsInFrame"] <- NaN
          sgRNAs_selected[x,"TranscriptsOutOfFrame"] <- ""
          sgRNAs_selected[x,"nTranscriptsOutOfFrame"] <- NaN

          if(is.na(sgRNAs_selected[x,"produces_frameshift"]) & sgRNAs_selected[x,"first_sgRNA_sequence_23mer"] != sgRNAs_selected[x,"second_sgRNA_sequence_23mer"]){
            
            input_position <- sgRNAs_selected$first_sgRNA_position_30mer[x]
            input_orientation <-  stringr::str_split(input_position, pattern = "[()]")[[1]][2]
            input_chr <-  str_replace(stringr::str_split(input_position, pattern = "[-:(]")[[1]][1], pattern = "chr", replacement = "")
            input_start <-  ifelse(input_orientation=="+", stringr::str_split(input_position, pattern = "[-:(]")[[1]][2], stringr::str_split(input_position, pattern = "[-:(]")[[1]][3])
            input_end <-  ifelse(input_orientation=="+", stringr::str_split(input_position, pattern = "[-:(]")[[1]][3], stringr::str_split(input_position, pattern = "[-:(]")[[1]][2])
            input_genomic_cutting_position = ifelse(input_orientation=="+", as.numeric(input_end) - 3 - 3 - 3, as.numeric(input_end) + 3 + 3 + 3)
            
            combination_cutting_site <- sgRNAs_selected$second_sgRNA_genomic_cutting_position[x]
            
            
            #get EnsDb object from current chromosome (only if chromosome changes)
            if(chr_buff != input_chr){
              chr_buff <- input_chr
              if(species=="human"){
                edbx <- ensembldb::filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == input_chr)
              }else{
                edbx <- ensembldb::filter(EnsDb.Mmusculus.v79, filter = ~ seq_name == input_chr)
              }
            }
            
            #make GRanges object with both sgRNA cutting sites
            gnm <- c(GRanges(paste0(input_chr, ":", input_genomic_cutting_position, "-", input_genomic_cutting_position+1)), GRanges(paste0(input_chr, ":", combination_cutting_site, "-", combination_cutting_site+1)))
            
            ## Map genomic coordinates to within-transcript coordinates
            gnm_tx <- genomeToTranscript(gnm, edbx)
            
            #get within-transcript distances from transcript pairs
            res <- gnm_tx %>% as.data.frame() %>%
              group_by(tx_id) %>%
              mutate(distance = start - last(start))
            
            dist_tx <- res %>% 
              ungroup %>%
              dplyr::filter(distance !=0) %>%
              dplyr::select(distance, names)
            
            if(!is.null(dist_tx) & nrow(dist_tx) >= 1){
            
              first_sgRNA_transcripts <- res[res$group==1,"names"] %>% .$names
              nfirst_sgRNA_transcripts <- first_sgRNA_transcripts %>% length
              second_sgRNA_transcripts <- res[res$group==2,"names"] %>% .$names
              nsecond_sgRNA_transcripts <- second_sgRNA_transcripts %>% length
              most_common_dist <- sort(table(dist_tx$distance),decreasing=TRUE)[1] %>% names %>% as.numeric
              tx_frames <- dist_tx %>% mutate(frame = distance %% 3)
              TranscriptsInFrame <- tx_frames %>% dplyr::filter(frame == 0) %>% .$names
              nTranscriptsInFrame <- TranscriptsInFrame %>% length
              TranscriptsOutOfFrame <- tx_frames %>% dplyr::filter(frame != 0) %>% .$names
              nTranscriptsOutOfFrame <- TranscriptsOutOfFrame %>% length
              
              #if the most common distance is not in frame
              if(most_common_dist %% 3 !=0){
                sgRNAs_selected[x,"produces_frameshift_in_most_transcripts"] <- TRUE
              }else{
                sgRNAs_selected[x,"produces_frameshift_in_most_transcripts"] <- FALSE
              }
              #if all distances are not in frame
              if(nTranscriptsInFrame == 0){
                sgRNAs_selected[x,"produces_frameshift"] <- TRUE
                sgRNAs_selected[x,"produces_frameshift_in_all_transcripts"] <- TRUE
              }else{
                sgRNAs_selected[x,"produces_frameshift"] <- FALSE
                sgRNAs_selected[x,"produces_frameshift_in_all_transcripts"] <- FALSE
              }
              
              sgRNAs_selected[x,"most_common_dist"] <- most_common_dist
              sgRNAs_selected[x,"first_sgRNA_targeted_transcripts"] <- paste(first_sgRNA_transcripts, collapse = ",")
              sgRNAs_selected[x,"nFirst_sgRNA_targeted_transcripts"] <- nfirst_sgRNA_transcripts
              sgRNAs_selected[x,"second_sgRNA_targeted_transcripts"] <- paste(second_sgRNA_transcripts, collapse = ",")
              sgRNAs_selected[x,"nSecond_sgRNA_targeted_transcripts"] <- nsecond_sgRNA_transcripts
              sgRNAs_selected[x,"TranscriptsInFrame"] <- paste(TranscriptsInFrame, collapse = ",")
              sgRNAs_selected[x,"nTranscriptsInFrame"] <- nTranscriptsInFrame
              sgRNAs_selected[x,"TranscriptsOutOfFrame"] <- paste(TranscriptsOutOfFrame, collapse = ",")
              sgRNAs_selected[x,"nTranscriptsOutOfFrame"] <- nTranscriptsOutOfFrame
            } 
          }
        }
        
        sgRNAs_selected <- sgRNAs_selected %>%
          arrange(EntrezID, desc(proximity_1kb), desc(produces_frameshift), final_rank)
          
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
  }
  
  # if(length(sgRNAs_selected$maps_to_genome %>% unique)==1){
  #   sgRNAs_selected <- sgRNAs_selected %>%
  #     select(-maps_to_genome)
  # }
  
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
