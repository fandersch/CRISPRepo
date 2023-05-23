# ----------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------

libTable <- reactive({
  con %>% 
    tbl("features") %>% 
    dplyr::filter(library_id %in% local(input$libSelectLibrary)) %>% 
    collect() %>%
    mutate(Length = nchar(sequence), Length_matching = nchar(sequence_matching), order_entrezID = as.numeric(entrez_id), 
           sgRNA_23_mer = ifelse(!gene_id %in% c("AMBIGUOUS", "UNMAPPED", "NOFEATURE", "SAFETARGETING", "NONTARGETING") & !is.na(context), paste0(gene_id, "_", substr(context, 5, nchar(context)-3)), NA)) %>%
    dplyr::arrange(ifelse(gene_id %in% c("AMBIGUOUS", "UNMAPPED", "NOFEATURE", "SAFETARGETING", "NONTARGETING"), 1, 0), order_entrezID, guide_id) %>%
    dplyr::select("Guide-ID" = guide_id, "Entrez-ID" = entrez_id, "Gene-Symbol" = symbol, "Ensembl-ID" = ensembl_id, Sequence = sequence, Length, Sequence_matching = sequence_matching, Length_matching,
           "sgRNA-ID-23-mer" = sgRNA_23_mer,
           "Library-ID" = library_id, Chromosome = chromosome, Strand = strand, "Genomic-Start-Position" = start, 
           "Genomic-End-Position" = end, "Perfect-matches-(PM)-total" = pm_total, "PM-with-PAM" = pm_pam, 
           "PM-with-PAM-in-CDS" = pm_pam_cds, "PM-with-PAM-in-CDS-unique" = pm_pam_cds_unique, 
           "30-mer-Genomic-Context" = context, "Legacy-ID" = legacy_id, "Legacy-Gene-Annotation" = legacy_group, Class = class) %>%
    discard(~all(is.na(.) | . ==""))
  
  
})

output$libTableOutput <- renderDataTable({
  libTable()
}, 
rownames= FALSE,
extensions = c('FixedColumns','FixedHeader'),
options = list(autoWidth = FALSE,
               headerCallback = JS(headerCallback),
               fixedColumns = list(leftColumns = 1),
               scrollX=TRUE,
               columnDefs = list(list(className = 'dt-center', targets = "_all")),
               pageLength = 25,
               lengthMenu = c(25, 50, 100, 200),
               searchHighlight = TRUE
               ),
               filter = list(position = 'top', clear = FALSE)
)

output$libBoxGuidesTotal <- renderInfoBox({
  infoBox(title = "Guides", value = libTable() %>% nrow())
})

output$libBoxGenesTotal <- renderInfoBox({
  infoBox(title = "Genes", value = libTable() %>% dplyr::select("Entrez-ID") %>% distinct() %>% nrow())
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

libLibraryList <- reactive({
  if(class(libraries)[1] == "tbl_SQLiteConnection"){
    libraries <<- libraries %>%
      collect()
  }
  
  if(input$libSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$libSpeciesSelect
  }
  
  libraries %>%
    dplyr::filter(species %in% speciesList) %>%
    dplyr::select(library_id) %>%
    distinct %>%
    arrange(library_id) %>%
    .$library_id
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------

observeEvent(input$libSpeciesSelect, {
  #update other species selects
  updateSpecies(input$libSpeciesSelect)
  #update library selectbox
  updateSelectizeInput(session, 'libSelectLibrary', choices = libLibraryList(), server = TRUE)
})

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$libButtonDownload <- downloadHandler(
  filename = function() {
    paste0(local(input$libSelectLibrary), ".txt")
  },
  content = function(file) {
    libTable() %>% write_tsv(file)
  }
)