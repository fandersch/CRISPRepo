# ----------------------------------------------------------------------------
# Essentialome
# ----------------------------------------------------------------------------

essentialomeTable <- reactive({
  if(input$essentialomeSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$essentialomeSpeciesSelect
  }
  
  essentialome %>% 
    filter(source %in% local(input$essentialomeSelectSource), species %in% speciesList, class %in% input$essentialomeDisplay) %>% 
    select(-source, -species)
})

output$essentialomeTableOutput <- renderDataTable({
  essentialomeTable()
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

output$essentialomeEssentialGenesTotal <- renderInfoBox({
  infoBox(title = "Number of essential genes", value = essentialomeTable() %>% filter(class == "essential") %>% nrow())
})

output$essentialomeNonessentialGenesTotal <- renderInfoBox({
  infoBox(title = "Number of nonessential genes", value = essentialomeTable() %>% filter(class == "nonessential") %>% nrow())
})

# ----------------------------------------------------------------------------
# Selectbox Lists
# ----------------------------------------------------------------------------

essentialomeSourceList <- reactive({
  
  if(input$essentialomeSpeciesSelect == "all"){
    speciesList <- c("human", "mouse")
  }else{
    speciesList <- input$essentialomeSpeciesSelect
  }
  
  essentialome %>% 
    filter(species %in% speciesList) %>% 
    select(source) %>%
    distinct %>%
    .$source
})

# ----------------------------------------------------------------------------
# Observers
# ----------------------------------------------------------------------------

observeEvent(input$essentialomeSpeciesSelect, {
  #update library selectbox
  updateSelectizeInput(session, 'essentialomeSelectSource', choices = essentialomeSourceList(), server = TRUE)
})

# ----------------------------------------------------------------------------
# Download handler
# ----------------------------------------------------------------------------

output$essentialomeButtonDownload <- downloadHandler(
  filename = function() {
    paste0(local(input$essentialomeSelectSource), "essential_nonessential_genes.txt")
  },
  content = function(file) {
    essentialomeTable() %>% write_tsv(file)
  }
)