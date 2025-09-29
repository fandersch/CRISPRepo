sgRNAsSidebarTab <- tabItem(
  tabName = "sgRNAsSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#sgRNAsInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "sgRNAsInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("sgRNAsTableOutput"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "sgRNAsSpeciesSelect",
          "Species:",
          choices = list("Human" = "human", "Mouse" = "mouse"),
          selected = "human",
          inline = TRUE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "sgRNAsGeneSelect",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        fileInput(
          "sgRNAsGeneInputFile",
          "Upload list of genes:",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain"
          )
        ),
        disabled(actionButton("sgRNAsLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "sgRNAsButtonDownload",
          label = "Download displayed table"
        )
      )
    )
  )
)

sgRNAsSidebarTab
