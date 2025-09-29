dualSgRNAsTopCombinationsSidebarTab <- tabItem(
  tabName = "dualSgRNAsTopCombinationsSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#dualSgRNAsTopCombinationsInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "dualSgRNAsTopCombinationsInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("dualSgRNAsTopCombinationsTableOutput"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "dualSgRNAsTopCombinationsSpeciesSelect",
          "Species:",
          choices = list("Human" = "human", "Mouse" = "mouse"),
          selected = "human",
          inline = TRUE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        sliderInput(
          "dualSgRNAsTopCombinationsnOutput",
          "Limit the number of reported dual-sgRNA-combinations per gene to:",
          min = 0,
          max = 25,
          value = 5
        ),
        selectizeInput(
          "dualSgRNAsTopCombinationsGeneSelect",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        fileInput(
          "dualSgRNAsTopCombinationsGeneInputFile",
          "Upload list of genes:",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain"
          )
        ),
        disabled(actionButton("dualSgRNAsTopCombinationsLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "dualSgRNAsTopCombinationsButtonDownload",
          label = "Download displayed table"
        )
      )
    )
  )
)

dualSgRNAsTopCombinationsSidebarTab
