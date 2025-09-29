dualSgRNAsPredictCombinationsSidebarTab <- tabItem(
  tabName = "dualSgRNAsPredictCombinationsSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#dualSgRNAsInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "dualSgRNAsInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("dualSgRNAsTableOutput"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxInput(
          "dualSgRNAs_LimitOutput",
          "Limit the number of reported dual-sgRNA-combinations per gene (check to activate)",
          value = FALSE
        ),
        disabled(
          sliderInput(
            "dualSgRNAs_nOutput",
            "Integer:",
            min = 0,
            max = 50,
            value = 25
          )
        ),
        fileInput(
          "dualSgRNAs_inputFile",
          "Choose CSV File",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        actionButton("dualSgRNALoadButton", "Load data!")
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "dualSgRNAsButtonDownload",
          label = "Download displayed table"
        )
      )
    )
  )
)

dualSgRNAsPredictCombinationsSidebarTab
