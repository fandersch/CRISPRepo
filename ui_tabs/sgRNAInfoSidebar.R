sgRNAInfoSidebarTab <- tabItem(
  tabName = "sgRNAInfoSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#sgRNAInfoInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        textOutput(outputId = "sgRNAInfoInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      h4("SgRNA values in all used screens:"),
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("sgRNAInfoTableOutputScreens"))
      ),
      h4("SgRNA prediction scores:"),
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("sgRNAInfoTableOutputPredictions"))
      ),
      h4("SgRNA validation scores:"),
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("sgRNAInfoTableOutputValidations"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "sgRNAInfoSpeciesSelect",
          "Species:",
          choices = list(
            "Human" = "human",
            "Mouse" = "mouse",
            "All" = "all"
          ),
          selected = "human",
          inline = TRUE
        ),
        radioButtons(
          "sgRNAInfoIndexRadio",
          "Display:",
          choices = list(
            "Log-fold change" = "lfc",
            "Effect" = "effect_essentialome"
          ),
          selected = "lfc",
          inline = TRUE
        ),
        selectInput(
          "sgRNAInfoSelectGene",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        selectInput(
          "sgRNAInfoSelectGuide",
          "Guide:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "sgRNAInfoCheckGuideAll",
          "Select all guides",
          value = FALSE
        ),
        disabled(actionButton("sgRNAInfoLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "sgRNAInfoButtonDownload",
          label = "Download displayed tables"
        )
      )
    )
  )
)

sgRNAInfoSidebarTab
