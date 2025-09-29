expressionDataSidebarTab <- tabItem(
  tabName = "expressionDataSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#expressionDataInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "expressionDataInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("expressionDataTable"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "expressionDataSpeciesSelect",
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
          "expressionDataUnitSelect",
          "Expression metric:",
          choices = list(""),
          selected = "",
          inline = FALSE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "expressionDataTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "expressionDataCheckTissueAll",
          "Search All Tissues",
          value = FALSE
        ),
        disabled(
          selectizeInput(
            "expressionDataCellLineSelect",
            "Cell Line:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "expressionDataCheckCellLineAll",
            "Search All Cell Lines",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "expressionDataGeneSelect",
            "Gene:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          fileInput(
            "expressionData_inputFile",
            "Upload list of genes:",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain"
            )
          )
        ),
        disabled(
          checkboxInput(
            "expressionDataCheckGeneAll",
            "Search All Genes",
            value = FALSE
          )
        ),
        disabled(actionButton("expressionDataLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "expressionDataButtonDownload",
          label = "Download displayed table"
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "expressionDataDownloadPrimaryTablesCheck",
          "Download primary data (expression value table of all tissues):",
          choices = list(
            "Human" = "Human",
            "Mouse" = "Mouse"
          ),
          selected = c("Human", "Mouse"),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "expressionDataButtonDownloadPrimaryTables",
          label = "Download primary data"
        )
      )
    )
  )
)

expressionDataSidebarTab
