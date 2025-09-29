cellLineSidebarTab <- tabItem(
  tabName = "cellLineSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#cellLineInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "cellLineInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Cellline meta data",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineDataTableMeta")
        )
      ),
      box(
        title = "Gene mutations ",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineDataTableMutations")
        )
      ),
      box(
        title = "Gene fusions ",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineDataTableFusions")
        )
      ),
      box(
        title = "Gene CNVs",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineDataTableCNVs")
        )
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "cellLineTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "cellLineCheckTissueAll",
          "Search All Tissues",
          value = FALSE
        ),
        disabled(
          selectizeInput(
            "cellLineCellLineSelect",
            "Cell Line:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "cellLineCheckCellLineAll",
            "Search All Cell Lines",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "cellLineGeneSelect",
            "Gene:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          fileInput(
            "cellLine_inputFile",
            "Upload list of genes:",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain"
            )
          )
        ),
        disabled(
          checkboxInput(
            "cellLineCheckGeneAll",
            "Search All Genes",
            value = FALSE
          )
        ),
        disabled(actionButton("cellLineLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "cellLineDownloadCheck",
          "Download result tables:",
          choices = list(
            "Cellline meta data" = "Cellline meta data",
            "Gene mutations" = "Gene mutations",
            "Gene fusions" = "Gene fusions",
            "Gene CNVs" = "Gene CNVs"
          ),
          selected = c(
            "Cellline meta data",
            "Gene mutations",
            "Gene fusions",
            "Gene CNVs"
          ),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "cellLineButtonDownload",
          label = "Download"
        )
      )
    )
  )
)

cellLineSidebarTab
