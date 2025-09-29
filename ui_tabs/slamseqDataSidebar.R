slamseqDataSidebarTab <- tabItem(
  tabName = "slamseqDataSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#slamseqDataInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "slamseqDataInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Expression Levels",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "slamseqTotalDataTable"))
      ),
      box(
        title = "Expression Dyamics",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "slamseqDynamicsDataTable"))
      ),
      box(
        title = "Sample Meta Data",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "slamseqSampleMetaTable"))
      ),
      box(
        title = "DEA results",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "slamseqDeseq2Table"))
      ),
      box(
        title = "DEA Meta Data",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "slamseqDeseq2MetaTable"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "slamseqDataIncludeDataset",
          "Include datasets:",
          choices = list("Quant-seq" = "quant", "SLAM-seq" = "slam"),
          selected = c("quant", "slam"),
          inline = FALSE
        ),
        radioButtons(
          "slamseqDataSpeciesSelect",
          "Species:",
          choices = list("Human" = "human", "Mouse" = "mouse"),
          selected = "human",
          inline = TRUE
        ),
        radioButtons(
          "slamseqDataUnitSelect",
          "Expression metric:",
          choices = list(""),
          selected = "",
          inline = TRUE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxInput(
          "slamseqDataCheckOnlyIncludeBaseline",
          "Only include baseline reference samples",
          value = FALSE
        ),
        selectizeInput(
          "slamseqDataTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        disabled(
          checkboxInput(
            "slamseqDataCheckTissueAll",
            "Search All Tissues",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "slamseqDataCellLineSelect",
            "Cell Line:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "slamseqDataCheckCellLineAll",
            "Search All Cell Lines",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "slamseqDataGeneSelect",
            "Gene:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          fileInput(
            "slamseqData_inputFile",
            "Upload list of genes:",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain"
            )
          )
        ),
        disabled(
          checkboxInput(
            "slamseqDataCheckGeneAll",
            "Search All Genes",
            value = FALSE
          )
        ),
        disabled(actionButton("slamseqDataLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "slamseqExpressionDataButtonDownload",
          label = "Download Expression Level table"
        ),
        downloadButton(
          width = NULL,
          outputId = "slamseqExpressionDynamicsDataButtonDownload",
          label = "Download Expression Dynamics table"
        ),
        tags$br(),
        downloadButton(
          width = NULL,
          outputId = "slamseqSampleMetaDataButtonDownload",
          label = "Download Sample metadata table"
        ),
        tags$br(),
        downloadButton(
          width = NULL,
          outputId = "slamseqDeseq2DataButtonDownload",
          label = "Download DEA results table"
        ),
        tags$br(),
        downloadButton(
          width = NULL,
          outputId = "slamseqDeseq2MetaDataButtonDownload",
          label = "Download DEA metadata table"
        )
      )
    )
  )
)

slamseqDataSidebarTab
