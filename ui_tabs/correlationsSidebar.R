correlationsSidebarTab <- tabItem(
  tabName = "correlationsSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#correlationsInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "correlationsInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Dependency <> Expression",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        withSpinner(
          dataTableOutput(outputId = "correlationsDependencyExpressionTableOutput")
        )
      ),
      box(
        title = "Expression <> Dependency",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        withSpinner(
          dataTableOutput(outputId = "correlationsExpressionDependencyTableOutput")
        )
      ),
      box(
        title = "Co-essentiality",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        withSpinner(
          dataTableOutput(outputId = "correlationsCoEssentialityTableOutput")
        )
      ),
      box(
        title = "Co-expression",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        withSpinner(
          dataTableOutput(outputId = "correlationsCoExpressionTableOutput")
        )
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "correlationsTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = FALSE,
          selected = NULL
        ),
        selectizeInput(
          "correlationsGeneSelect",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        fileInput(
          "correlationsGeneInputFile",
          "Upload list of genes:",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain"
          )
        ),
        sliderInput(
          "correlationsSliderCoeff",
          "Show top 20 correlations per gene or correlations with a pearson coeff above:",
          min = 0.3,
          max = 0.8,
          value = 0.6
        ),
        sliderInput(
          "correlationsSliderDatapoints",
          "Minimum amount of available datapoints for correlation:",
          min = 10,
          max = 500,
          value = 150,
          step = 1
        ),
        disabled(actionButton("correlationsLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "correlationsDownloadCheck",
          "Download result tables:",
          choices = list(
            "Dependency <> Expression" = "Dependency <> Expression",
            "Expression <> Dependency" = "Expression <> Dependency",
            "Co-Essentiality" = "Co-Essentiality",
            "Co-Expression" = "Co-Expression"
          ),
          selected = c(
            "Dependency <> Expression",
            "Expression <> Dependency",
            "Co-Expression",
            "Co-Essentiality"
          ),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "correlationsButtonDownload",
          label = "Download"
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "correlationsDownloadPrimaryTablesCheck",
          "Download primary data (top 20 correlation per gene or correlations with a pearson coeff > 0.6):",
          choices = list(
            "Dependency <> Expression" = "Dependency <> Expression",
            "Expression <> Dependency" = "Expression <> Dependency",
            "Co-Essentiality" = "Co-Essentiality",
            "Co-Expression" = "Co-Expression"
          ),
          selected = c(
            "Dependency <> Expression",
            "Expression <> Dependency",
            "Co-Expression",
            "Co-Essentiality"
          ),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "correlationsButtonDownloadPrimaryTables",
          label = "Download primary data"
        )
      )
    )
  )
)

correlationsSidebarTab
