cellLineSelectorSidebarTab <- tabItem(
  tabName = "cellLineSelectorSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#cellLineSelectorInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "cellLineSelectorInfo")
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
          dataTableOutput(outputId = "cellLineSelectorDataTableMeta")
        )
      ),
      box(
        title = "Screen results",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineSelectorDataTableScreens")
        )
      ),
      box(
        title = "Expression levels",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineSelectorDataTableExpression")
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
          dataTableOutput(outputId = "cellLineSelectorDataTableMutations")
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
          dataTableOutput(outputId = "cellLineSelectorDataTableFusions")
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
          dataTableOutput(outputId = "cellLineSelectorDataTableCNVs")
        )
      ),
      box(
        title = "HLA types",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "cellLineSelectorDataTableHLAs")
        )
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "cellLineSelectorTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "cellLineSelectorCheckTissueAll",
          "Search All Tissues",
          value = FALSE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorGeneMutationSelect",
            "Gene mutation:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          selectizeInput(
            "cellLineSelectorGeneMutationProteinSelect",
            "Protein mutation:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorGeneDependencySelect",
            "Gene dependency:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        sliderInput(
          "cellLineSelectorGeneDependencySlider",
          "Filter for cell lines with scaled gene dependency (-1 equals to mean essential gene dropout):",
          min = -1,
          max = 0,
          value = -0.5
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorGeneExpressionSelect",
            "Gene expression:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        sliderInput(
          "cellLineSelectorGeneExpressionSlider",
          "Filter for cell lines with gene expression (log2-TPMs):",
          min = 0,
          max = 15,
          value = 1
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorHLAtypeSelect",
            "HLA-type:",
            choices = c(
              "-" = "none",
              "HLA-A" = "HLA-A",
              "HLA-B" = "HLA-B",
              "HLA-C" = "HLA-C",
              "HLA-DQA1" = "HLA-DQA1",
              "HLA-DQB1" = "HLA-DQB1",
              "HLA-DRB1" = "HLA-DRB1"
            ),
            multiple = FALSE,
            selected = NULL
          )
        ),
        disabled(
          selectizeInput(
            "cellLineSelectorHLAalleleSelect",
            "HLA-allele:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        sliderInput(
          "cellLineSelectorHLAExpressionSlider",
          "Filter for cell lines with HLA-type expression (RPKMs):",
          min = 0,
          max = 1000,
          value = 5
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorGeneFusion3primeSelect",
            "3' gene fusion :",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          selectizeInput(
            "cellLineSelectorGeneFusion5primeSelect",
            "5' gene fusion :",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(
          selectizeInput(
            "cellLineSelectorGeneCNVSelect",
            "Gene CNV:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          selectizeInput(
            "cellLineSelectorGeneCNVCategorySelect",
            "Gene CNV category:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(actionButton("cellLineSelectorLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "cellLineSelectorDownloadCheck",
          "Download result tables:",
          choices = list(
            "Cellline meta data" = "Cellline meta data",
            "Screen results" = "Screen results",
            "Gene mutations" = "Gene mutations",
            "Gene fusions" = "Gene fusions",
            "Gene CNVs" = "Gene CNVs",
            "HLA types" = "HLA types"
          ),
          selected = c(
            "Cellline meta data",
            "Screen results",
            "Gene mutations",
            "Gene fusions",
            "Gene CNVs",
            "HLA types" = "HLA types"
          ),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "cellLineSelectorButtonDownload",
          label = "Download"
        )
      )
    )
  )
)

cellLineSelectorSidebarTab
