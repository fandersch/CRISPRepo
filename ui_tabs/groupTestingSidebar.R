groupTestingSidebarTab <- tabItem(
  tabName = "groupTestingSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#groupTestingInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "groupTestingInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Group testing plot",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          plotlyOutput(outputId = "groupTestingPlot", height = "auto")
        )
      ),
      box(
        title = "Group testing results",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "groupTestingDataTableResults")
        )
      ),
      box(
        title = "Group testing contrasts",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "groupTestingDataTableContrasts")
        )
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "groupTestingSpeciesSelect",
          "Species:",
          choices = list("Human" = "human"),
          selected = "human",
          inline = TRUE
        ),
        radioButtons(
          "groupTestingDatasetSelect",
          "Dataset:",
          choices = list(
            "Gene dependency (dropout screens)" = "dependency",
            "Gene expression" = "expression"
          ),
          selected = "dependency",
          inline = FALSE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "groupTestingSearchRadio",
          "Compare:",
          choices = list("Guide-level" = "guide", "Gene-level" = "gene"),
          selected = "gene",
          inline = TRUE
        ),
        selectizeInput(
          "groupTestingLibrarySelect",
          "sgRNA library:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "groupTestingCheckLibraryAll",
          " Consider all libraries",
          value = TRUE
        )
      ),
      box(
        title = "Target group",
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "groupTestingFilteringGroupSelect",
          "Tissue/Cancer type filtering groups:",
          choices = c(
            "Tissue" = "tissue_name",
            "Cancer Type (Cell-Model-Passports)" = "SangerCellModelPassports_cancer_type",
            "Cancer Type Detail (Cell-Model-Passports)" = "SangerCellModelPassports_cancer_type_detail",
            "Oncotree Subtype (DepMap)" = "BroadDepMap_OncotreeSubtype",
            "Oncrotree Primary Disease (DepMap)" = "BroadDepMap_OncotreePrimaryDisease"
          ),
          multiple = TRUE,
          selected = c("tissue_name")
        ),
        selectizeInput(
          "groupTestingFilteringValueSelect",
          "Tissue/Cancer type filtering values:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "groupTestingCheckTissueAll",
          "Consider all tissues/cancer types",
          value = TRUE
        ),
        selectizeInput(
          "groupTestingGeneDependencySelect",
          "Gene dependency:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        sliderInput(
          "groupTestingGeneDependencySlider",
          "Filter for cell lines with scaled gene dependency stronger or euqal than (-1 equals to mean essential gene dropout):",
          min = -1,
          max = 0,
          value = c(-0.66)
        ),
        selectizeInput(
          "groupTestingGeneExpressionSelect",
          "Gene expression:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        sliderInput(
          "groupTestingGeneExpressionSlider",
          "Filter for cell lines with gene expression higher or equal than (log2-TPMs):",
          min = 0,
          max = 15,
          value = c(5)
        ),
        selectizeInput(
          "groupTestingCellLineSelect",
          "Cell Line:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "groupTestingCheckCellLineAll",
          "Consider all cell lines",
          value = FALSE
        )
      ),
      disabled(
        box(
          id = "groupTestingRest",
          title = "Control group",
          width = NULL,
          solidHeader = TRUE,
          selectizeInput(
            "groupTestingRestFilteringGroupSelect",
            "Other tissue/cancer type filtering groups:",
            choices = c(
              "Tissue" = "tissue_name",
              "Cancer Type (Cell-Model-Passports)" = "SangerCellModelPassports_cancer_type",
              "Cancer Type Detail (Cell-Model-Passports)" = "SangerCellModelPassports_cancer_type_detail",
              "Oncotree Subtype (DepMap)" = "BroadDepMap_OncotreeSubtype",
              "Oncrotree Primary Disease (DepMap)" = "BroadDepMap_OncotreePrimaryDisease"
            ),
            multiple = TRUE,
            selected = c("tissue_name")
          ),
          selectizeInput(
            "groupTestingRestFilteringValueSelect",
            "Other tissue/cancer type filtering values:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          ),
          checkboxInput(
            "groupTestingRestCheckTissueAll",
            "Consider all other tissues/cancer types",
            value = FALSE
          ),
          selectizeInput(
            "groupTestingRestGeneDependencySelect",
            "Gene dependency:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          ),
          sliderInput(
            "groupTestingRestGeneDependencySlider",
            "Filter for cell lines with scaled gene dependency weaker or euqal than (-1 equals to mean essential gene dropout):",
            min = -1,
            max = 0,
            value = c(-0.33)
          ),
          selectizeInput(
            "groupTestingRestGeneExpressionSelect",
            "Gene expression:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          ),
          sliderInput(
            "groupTestingRestGeneExpressionSlider",
            "Filter for cell lines with gene expression lower or equal than (log2-TPMs):",
            min = 0,
            max = 15,
            value = c(0)
          ),
          selectizeInput(
            "groupTestingRestCellLineSelect",
            "Other Cell Line:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          ),
          checkboxInput(
            "groupTestingRestCheckCellLineAll",
            "Consider all other cell lines",
            value = FALSE
          )
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        disabled(actionButton("groupTestingLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "groupTestingButtonDownload",
          label = "Download"
        )
      )
    )
  )
)

groupTestingSidebarTab
