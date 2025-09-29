gwsGeneTab <- tabItem(
  tabName = "gwsGeneTab",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#gwsGeneInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "gwsGeneInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Screens",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput("gwsGeneTable"))
      ),
      box(
        title = "Contrasts",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "gwsGeneContrastTable"))
      ),
      box(
        title = "Samples",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "gwsGeneSampleTable"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "gwsGeneSearchRadio",
          "Search:",
          choices = list("Genes" = "gene_id", "Guides" = "guide_id"),
          selected = "gene_id",
          inline = TRUE
        ),
        radioButtons(
          "gwsGeneSpeciesSelect",
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
          "gwsGeneIndexRadio",
          "Display data as:",
          choices = list(
            "Log-fold change" = "lfc",
            "Effect" = "effect_essentialome",
            "FDR-adjusted effect" = "adjusted_effect_essentialome"
          ),
          selected = "adjusted_effect_essentialome",
          inline = FALSE
        ),
        radioButtons(
          "gwsGeneDisplayName",
          "Display screen IDs as:",
          choices = list(
            "Unique ID" = "long",
            "Short quality-control ID" = "short"
          ),
          selected = "long",
          inline = TRUE
        ),
        checkboxGroupInput(
          "gwsGeneInclude",
          "Include gene-level statistics:",
          choices = list(
            "P-value" = "p",
            "FDR" = "fdr",
            "Guides-good" = "guides_good",
            "Guides-total" = "guides"
          ),
          selected = NULL,
          inline = FALSE
        ),
        sliderInput(
          "gwsGeneQuality",
          "Minimal dynamic range of screens (quality filtering):",
          min = 0,
          max = 8,
          value = 1.5,
          step = 0.1
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "gwsGeneDatasetSelect",
          "Dataset:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsGeneCheckDatasetAll",
          "Check all dataset types",
          value = FALSE
        ),
        selectizeInput(
          "gwsGeneReferenceSelect",
          "Reference type:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsGeneCheckReferenceAll",
          "Check all reference types",
          value = FALSE
        ),
        selectizeInput(
          "gwsGeneTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsGeneCheckTissueAll",
          "Search All Tissues",
          value = FALSE
        ),
        selectizeInput(
          "gwsGeneCellLineSelect",
          "Cell line:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsGeneCheckCellLineAll",
          "Browse all cell lines",
          value = FALSE
        ),
        disabled(
          selectizeInput(
            "gwsGeneLibrarySelect",
            "Library:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "gwsGeneCheckLibraryAll",
            "Search All Libraries",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "gwsGeneContrastSelect",
            "Contrast:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "gwsGeneCheckContrastAll",
            "Search All Contrasts",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "gwsGeneGeneSelect",
            "Gene:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          fileInput(
            "gwsGeneGeneInputFile",
            "Upload list of genes:",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain"
            )
          )
        ),
        disabled(actionButton("gwsGeneLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "gwsGeneButtonDownload",
          label = "Download displayed table"
        )
      )
    )
  )
)

gwsGeneTab
