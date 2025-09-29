geneOfInterestSidebarTab <- tabItem(
  tabName = "geneOfInterestSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#geneOfInterestInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "geneOfInterestInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "CRISPR/Cas9 screen hits",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput("geneOfInterestScreenHitsTable"))
      ),
      box(
        title = "Differental expression analysis hits",
        status = "info",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput("geneOfInterestDEAHitsTable"))
      ),
      box(
        title = "CRISPR/Cas9 screen meta data",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput("geneOfInterestScreenMetaTable"))
      ),
      box(
        title = "Differental expression analysis meta data",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput("geneOfInterestDEAMetaTable"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "geneOfInterestSpeciesSelect",
          "Species:",
          choices = list("Human" = "human", "Mouse" = "mouse"),
          selected = "human",
          inline = TRUE
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "geneOfInterestSearchRadio",
          "Screens filtering level:",
          choices = list("Genes" = "gene_stats", "Guides" = "guide_stats"),
          selected = "gene_stats",
          inline = TRUE
        ),
        checkboxGroupInput(
          "geneOfInterestScreenInclude",
          "Include screen gene hits:",
          choices = list(
            "Depleting" = "depletion",
            "Enriching" = "enrichment"
          ),
          selected = c("depletion", "enrichment"),
          inline = TRUE
        ),
        sliderInput(
          "geneOfInterestScreenDepletingLFC",
          "Minimal depleting LFC in CRISPR/Cas9 screen:",
          min = -8,
          max = 0,
          value = -2,
          step = 0.1
        ),
        sliderInput(
          "geneOfInterestScreenEnrichingLFC",
          "Minimal enriching LFC in CRISPR/Cas9 screen:",
          min = 0,
          max = 6,
          value = 2,
          step = 0.1
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "geneOfInterestDEAInclude",
          "Include DEA gene hits:",
          choices = list(
            "Depleting" = "depletion",
            "Enriching" = "enrichment"
          ),
          selected = c("depletion", "enrichment"),
          inline = TRUE
        ),
        sliderInput(
          "geneOfInterestDEADepletingLFC",
          "Minimal depleting LFC in differential expression analysis:",
          min = -3,
          max = 0,
          value = -1,
          step = 0.1
        ),
        sliderInput(
          "geneOfInterestDEAEnrichingLFC",
          "Minimal enriching LFC in differential expression analysis:",
          min = 0,
          max = 3,
          value = 1,
          step = 0.1
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        selectizeInput(
          "geneOfInterestGeneSelect",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        fileInput(
          "geneOfInterestGeneInputFile",
          "Upload list of genes:",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain"
          )
        ),
        disabled(actionButton("geneOfInterestLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "geneOfInterestScreenHitsButtonDownload",
          label = "Download screen hits table"
        ),
        downloadButton(
          width = NULL,
          outputId = "geneOfInterestDEAHitsButtonDownload",
          label = "Download DEA hits table"
        ),
        downloadButton(
          width = NULL,
          outputId = "geneOfInterestScreenMetaButtonDownload",
          label = "Download screen metadata table"
        ),
        downloadButton(
          width = NULL,
          outputId = "geneOfInterestDEAMetaButtonDownload",
          label = "Download DEA metadata table"
        )
      )
    )
  )
)

geneOfInterestSidebarTab
