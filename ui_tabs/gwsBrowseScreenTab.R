gwsBrowseScreenTab <- tabItem(
  tabName = "gwsBrowseScreenTab",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML(
          "#gwsBrowseScreenInfo{
                                 color:tomato;
                                 font-weight:
                                 bold;
                                }
                                [data-title] {
                                  position: relative;
                                }
                                [data-title]:after {
                                  content: attr(data-title);
                                  background-color: #fed8b1;
                                  color: black;
                                  font-size: 100%;
                                  padding: 1px 5px 2px 5px;
                                  bottom: -2.5em;
                                  border: 4px solid rgb(255, 255, 255);
                                  border-radius: 5px;
                                  left: -100%;
                                  position: absolute;
                                  z-index: 99999;
                                  visibility: hidden;
                                  width: min-content;
                                }
                                [data-title]:hover:after {
                                  opacity: 1;
                                  transition: all 0.1s ease 0.5s;
                                  visibility: visible;
                                }"
        )
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        textOutput(outputId = "gwsBrowseScreenInfo")
      )
    )
  ),
  fluidRow(
    useShinyjs(),
    column(
      width = 9,
      box(
        title = "Screens",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "gwsBrowseScreenTable"))
      ),
      box(
        title = "Contrasts",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "gwsBrowseScreenContrastTable"))
      ),
      box(
        title = "Samples",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(dataTableOutput(outputId = "gwsBrowseScreenSampleTable"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "gwsBrowseScreenSearchRadio",
          "Search:",
          choices = list("Genes" = "gene_id", "Guides" = "guide_id"),
          selected = "gene_id",
          inline = TRUE
        ),
        radioButtons(
          "gwsBrowseScreenSpeciesSelect",
          "Species:",
          choices = list(
            "Human" = "human",
            "Mouse" = "mouse",
            "All" = "all"
          ),
          selected = "human",
          inline = TRUE
        ),
        bsTooltip(
          "gwsBrowseScreenIndexRadio",
          HTML(
            "Effect: scaled LFCs between -1 (essential) and 0 (nonessential) <br/><br/> FDR-adjusted effect: scaled LFCs gets adjusted (weakened) for false positive (high FDR) genes. Scaled LFCs incoorpoprate its gene-level FDR and good-guides ratio (calculated by mageck) to reduce amount of false positives."
          ),
          placement = "left",
          trigger = "hover"
        ),
        radioButtons(
          "gwsBrowseScreenIndexRadio",
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
          "gwsBrowseScreenDisplayName",
          "Display screen IDs as:",
          choices = list(
            "Unique ID" = "long",
            "Short quality-control ID" = "short"
          ),
          selected = "long",
          inline = TRUE
        ),
        checkboxGroupInput(
          "gwsBrowseScreenInclude",
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
          "gwsBrowseScreenQuality",
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
          "gwsBrowseScreenDatasetSelect",
          "Dataset:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsBrowseScreenCheckDatasetAll",
          "Check all dataset types",
          value = FALSE
        ),
        selectizeInput(
          "gwsBrowseScreenReferenceSelect",
          "Reference type:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsBrowseScreenCheckReferenceAll",
          "Check all reference types",
          value = FALSE
        ),
        selectizeInput(
          "gwsBrowseScreenTissueSelect",
          "Tissue:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsBrowseScreenCheckTissueAll",
          "Browse all tissues",
          value = FALSE
        ),
        selectizeInput(
          "gwsBrowseScreenCellLineSelect",
          "Cell line:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "gwsBrowseScreenCheckCellLineAll",
          "Browse all cell lines",
          value = FALSE
        ),
        disabled(
          selectizeInput(
            "gwsBrowseScreenLibrarySelect",
            "Library:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "gwsBrowseScreenCheckLibraryAll",
            "Browse all libraries",
            value = FALSE
          )
        ),
        disabled(
          selectizeInput(
            "gwsBrowseScreenContrastSelect",
            "Screen:",
            choices = NULL,
            multiple = TRUE,
            selected = NULL
          )
        ),
        disabled(
          checkboxInput(
            "gwsBrowseScreenCheckContrastAll",
            "Browse all Contrasts",
            value = FALSE
          )
        ),
        disabled(actionButton("gwsBrowseScreenLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        downloadButton(
          width = NULL,
          outputId = "gwsBrowseScreenButtonDownload",
          label = "Download displayed table"
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "gwsBrowseScreenDownloadPrimaryTablesCheck",
          "Download HQ dropout screen data (scaled LFCs table):",
          choices = list("Human" = "Human", "Mouse" = "Mouse"),
          selected = c("Human", "Mouse"),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "gwsBrowseScreenButtonDownloadPrimaryTables",
          label = "Download essentiality data"
        )
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "gwsBrowseScreenDownloadAdjustedPrimaryTablesCheck",
          "Download fdr-adjusted HQ dropout screen data (fdr-adjusted scaled LFCs table):",
          choices = list("Human" = "Human", "Mouse" = "Mouse"),
          selected = c("Human", "Mouse"),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "gwsBrowseScreenButtonDownloadAdjustedPrimaryTables",
          label = "Download fdr-adjusted essentiality data"
        )
      )
    )
  )
)

gwsBrowseScreenTab
