patientMutationSidebarTab <- tabItem(
  tabName = "patientMutationSidebar",
  width = NULL,
  fluidRow(
    tags$head(
      tags$style(
        HTML("#patientMutationInfo{color:tomato; font-weight: bold;}")
      )
    ),
    column(
      width = 12,
      box(
        width = NULL,
        solidHeader = TRUE,
        htmlOutput(outputId = "patientMutationInfo")
      )
    )
  ),
  fluidRow(
    column(
      width = 9,
      box(
        title = "Gene alterations heatmap",
        status = "success",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          imageOutput(
            outputId = "patientMutationHeatmapAlterations",
            width = "100%",
            height = "20px"
          )
        )
      ),
      box(
        title = "Gene alterations",
        status = "warning",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "patientMutationDataTableAlterations")
        )
      ),
      box(
        title = "Gene mutations ",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "patientMutationDataTableMutations")
        )
      ),
      box(
        title = "Gene fusions ",
        status = "danger",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        withSpinner(
          dataTableOutput(outputId = "patientMutationDataTableFusions")
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
          dataTableOutput(outputId = "patientMutationDataTableCNVs")
        )
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "patientMutationDataSourceSelect",
          "Data sources:",
          choices = list(
            "TCGA PanCancer Atlas Studies" = "tcga_pan_cancer",
            "All Non-redundant Studies" = "all",
            "Pediatric Cancer Studies" = "pediatric",
            "Non-pediatric Cancer Studies" = "non_pediatric"
          ),
          selected = "tcga_pan_cancer",
          inline = FALSE
        ),
        radioButtons(
          "patientMutationDataTypeSelect",
          "Focus on:",
          choices = list(
            "Alterations (any gene modification)" = "alterations",
            "Deletions" = "deletions",
            "Amplifications" = "amplifications",
            "Mutations" = "mutations",
            "Fusions" = "fusions"
          ),
          selected = "alterations",
          inline = FALSE
        ),
        radioButtons(
          "patientMutationDataGrouping",
          "Grouping by:",
          choices = list(
            "Cancer type" = "cancer_type",
            "Study" = "study_name"
          ),
          selected = "cancer_type",
          inline = TRUE
        ),
        selectizeInput(
          "patientMutationGroupSelect",
          "Group:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        checkboxInput(
          "patientMutationCheckGroupAll",
          "Search all groups",
          value = FALSE
        ),
        selectizeInput(
          "patientMutationGeneSelect",
          "Gene:",
          choices = NULL,
          multiple = TRUE,
          selected = NULL
        ),
        fileInput(
          "patientMutationGeneInputFile",
          "Upload list of genes:",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain"
          )
        ),
        checkboxInput(
          "patientMutationCheckAgeFiltering",
          "Filter patients by Age",
          value = FALSE
        ),
        disabled(
          sliderInput(
            "patientMutationAgeFiltering",
            "Age of patients:",
            min = 0,
            max = 100,
            value = c(0, 18)
          )
        ),
        sliderInput(
          "patientMutationMinPatients",
          "Minimum amount of patients in group:",
          min = 10,
          max = 500,
          value = 25
        ),
        sliderInput(
          "patientMutationMaxHeatmapColour",
          "Group percentage threshold for maximum color:",
          min = 2,
          max = 50,
          value = 5
        ),
        disabled(actionButton("patientMutationLoadButton", "Load data!"))
      ),
      box(
        width = NULL,
        solidHeader = TRUE,
        checkboxGroupInput(
          "patientMutationDownloadCheck",
          "Download result tables:",
          choices = list(
            "Gene alterations heatmap" = "Gene alterations heatmap",
            "Gene alterations" = "Gene alterations",
            "Gene mutations" = "Gene mutations",
            "Gene fusions" = "Gene fusions",
            "Gene CNVs" = "Gene CNVs"
          ),
          selected = c(
            "Gene alterations heatmap",
            "Gene alterations",
            "Gene mutations",
            "Gene fusions",
            "Gene CNVs"
          ),
          inline = FALSE
        ),
        downloadButton(
          width = NULL,
          outputId = "patientMutationButtonDownload",
          label = "Download"
        )
      )
    )
  )
)

patientMutationSidebarTab
