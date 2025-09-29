essentialomeSidebarTab <- tabItem(
  tabName = "essentialomeSidebar",
  width = NULL,
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("essentialomeTableOutput"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "essentialomeSpeciesSelect",
          "Species:",
          choices = list(
            "Human" = "human",
            "Mouse" = "mouse",
            "All" = "all"
          ),
          selected = "human",
          inline = TRUE
        ),
        checkboxGroupInput(
          "essentialomeDisplay",
          "Display:",
          choices = list(
            "Essential genes" = "essential",
            "Nonessential genes" = "nonessential"
          ),
          selected = c("essential", "nonessential"),
          inline = TRUE
        ),
        selectInput(
          "essentialomeSelectSource",
          "Source:",
          choices = NULL,
          multiple = TRUE,
          selectize = TRUE
        ),
        checkboxInput(
          "essentialomeCheckSourceAll",
          "Display all essentialomes",
          value = TRUE
        )
      ),
      infoBoxOutput(width = NULL, "essentialomeBoxEssentialGenesTotal"),
      infoBoxOutput(width = NULL, "essentialomeBoxNonessentialGenesTotal"),
      downloadButton(
        width = NULL,
        outputId = "essentialomeButtonDownload",
        label = "Download"
      )
    )
  )
)

essentialomeSidebarTab
