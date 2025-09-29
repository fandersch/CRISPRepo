libSidebarTab <- tabItem(
  tabName = "libSidebar",
  width = NULL,
  fluidRow(
    column(
      width = 9,
      box(
        width = NULL,
        solidHeader = TRUE,
        withSpinner(dataTableOutput("libTableOutput"))
      )
    ),
    column(
      width = 3,
      box(
        width = NULL,
        solidHeader = TRUE,
        radioButtons(
          "libSpeciesSelect",
          "Species:",
          choices = list(
            "Human" = "human",
            "Mouse" = "mouse",
            "All" = "all"
          ),
          selected = "human",
          inline = TRUE
        ),
        selectInput(
          "libSelectLibrary",
          "Library:",
          choices = NULL,
          selectize = TRUE
        )
      ),
      infoBoxOutput(width = NULL, "libBoxGuidesTotal"),
      infoBoxOutput(width = NULL, "libBoxGenesTotal"),
      downloadButton(
        width = NULL,
        outputId = "libButtonDownload",
        label = "Download displayed table"
      )
    )
  )
)

libSidebarTab
