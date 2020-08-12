# crisprepo-shiny

# Copyright (c) 2018 Tobias Neumann, Jesse Lipp.
# 
# crisprepo-shiny is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# crisprepo-shiny is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


header <- dashboardHeader(
  title = "CRISPRepo"
)

sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs",
    menuItem(text = "Genome-wide Screens", 
             tabName = "gwsSidebar",startExpanded=TRUE,
             menuSubItem(text = "Browse Screen", tabName = "gwsBrowseScreenTab", selected = TRUE), 
             menuSubItem(text = "Gene Search", tabName = "gwsGeneTab")),
    menuItem(text = "SgRNA info", tabName = "sgRNAInfoSidebar"),
    menuItem(text = "Libraries", tabName = "libSidebar"),
    menuItem(text = "Genome-wide sgRNA predictions", tabName = "sgRNAsSidebar"),
    menuItem(text = "Dual sgRNA design", tabName = "dualSgRNAsSidebar"),
    menuItem(text = "Expression data", tabName = "expressionDataSidebar")
  )
)

body <- dashboardBody(

  tabItems(
    
    # Genome-wide Screens
  
    #Browse Screen
    tabItem(tabName = "gwsBrowseScreenTab", width = NULL,
            fluidRow(tags$head(tags$style(HTML('#gwsBrowseScreenInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, textOutput(outputId="gwsBrowseScreenInfo")))),
            fluidRow(
              useShinyjs(),
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput(outputId="gwsBrowseScreenTable")))),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         radioButtons(
                           "gwsBrowseScreenSearchRadio",
                           label = "Search:",
                           choices = list("Genes" = "gene_id", "Guides" = "guide_id"),
                           selected = "gene_id",
                           inline = T
                         ),

                         radioButtons(
                           "gwsBrowseScreenSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         ),

                         radioButtons(
                           "gwsBrowseScreenIndexRadio",
                           label = "Display:",
                           choices = list("Log-fold change" = "lfc", "Effect" = "effect"),
                           selected = "lfc",
                           inline = T
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "gwsBrowseScreenDatasetSelect",
                           label = "Dataset:",
                           choices = list("dropout" = "dropout", "drug_modifier" = "synthetic", "facs_based" = "facs"),
                           multiple = FALSE,
                           selected = "dropout"
                         ),
                         selectizeInput(
                           inputId = "gwsBrowseScreenTissueSelect",
                           label = "Tissue:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "gwsBrowseScreenCheckTissueAll",
                           label = "Browse all tissues",
                           value = FALSE
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "gwsBrowseScreenLibrarySelect",
                             label = "Library:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "gwsBrowseScreenCheckLibraryAll",
                             label = "Browse all libraries",
                             value = FALSE
                           )
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "gwsBrowseScreenContrastSelect",
                             label = "Screen:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "gwsBrowseScreenCheckContrastAll",
                             label = "Browse all Contrasts",
                             value = FALSE
                           )
                         ),
                         disabled(
                           actionButton(inputId = "gwsBrowseScreenLoadButton", 
                                        label = "Load data!"
                           )
                         )
                     ),

                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "gwsBrowseScreenButtonDownload",
                         label = "Download"
                       )
                     )
              )
            )),
    
    
    #Gene Search
    tabItem(tabName = "gwsGeneTab", width = NULL,
            fluidRow(tags$head(tags$style(HTML('#gwsGeneInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, textOutput(outputId="gwsGeneInfo")))),
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("gwsGeneTable")))),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         
                         radioButtons(
                           "gwsGeneSearchRadio",
                           label = "Search:",
                           choices = list("Genes" = "gene_id", "Guides" = "guide_id"),
                           selected = "gene_id",
                           inline = T
                         ),
                         
                         radioButtons(
                           "gwsGeneSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         ),
                         
                         radioButtons(
                           "gwsGeneIndexRadio",
                           label = "Display:",
                           choices = list("Log-fold change" = "lfc", "Effect" = "effect"),
                           selected = "lfc",
                           inline = T
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "gwsGeneDatasetSelect",
                           label = "Dataset:",
                           choices = list("dropout" = "dropout", "drug_modifier" = "synthetic", "facs_based" = "facs"),
                           multiple = FALSE,
                           selected = "dropout"
                         ),
                         selectizeInput(
                           inputId = "gwsGeneTissueSelect",
                           label = "Tissue:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "gwsGeneCheckTissueAll",
                           label = "Search All Tissues",
                           value = FALSE
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "gwsGeneLibrarySelect",
                             label = "Library:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "gwsGeneCheckLibraryAll",
                             label = "Search All Libraries",
                             value = FALSE
                             
                           )
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "gwsGeneContrastSelect",
                             label = "Contrast:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "gwsGeneCheckContrastAll",
                             label = "Search All Contrasts",
                             value = FALSE
                           )
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "gwsGeneGeneSelect",
                             label = "Gene:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           actionButton(inputId = "gwsGeneLoadButton", 
                                        label = "Load data!"
                           )
                         )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "gwsGeneButtonDownload",
                         label = "Download"
                       )
                     )
              )
            )), 
    
     
    # sgRNA info
    tabItem(tabName = "sgRNAInfoSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#sgRNAInfoInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, textOutput(outputId="sgRNAInfoInfo")))),
            fluidRow(
              column(width = 9, 
                     h4("SgRNA values in all used screens:"),
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("sgRNAInfoTableOutputScreens"))),
                     h4("SgRNA prediction scores:"),
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("sgRNAInfoTableOutputPredictions"))),
                     h4("SgRNA validation scores:"),
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("sgRNAInfoTableOutputValidations")))
              ), 
              column(width = 3, 
                     box(width = NULL, solidHeader = TRUE, 
                         radioButtons(
                           "sgRNAInfoSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         ),
                         radioButtons(
                           "sgRNAInfoIndexRadio",
                           label = "Display:",
                           choices = list("Log-fold change" = "lfc", "Effect" = "effect"),
                           selected = "lfc",
                           inline = T
                         ),
                         selectInput(inputId = "sgRNAInfoSelectGene",
                                     label = "Gene:",
                                     choices = NULL,
                                     multiple = TRUE,
                                     selected = NULL
                                     ),
                         selectInput(inputId = "sgRNAInfoSelectGuide",
                                     label = "Guide:",
                                     choices = NULL,
                                     multiple = TRUE,
                                     selected = NULL),
                         checkboxInput(
                           inputId = "sgRNAInfoCheckGuideAll",
                           label = "Select all guides",
                           value = FALSE
                         ),
                         disabled(
                           actionButton(inputId = "sgRNAInfoLoadButton", 
                                        label = "Load data!"
                           )
                         )),
                     downloadButton(width = NULL, 
                                    outputId = "sgRNAButtonDownload",
                                    label = "Download")
              )
            )
    ),
    
    # Library
    tabItem(tabName = "libSidebar", width = NULL, 
            fluidRow(
              column(width = 9, 
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("libTableOutput")))
              ), 
              column(width = 3, 
                     box(width = NULL, solidHeader = TRUE, 
                         radioButtons(
                           "libSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         ),
                         selectInput(inputId = "libSelectLibrary",
                                     label = "Library:",
                                     choices = NULL,
                                     selectize = TRUE)), 
                     infoBoxOutput(width = NULL, "libBoxGuidesTotal"),
                     infoBoxOutput(width = NULL, "libBoxGenesTotal"),
                     downloadButton(width = NULL, 
                                    outputId = "libButtonDownload",
                                    label = "Download")
              )
            )
    ),
    
    # sgRNAs
    tabItem(tabName = "sgRNAsSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#sgRNAsInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, textOutput(outputId="sgRNAsInfo")))),
            fluidRow(
              column(width = 9, 
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("sgRNAsTableOutput")))
              ), 
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         
                         radioButtons(
                           "sgRNAsSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         
                           selectizeInput(
                             inputId = "sgRNAsGeneSelect",
                             label = "Gene:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "sgRNAsButtonDownload",
                         label = "Download"
                       )
                     )
              )
            )),
    
    # dual sgRNAs
    tabItem(tabName = "dualSgRNAsSidebar", width = NULL,
            fluidRow(tags$head(tags$style(HTML('#dualSgRNAsInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput(outputId="dualSgRNAsInfo")))),
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("dualSgRNAsTableOutput")))
              ),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         checkboxInput(
                           inputId = "dualSgRNAs_LimitOutput",
                           label = "Limit the number of reported dual-sgRNA-combinations per gene (check to activate)",
                           value = FALSE
                         ),
                         disabled(
                          sliderInput("dualSgRNAs_nOutput", "Integer:",
                                     min = 0, max = 50,
                                     value = 25)
                         ),
                         fileInput("dualSgRNAs_inputFile", "Choose CSV File",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                         ),
                         
                           actionButton(inputId = "dualSgRNALoadButton", 
                                        label = "Load data!"
                           )
                         
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "dualSgRNAsButtonDownload",
                         label = "Download"
                       )
                     )
              )
            )),
    
    # Expression Data
    tabItem(tabName = "expressionDataSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#expressionDataInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, textOutput(outputId="expressionDataInfo")))),
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("expressionDataTable")))),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         
                         radioButtons(
                           "expressionDataSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All" = "all"),
                           selected = "human",
                           inline = T
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "expressionDataTissueSelect",
                           label = "Tissue:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "expressionDataCheckTissueAll",
                           label = "Search All Tissues",
                           value = FALSE
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "expressionDataCellLineSelect",
                             label = "Cell Line:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "expressionDataCheckCellLineAll",
                             label = "Search All Cell Lines",
                             value = FALSE
                           )
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "expressionDataGeneSelect",
                             label = "Gene:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "expressionDataCheckGeneAll",
                             label = "Search All Genes",
                             value = FALSE
                           )
                         ),
                         disabled(
                           actionButton(inputId = "expressionDataLoadButton", 
                                        label = "Load data!"
                           )
                         )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "expressionDataButtonDownload",
                         label = "Download"
                       )
                     )
              )
            ) 
    )
  )
)
  
dashboardPage(
  header,
  sidebar, 
  body
)