# crisprepo-shiny

# Copyright (c) 2018 Tobias Neumann, Jesse Lipp, Florian Andersch.
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
  useShinyjs(),
  sidebarMenu(id="tabs", 
              tags$head(tags$style(".inactiveLink { 
                                      visibility: hidden;
                                      }")),
              menuItem(text = "Genome-wide Screens", 
                       tabName = "gwsSidebar",startExpanded=TRUE,
                       menuSubItem(text = "Browse Screen", tabName = "gwsBrowseScreenTab", selected = TRUE), 
                       menuSubItem(text = "Gene Search", tabName = "gwsGeneTab")),
              menuItem(text = "Libraries", tabName = "libSidebar"),
              menuItemOutput("sgRNAInfoSidebar"),
              menuItemOutput("sgRNAsSidebar"),
              menuItemOutput("dualSgRNAsSidebar"),
              menuItemOutput("expressionDataSidebar"),
              menuItemOutput("essentialomeSidebar"),
              menuItemOutput("correlationsSidebar"),
              menuItemOutput("cellLineSidebar")
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
                     box(title = "Screens", status = "success", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput(outputId="gwsBrowseScreenTable"))),
                     box(title = "Contrasts", status = "warning", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, withSpinner(dataTableOutput(outputId="gwsBrowseScreenContrastTable"))),
                     box(title = "Samples", status = "danger", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, withSpinner(dataTableOutput(outputId="gwsBrowseScreenSampleTable")))),
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
                         
                         bsTooltip("gwsBrowseScreenIndexRadio", HTML("Effect: scaled LFCs between -1 (essential) and 0 (nonessential) <br/><br/> FDR-adjusted effect: scaled LFCs gets adjusted (weakened) for false positive (high FDR) genes. Scaled LFCs incoorpoprate its gene-level FDR and good-guides ratio (calculated by mageck) to reduce amount of false positives."), placement = "left", trigger = "hover"),

                         radioButtons(
                           "gwsBrowseScreenIndexRadio",
                           label = "Display data as:",
                           choices = list("Log-fold change" = "lfc", "Effect" = "effect", "FDR-adjusted effect" = "adjusted_effect"),
                           selected = "adjusted_effect",
                           inline = F
                         ),
                         
                         radioButtons(
                           "gwsBrowseScreenDisplayName",
                           label = "Display screen IDs as:",
                           choices = list("Short quality-control ID" = "short", "Unique ID" = "long"),
                           selected = "short",
                           inline = T
                         ),
                         
                         checkboxGroupInput(
                           "gwsBrowseScreenInclude",
                           label = "Include gene-level statistics:",
                           choices = list("P-value" = "p", "FDR" = "fdr", "Guides-good" = "guides_good", "Guides-total" = "guides"),
                           selected = NULL,
                           inline = F
                         )
                         ,
                         sliderInput("gwsBrowseScreenQuality", "Minimal dynamic range of screens (quality filtering):",
                                     min = 0, 
                                     max = 8, 
                                     value = 1.5,
                                     step = 0.1
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "gwsBrowseScreenDatasetSelect",
                           label = "Dataset:",
                           choices = dataset_selection_all,
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
                         selectizeInput(
                           inputId = "gwsBrowseScreenCellLineSelect",
                           label = "Cell line:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "gwsBrowseScreenCheckCellLineAll",
                           label = "Browse all cell lines",
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
                         label = "Download displayed table"
                       )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       checkboxGroupInput(
                         "gwsBrowseScreenDownloadPrimaryTablesCheck",
                         label = "Download HQ dropout screen data (scaled LFCs table):",
                         choices = list("Human" = "Human", "Mouse" = "Mouse"),
                         selected = c("Human", "Mouse"),
                         inline = F
                       ),
                       downloadButton(
                         width = NULL,
                         outputId = "gwsBrowseScreenButtonDownloadPrimaryTables",
                         label = "Download essentiality data"
                       )),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       checkboxGroupInput(
                         "gwsBrowseScreenDownloadAdjustedPrimaryTablesCheck",
                         label = "Download fdr-adjusted HQ dropout screen data (fdr-adjusted scaled LFCs table):",
                         choices = list("Human" = "Human", "Mouse" = "Mouse"),
                         selected = c("Human", "Mouse"),
                         inline = F
                       ),
                       downloadButton(
                         width = NULL,
                         outputId = "gwsBrowseScreenButtonDownloadAdjustedPrimaryTables",
                         label = "Download fdr-adjusted essentiality data"
                       ))
                     
              )
            )),
    
    
    #Gene Search
    tabItem(tabName = "gwsGeneTab", width = NULL,
            fluidRow(tags$head(tags$style(HTML('#gwsGeneInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput(outputId="gwsGeneInfo")))),
            fluidRow(
              column(width = 9,
                     box(title = "Screens", status = "success", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput("gwsGeneTable"))),
                     box(title = "Contrasts", status = "warning", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, withSpinner(dataTableOutput(outputId="gwsGeneContrastTable"))),
                     box(title = "Samples", status = "danger", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=TRUE, withSpinner(dataTableOutput(outputId="gwsGeneSampleTable")))),
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
                           label = "Display data as:",
                           choices = list("Log-fold change" = "lfc", "Effect" = "effect", "FDR-adjusted effect" = "adjusted_effect"),
                           selected = "adjusted_effect",
                           inline = F
                         ),
                         
                         radioButtons(
                           "gwsGeneDisplayName",
                           label = "Display screen IDs as:",
                           choices = list("Short quality-control ID" = "short", "Unique ID" = "long"),
                           selected = "short",
                           inline = T
                         ),
                         
                         checkboxGroupInput(
                           "gwsGeneInclude",
                           label = "Include gene-level statistics:",
                           choices = list("P-value" = "p", "FDR" = "fdr", "Guides-good" = "guides_good", "Guides-total" = "guides"),
                           selected = NULL,
                           inline = F
                         ),
                         sliderInput("gwsGeneQuality", "Minimal dynamic range of screens (quality filtering):",
                                     min = 0, 
                                     max = 8, 
                                     value = 1.5,
                                     step = 0.1
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "gwsGeneDatasetSelect",
                           label = "Dataset:",
                           choices = dataset_selection_all,
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
                         selectizeInput(
                           inputId = "gwsGeneCellLineSelect",
                           label = "Cell line:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "gwsGeneCheckCellLineAll",
                           label = "Browse all cell lines",
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
                           fileInput(
                             "gwsGeneGeneInputFile", "Upload list of genes:",
                             accept = c("text/csv","text/comma-separated-values,text/plain")
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
                         label = "Download displayed table"
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
                     box(width = NULL, solidHeader = TRUE,
                         downloadButton(width = NULL, 
                                        outputId = "sgRNAInfoButtonDownload",
                                        label = "Download displayed tables"
                        )
                     )
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
                                    label = "Download displayed table")
              )
            )
    ),
    
    # sgRNAs
    tabItem(tabName = "sgRNAsSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#sgRNAsInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput(outputId="sgRNAsInfo")))),
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
                           selected = "",
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
                           ), 
                           fileInput(
                             "sgRNAsGeneInputFile", "Upload list of genes:",
                             accept = c("text/csv","text/comma-separated-values,text/plain")
                           ),
                          disabled(
                             actionButton(inputId = "sgRNAsLoadButton", 
                                          label = "Load data!"
                             )
                           )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "sgRNAsButtonDownload",
                         label = "Download displayed table"
                       )
                     )
              )
            )),
    
    # dual sgRNAs predictions
    tabItem(tabName = "dualSgRNAsPredictCombinationsSidebar", width = NULL,
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
                         label = "Download displayed table"
                       )
                     )
              )
            )),
    
    # dual sgRNAs top predicted Combinations
    tabItem(tabName = "dualSgRNAsTopCombinationsSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#dualSgRNAsTopCombinationsInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput(outputId="dualSgRNAsTopCombinationsInfo")))),
            fluidRow(
              column(width = 9, 
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("dualSgRNAsTopCombinationsTableOutput")))
              ), 
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         
                         radioButtons(
                           "dualSgRNAsTopCombinationsSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "",
                           inline = T
                         )
                     ),
                     box(width = NULL, solidHeader = TRUE,
                         
                         selectizeInput(
                           inputId = "dualSgRNAsTopCombinationsGeneSelect",
                           label = "Gene:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ), 
                         fileInput(
                           "dualSgRNAsTopCombinationsGeneInputFile", "Upload list of genes:",
                           accept = c("text/csv","text/comma-separated-values,text/plain")
                         ),
                         disabled(
                           actionButton(inputId = "dualSgRNAsTopCombinationsLoadButton", 
                                        label = "Load data!"
                           )
                         )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       downloadButton(
                         width = NULL,
                         outputId = "dualSgRNAsTopCombinationsButtonDownload",
                         label = "Download displayed table"
                       )
                     )
              )
            )),
    
    # Expression Data
    tabItem(tabName = "expressionDataSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#expressionDataInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput((outputId="expressionDataInfo"))))),
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("expressionDataTable")))),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         
                         radioButtons(
                           "expressionDataSpeciesSelect",
                           label = "Species:",
                           choices = list(""),
                           selected = "",
                           inline = T
                         ),
                         radioButtons(
                           "expressionDataUnitSelect",
                           label = "Expression metric:",
                           choices = list(""),
                           selected = "",
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
                         disabled(fileInput("expressionData_inputFile", "Upload list of genes:",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain")
                         )),
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
                         label = "Download displayed table"
                       )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       checkboxGroupInput(
                         "expressionDataDownloadPrimaryTablesCheck",
                         label = "Download primary data (expression value table of all tissues):",
                         choices = list("Human" = "Human", "Mouse" = "Mouse"),
                         selected = c("Human", "Mouse"),
                         inline = F
                       ),
                       downloadButton(
                         width = NULL,
                         outputId = "expressionDataButtonDownloadPrimaryTables",
                         label = "Download primary data"
                       ))
              )
            ) 
    ),
    # Essentialome
    tabItem(tabName = "essentialomeSidebar", width = NULL, 
            fluidRow(
              column(width = 9, 
                     box(width = NULL, solidHeader = TRUE, withSpinner(dataTableOutput("essentialomeTableOutput")))
              ), 
              column(width = 3, 
                     box(width = NULL, solidHeader = TRUE, 
                         radioButtons(
                           "essentialomeSpeciesSelect",
                           label = "Species:",
                           choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"),
                           selected = "human",
                           inline = T
                         ),
                         checkboxGroupInput(
                           "essentialomeDisplay",
                           label = "Display:",
                           choices = list("Essential genes" = "essential", "Nonessential genes" = "nonessential"),
                           selected = c("essential", "nonessential"),
                           inline = T
                         ),
                         selectInput(inputId = "essentialomeSelectSource",
                                     label = "Source:",
                                     choices = NULL,
                                     selectize = TRUE)), 
                     infoBoxOutput(width = NULL, "essentialomeBoxEssentialGenesTotal"),
                     infoBoxOutput(width = NULL, "essentialomeBoxNonessentialGenesTotal"),
                     downloadButton(width = NULL, 
                                    outputId = "essentialomeButtonDownload",
                                    label = "Download")
              )
            )
    ),
    # correlations
    tabItem(tabName = "correlationsSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#correlationsInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput(outputId="correlationsInfo")))),
            fluidRow(
              column(width = 9, 
                     box(title = "Dependency <> Expression", status = "success", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput(outputId="correlationsDependencyExpressionTableOutput"))),
                     box(title = "Expression <> Dependency", status = "success", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput(outputId="correlationsExpressionDependencyTableOutput"))),
                     box(title = "Co-essentiality", status = "warning", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput(outputId="correlationsCoEssentialityTableOutput"))),
                     box(title = "Co-expression", status = "danger", width = NULL, solidHeader = TRUE, collapsible = TRUE, withSpinner(dataTableOutput(outputId="correlationsCoExpressionTableOutput")))
              ), 
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                       selectizeInput(
                         inputId = "correlationsTissueSelect",
                         label = "Tissue:",
                         choices = NULL,
                         multiple = FALSE,
                         selected = NULL
                       ),
                       selectizeInput(
                         inputId = "correlationsGeneSelect",
                         label = "Gene:",
                         choices = NULL,
                         multiple = TRUE,
                         selected = NULL
                       ),
                       fileInput(
                         "correlationsGeneInputFile", "Upload list of genes:",
                         accept = c("text/csv","text/comma-separated-values,text/plain")
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
                       disabled(
                         actionButton(inputId = "correlationsLoadButton", 
                                      label = "Load data!"
                         )
                       )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       checkboxGroupInput(
                         "correlationsDownloadCheck",
                         label = "Download result tables:",
                         choices = list("Dependency <> Expression" = "Dependency <> Expression", "Expression <> Dependency" = "Expression <> Dependency", "Co-Essentiality" = "Co-Essentiality", "Co-Expression" = "Co-Expression"),
                         selected = c("Dependency <> Expression", "Expression <> Dependency", "Co-Expression", "Co-Essentiality"),
                         inline = F
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
                         label = "Download primary data (top 20 correlation per gene or correlations with a pearson coeff > 0.6):",
                         choices = list("Dependency <> Expression" = "Dependency <> Expression","Expression <> Dependency" = "Expression <> Dependency", "Co-Essentiality" = "Co-Essentiality", "Co-Expression" = "Co-Expression"),
                         selected = c("Dependency <> Expression", "Expression <> Dependency", "Co-Expression", "Co-Essentiality"),
                         inline = F
                       ),
                       downloadButton(
                         width = NULL,
                         outputId = "correlationsButtonDownloadPrimaryTables",
                         label = "Download primary data"
                     ))
              )
          )
      ),
    # Cell line meta data 
    tabItem(tabName = "cellLineSidebar", width = NULL, 
            fluidRow(tags$head(tags$style(HTML('#cellLineInfo{color:tomato; font-weight: bold;}'))),
                     column(width = 12,
                            box(width = NULL, solidHeader = TRUE, htmlOutput((outputId="cellLineInfo"))))),
            fluidRow(
              column(width = 9,
                     box(title = "Cellline meta data", status = "info", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=F, withSpinner(dataTableOutput(outputId="cellLineDataTableMeta"))),
                     box(title = "Gene mutations ", status = "success", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=F, withSpinner(dataTableOutput(outputId="cellLineDataTableMutations"))),
                     box(title = "Gene fusions ", status = "warning", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=F, withSpinner(dataTableOutput(outputId="cellLineDataTableFusions"))),
                     box(title = "Gene CNVs", status = "danger", width = NULL, solidHeader = TRUE, collapsible = TRUE, collapsed=F, withSpinner(dataTableOutput(outputId="cellLineDataTableCNVs")))),
              column(width = 3,
                     box(width = NULL, solidHeader = TRUE,
                         selectizeInput(
                           inputId = "cellLineTissueSelect",
                           label = "Tissue:",
                           choices = NULL,
                           multiple = TRUE,
                           selected = NULL
                         ),
                         checkboxInput(
                           inputId = "cellLineCheckTissueAll",
                           label = "Search All Tissues",
                           value = FALSE
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "cellLineCellLineSelect",
                             label = "Cell Line:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(
                           checkboxInput(
                             inputId = "cellLineCheckCellLineAll",
                             label = "Search All Cell Lines",
                             value = FALSE
                           )
                         ),
                         disabled(
                           selectizeInput(
                             inputId = "cellLineGeneSelect",
                             label = "Gene:",
                             choices = NULL,
                             multiple = TRUE,
                             selected = NULL
                           )
                         ),
                         disabled(fileInput("cellLine_inputFile", "Upload list of genes:",
                                            accept = c(
                                              "text/csv",
                                              "text/comma-separated-values,text/plain")
                         )),
                         disabled(
                           checkboxInput(
                             inputId = "cellLineCheckGeneAll",
                             label = "Search All Genes",
                             value = FALSE
                           )
                         ),
                         disabled(
                           actionButton(inputId = "cellLineLoadButton", 
                                        label = "Load data!"
                           )
                         )
                     ),
                     box(
                       width = NULL,
                       solidHeader = TRUE,
                       checkboxGroupInput(
                         "cellLineDownloadCheck",
                         label = "Download result tables:",
                         choices = list("Cellline meta data" = "Cellline meta data", "Gene mutations" = "Gene mutations", "Gene fusions" = "Gene fusions", "Gene CNVs" = "Gene CNVs"),
                         selected = c("Cellline meta data", "Gene mutations", "Gene fusions", "Gene CNVs"),
                         inline = F
                       ),
                       downloadButton(
                         width = NULL,
                         outputId = "cellLineButtonDownload",
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