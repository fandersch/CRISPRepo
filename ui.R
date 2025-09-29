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
  width = 250,
  useShinyjs(),
  sidebarMenu(id="tabs", 
              tags$head(tags$style(".inactiveLink { visibility: hidden; }")),
              menuItem(text = "Genome-wide Screens", 
                       tabName = "gwsSidebar", startExpanded=TRUE,
                       menuSubItem(text = "Browse Screen", tabName = "gwsBrowseScreenTab", selected = TRUE), 
                       menuSubItem(text = "Gene Search", tabName = "gwsGeneTab"),
                       menuSubItem(text = "Libraries", tabName = "libSidebar"),
                       menuItemOutput("essentialomeSidebar"),
                       menuItemOutput("sgRNAInfoSidebar")
              ),
              menuItem(text = "Expression data", 
                       tabName = "expressionSidebar", startExpanded=TRUE,
                       menuItemOutput("expressionDataSidebar"), 
                       menuItemOutput("slamseqDataSidebar")
              ),
              menuItem(text = "Genome-wide sgRNA predictions", 
                       tabName = "gwsSgRNASidebar", startExpanded=TRUE,
                       menuItemOutput("sgRNAsSidebar"),
                       menuItemOutput("dualSgRNAsTopCombinationsSidebar"),
                       menuItemOutput("dualSgRNAsPredictCombinationsSidebar")
              ),
              menuItem(text = "Other tools",
                       tabName = "tools", startExpanded = TRUE,
                       menuItemOutput("patientMutationSidebar"),
                       menuItemOutput("cellLineSidebar"),
                       menuItemOutput("cellLineSelectorSidebar"),
                       menuItemOutput("geneOfInterestSidebar"),
                       menuItemOutput("groupTestingSidebar"),
                       menuItemOutput("correlationsSidebar")
              )
  )
)

# --- Source each tab (Option A) ---
gwsBrowseScreenTab           <- source("ui_tabs/gwsBrowseScreenTab.R",           local = TRUE)$value
gwsGeneTab                   <- source("ui_tabs/gwsGeneTab.R",                   local = TRUE)$value
libSidebarTab                <- source("ui_tabs/libSidebar.R",                   local = TRUE)$value
essentialomeSidebarTab       <- source("ui_tabs/essentialomeSidebar.R",          local = TRUE)$value
sgRNAInfoSidebarTab          <- source("ui_tabs/sgRNAInfoSidebar.R",             local = TRUE)$value
expressionDataSidebarTab     <- source("ui_tabs/expressionDataSidebar.R",        local = TRUE)$value
slamseqDataSidebarTab        <- source("ui_tabs/slamseqDataSidebar.R",           local = TRUE)$value
sgRNAsSidebarTab             <- source("ui_tabs/sgRNAsSidebar.R",                local = TRUE)$value
dualSgRNAsTopCombinationsSidebarTab <- source("ui_tabs/dualSgRNAsTopCombinationsSidebar.R", local = TRUE)$value
dualSgRNAsPredictCombinationsSidebarTab <- source("ui_tabs/dualSgRNAsPredictCombinationsSidebar.R", local = TRUE)$value
patientMutationSidebarTab    <- source("ui_tabs/patientMutationSidebar.R",       local = TRUE)$value
cellLineSidebarTab           <- source("ui_tabs/cellLineSidebar.R",              local = TRUE)$value
cellLineSelectorSidebarTab   <- source("ui_tabs/cellLineSelectorSidebar.R",      local = TRUE)$value
geneOfInterestSidebarTab     <- source("ui_tabs/geneOfInterestSidebar.R",        local = TRUE)$value
groupTestingSidebarTab       <- source("ui_tabs/groupTestingSidebar.R",          local = TRUE)$value
correlationsSidebarTab       <- source("ui_tabs/correlationsSidebar.R",          local = TRUE)$value

body <- dashboardBody(
  mainPanel(tags$head(tags$script(jscode))),
  tabItems(
    gwsBrowseScreenTab,
    gwsGeneTab,
    libSidebarTab,
    essentialomeSidebarTab,
    sgRNAInfoSidebarTab,
    expressionDataSidebarTab,
    slamseqDataSidebarTab,
    sgRNAsSidebarTab,
    dualSgRNAsTopCombinationsSidebarTab,
    dualSgRNAsPredictCombinationsSidebarTab,
    patientMutationSidebarTab,
    cellLineSidebarTab,
    cellLineSelectorSidebarTab,
    geneOfInterestSidebarTab,
    groupTestingSidebarTab,
    correlationsSidebarTab
  )
)

ui <- tagList(
  tags$head(
    tags$style(HTML("
      html { font-size: 80% !important; }
      .sidebar-menu li a { font-size: 0.8em !important; }
      .box-body { font-size: 0.8em !important; }
      .box-title { font-size: 1em !important; }
      .item { font-size: 0.8em !important; }
      .shiny-input-checkbox { transform: scale(0.8); }
      .form-control { padding: 0.4em 0.8em !important; height: auto !important; }
    "))
  ),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('updateImageHeight', function(message) {
      $('#patientMutationHeatmapAlterations').css('height', message.height + 'px');
    });
  ")),
  dashboardPage(header, sidebar, body)
)