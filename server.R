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

function(input, output, session) {
  
  output$sgRNAInfoSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "SgRNA info", tabName = "sgRNAInfoSidebar")
  })
  
  output$sgRNAsSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Genome-wide sgRNA predictions", tabName = "sgRNAsSidebar")
  })
  
  output$dualSgRNAsSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Dual sgRNA design", tabName = "dualSgRNAsSidebar", startExpanded=TRUE,
             menuSubItem(text = "Predict combinations", tabName = "dualSgRNAsPredictCombinationsSidebar", selected = TRUE), 
             menuSubItem(text = "Top combinations per gene", tabName = "dualSgRNAsTopCombinationsSidebar"))
  })
  
  output$expressionDataSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Expression data", tabName = "expressionDataSidebar")
  })
  
  output$slamseqDataSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "SLAM-seq data", tabName = "slamseqDataSidebar")
  })
  
  output$essentialomeSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Essentialome", tabName = "essentialomeSidebar")
  })
  
  output$correlationsSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Correlations", tabName = "correlationsSidebar")
  })
  
  output$cellLineSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Cellline meta data", tabName = "cellLineSidebar")
  })
  
  output$cellLineSelectorSidebar <- renderMenu({
    if(view == "internal")
      menuItem(text = "Cellline selector", tabName = "cellLineSelectorSidebar")
  })
  
  # ----------------------------------------------------------------------------
  # Browse Screen
  # ----------------------------------------------------------------------------
  
  source(file = "server/browse_screen.R", local = T)
  
  # ----------------------------------------------------------------------------
  # Gene Search
  # ----------------------------------------------------------------------------

  source(file = "server/gene_search.R", local = T)

  # ----------------------------------------------------------------------------
  # sgRNA Info
  # ----------------------------------------------------------------------------

  source(file = "server/sgRNA_info.R", local = T)

  # ----------------------------------------------------------------------------
  # Libraries
  # ----------------------------------------------------------------------------

  source(file = "server/libraries.R", local = T)

  # ----------------------------------------------------------------------------
  # genome-wide sgRNA predictions
  # ----------------------------------------------------------------------------

  source(file = "server/genomewide_sgRNA_predictions.R", local = T)

  # ----------------------------------------------------------------------------
  # dual sgRNA designs
  # ----------------------------------------------------------------------------

  source(file = "server/dual_sgRNA_design.R", local = T)
  source(file = "server/dual_sgRNA_design_top_pairs.R", local = T)

  # ----------------------------------------------------------------------------
  # ExpressionData
  # ----------------------------------------------------------------------------

  source(file = "server/expression_data.R", local = T)
  
  # ----------------------------------------------------------------------------
  # SlamseqData
  # ----------------------------------------------------------------------------
  
  source(file = "server/slamseq_data.R", local = T)

  # ----------------------------------------------------------------------------
  # Essentialome
  # ----------------------------------------------------------------------------

  source(file = "server/essentialome.R", local = T)

  # ----------------------------------------------------------------------------
  # Correlations
  # ----------------------------------------------------------------------------

  source(file = "server/correlations.R", local = T)
  updateSelectizeInput(session, 'correlationsGeneSelect', choices = gene_list_correlations, server = TRUE)
  updateSelectizeInput(session, 'correlationsTissueSelect', choices = tissue_list_correlations_tissue, server = TRUE)

  # ----------------------------------------------------------------------------
  # Cell line meta  data
  # ----------------------------------------------------------------------------

  source(file = "server/cell_line_meta_data.R", local = T)
  updateSelectizeInput(session, 'cellLineTissueSelect', choices = tissue_list_cellLine, server = TRUE)

  # ----------------------------------------------------------------------------
  # Cell line meta  data
  # ----------------------------------------------------------------------------

  source(file = "server/cell_line_selector.R", local = T)
  updateSelectizeInput(session, 'cellLineSelectorTissueSelect', choices = tissue_list_cellLine, server = TRUE)
  updateSelectizeInput(session, 'cellLineSelectorGeneMutationSelect', choices = gene_list_cellLine$gene_symbol, server = TRUE)


  # ----------------------------------------------------------------------------
  # Header callback
  # ----------------------------------------------------------------------------

  #rotate vertical
  headerCallback <- c(
    "function(thead, data, start, end, display){",
    "  var $ths = $(thead).find('th');",
    "  $ths.css({'vertical-align': 'bottom', 'white-space': 'nowrap'});",
    "  var betterCells = [];",
    "  $ths.each(function(){",
    "    var cell = $(this);",
    "    var newDiv = $('<div>', {height: 'auto', width: 'auto'});",
    "    var newInnerDiv = $('<div>', {text: cell.text()});",
    "    newDiv.css({margin: 'auto'});",
    "    newInnerDiv.css({",
    "      transform: 'rotate(180deg)',",
    "      'writing-mode': 'tb-rl',",
    "      'white-space': 'nowrap'",
    "    });",
    "    newDiv.append(newInnerDiv);",
    "    betterCells.push(newDiv);",
    "  });",
    "  $ths.each(function(i){",
    "    $(this).html(betterCells[i]);",
    "  });",
    "}"
  )

  # ----------------------------------------------------------------------------
  # update species select boxes
  # ----------------------------------------------------------------------------

  updateSpecies <- function(species){
    updateSelectizeInput(session, 'gwsBrowseScreenSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'gwsGeneSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'libSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'sgRNAInfoSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'sgRNAsSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'dualSgRNAsTopCombinationsSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'expressionDataSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'slamseqDataSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
    updateSelectizeInput(session, 'essentialomeSpeciesSelect', choices = list("Human" = "human", "Mouse" = "mouse", "All"="all"), selected = species, server = TRUE)
  }

}