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
  
  #Disable menuitem when the app loads
  if(view == "external"){
    addCssClass(selector = "a[data-value='sgRNAInfoSidebar']", class = "inactiveLink")
    addCssClass(selector = "a[data-value='sgRNAsSidebar']", class = "inactiveLink")
    addCssClass(selector = "a[data-value='dualSgRNAsSidebar']", class = "inactiveLink")
    addCssClass(selector = "a[data-value='expressionDataSidebar']", class = "inactiveLink")
    addCssClass(selector = "a[data-value='dualSgRNAsSidebar']", class = "inactiveLink")
  }

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

  # ----------------------------------------------------------------------------
  # ExpressionData
  # ----------------------------------------------------------------------------
  
  source(file = "server/expression_data.R", local = T)
  
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
  
}