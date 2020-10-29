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

library(shinydashboard)
library(tidyverse)
library(stringr)
library(DT)
library(rlang)
library(viridis)
library(forcats)
library(shinycssloaders)
library(shinyjs)
library(readxl)


con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "screen.db")

# con_facs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "screen_facs.db")

con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "expression_data.db")

con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "sgRNAs.db")


pheno <- con %>%
  tbl("pheno")

libraries <- pheno %>%
  select(library_id, cellline_name, tissue_name, species, type) %>%
  distinct 

species <- pheno %>%
  select(species) %>%
  distinct %>%
  .$species

features <- con %>%
  tbl("features") %>%
  select(guide_id, gene_id, symbol=hgnc_symbol, entrez_id, sequence, guide_id, context, library_id) %>%
  filter(gene_id != "AMBIGUOUS") %>%
  filter(gene_id != "UNMAPPED") %>%
  filter(gene_id != "NOFEATURE") %>%
  filter(gene_id != "SAFETARGETING") %>%
  filter(gene_id != "NONTARGETING") %>%
  filter(!is.na(gene_id)) %>%
  distinct

contrasts <- con %>%
  tbl("contrasts") %>%
  select(contrast_id, contrast_id_QC, library_id, cellline_name, tissue_name, species, type, dynamic_range)

gene_list_screens <- con %>%
  tbl("features") %>%
  filter(gene_id != "AMBIGUOUS") %>%
  filter(gene_id != "UNMAPPED") %>%
  filter(gene_id != "NOFEATURE") %>%
  filter(gene_id != "SAFETARGETING") %>%
  filter(gene_id != "NONTARGETING") %>%
  filter(!is.na(gene_id)) %>%
  select(symbol=hgnc_symbol, entrez_id, library_id) %>%
  distinct %>%
  arrange(symbol)

gene_list_human <- NULL
gene_list_mouse <- NULL 

if("hs_gw_zuber_v2" %in% (libraries %>% collect %>% .$library_id)){
  gene_list_human <- con_sgRNAs %>%
    tbl("genes_human") %>%
    select(Symbol, EntrezID) %>%
    distinct %>%
    arrange(Symbol) %>%
    collect
  
  gene_list_mouse <- con_sgRNAs %>%
    tbl("genes_mouse") %>%
    select(Symbol, EntrezID) %>%
    distinct %>%
    arrange(Symbol) %>%
    collect
}

cellline_list_expressionData <- con_expression %>%
  tbl("expression_data_meta_info") %>%
  select(sample_id, cell_line_name, tissue_name, species, unit) %>%
  distinct() %>%
  arrange(cell_line_name)

gene_list_expressionData <- con_expression %>%
  tbl("expression_data_genes")

loadExpressionDataTissueList <- F

#make dictionary
dict_joined <- read_tsv("dict/dict_joined.txt") %>%
  select(EntrezID_human=entrezID_human, Symbol_human, EntrezID_mouse=entrezID_mouse, Symbol_mouse)

#distinguish between internal and external verson
if("hs_gw_zuber_v2" %in% (libraries %>% collect %>% .$library_id)){
  view <- "internal"
  dataset_selection_all <- list("dropout" = "dropout", "drug_modifier" = "drug_modifier", "facs_based" = "facs_based")
}else{
  view <- "external"
  dataset_selection_all <- list("dropout" = "dropout")
}


