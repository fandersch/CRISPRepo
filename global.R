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
library(DBI)
library(zip)
library(shinyBS)

#set default packages for functions
renderDataTable <- DT::renderDataTable

con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")

con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm.db")

con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")

con_correlations <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/correlations.db")

con_correlations_tissue <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/correlations_tissue.db")

con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")

pheno <- con %>%
  tbl("pheno")

species <- pheno %>%
  dplyr::select(species) %>%
  distinct %>%
  .$species

features <- con %>%
  tbl("features") %>%
  dplyr::select(guide_id, gene_id, symbol, entrez_id, sequence, sequence_matching, sgRNA_23mer, context, library_id) %>%
  dplyr::filter(gene_id != "AMBIGUOUS") %>%
  dplyr::filter(gene_id != "UNMAPPED") %>%
  dplyr::filter(gene_id != "NOFEATURE") %>%
  dplyr::filter(gene_id != "SAFETARGETING") %>%
  dplyr::filter(gene_id != "NONTARGETING") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct

contrasts <- con %>%
  tbl("contrasts") %>%
  dplyr::select(contrast_id, contrast_id_QC, library_id, cellline_name, tissue_name, species, type, reference_type, dynamic_range, auc, treatment, control)

libraries <- contrasts %>%
  dplyr::select(library_id, cellline_name, tissue_name, species, type) %>%
  distinct 

gene_list_screens <- con %>%
  tbl("features") %>%
  dplyr::filter(gene_id != "AMBIGUOUS") %>%
  dplyr::filter(gene_id != "UNMAPPED") %>%
  dplyr::filter(gene_id != "NOFEATURE") %>%
  dplyr::filter(gene_id != "SAFETARGETING") %>%
  dplyr::filter(gene_id != "NONTARGETING") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  dplyr::select(gene_id, symbol, entrez_id, library_id) %>%
  distinct %>%
  arrange(symbol)

gene_list_human <- NULL
gene_list_mouse <- NULL 

if("hs_gw_zuber_v2" %in% (libraries %>% collect %>% .$library_id)){
  gene_list_human <- con_sgRNAs %>%
    tbl("genes_human") %>%
    dplyr::select(Symbol, EntrezID) %>%
    distinct %>%
    arrange(Symbol) %>%
    collect

  gene_list_mouse <- con_sgRNAs %>%
    tbl("genes_mouse") %>%
    dplyr::select(Symbol, EntrezID) %>%
    distinct %>%
    arrange(Symbol) %>%
    collect
}

cellline_list_expressionData <- con_expression %>%
  tbl("expression_data_meta_info") %>%
  dplyr::select(sample_id, cell_line_name = cell_line_name_stripped, tissue_name, species) %>%
  distinct() %>%
  arrange(cell_line_name)

gene_list_expressionData <- con_expression %>%
  tbl("expression_data_genes")

loadExpressionDataTissueList <- F

#cell line
cellline_list_cellLine <- con_cell_lines %>%
  tbl("cell_line_meta") %>%
  dplyr::select(cell_line_name, tissue_name, Sanger_model_ID) %>%
  distinct() %>%
  arrange(cell_line_name) %>%
  collect %>%
  dplyr::rename(model_id = Sanger_model_ID)

gene_list_cellLine <- con_cell_lines %>%
  tbl("cell_line_genes") %>%
  collect %>%
  arrange(gene_symbol)

tissue_list_cellLine <- cellline_list_cellLine %>%
  dplyr::select(tissue_name) %>%
  distinct() %>%
  arrange(tissue_name) %>%
  .$tissue_name

#correlations
gene_list_correlations <- con_correlations %>%
  tbl("genes") %>%
  collect %>%
  .$gene

gene_list_correlations_tissue <- con_correlations_tissue %>%
  tbl("genes") %>%
  collect %>%
  .$gene

tissue_list_correlations_tissue <- con_correlations_tissue %>%
  tbl("tissues") %>%
  collect %>%
  .$tissue

tissue_list_correlations_tissue <- c("All", tissue_list_correlations_tissue)

#load dictionary
dict_joined <- read_tsv("dict/dict_joined.txt") %>%
  dplyr::select(EntrezID_human=entrezID_human, Symbol_human, EntrezID_mouse=entrezID_mouse, Symbol_mouse)

#load essentialome
essentialome <- read_tsv("essentialome_file_crisprepo.txt", col_names = T)

all_types <- contrasts %>%
  collect %>%
  dplyr::filter(species %in% "human") %>%
  .$type %>%
  unique

default_df<-NULL
#distinguish between internal and external version

if("hs_gw_zuber_v2" %in% (libraries %>% collect %>% .$library_id)){
  view <- "internal"
  dataset_selection_all <- setNames(all_types, all_types)
  #load default df
  default_df <- read_tsv(file="essentiality_data/human_scaledLFC_fdr_adjusted_geneLevel_HQscreens.tsv", col_names = T)
}else{
  view <- "external"
  dataset_selection_all <- list("dropout" = "dropout")
}

displayed_table <- NULL