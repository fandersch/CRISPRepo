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
library(shinycssloaders)
library(shinyjs)
library(tidyverse)
library(DT)
library(rlang)
library(DBI)
library(zip)
library(shinyBS)
library(Biostrings)
library(BayesFactor)
library(plotly)
library(pheatmap)
library(EnhancedVolcano)

#back button click handling
jscode <- 'window.onbeforeunload = function() { return "Please use the button on the webpage"; };'

#activate with 'modify_stop_propagation(menuItem(...))' in ui.R
modify_stop_propagation <- function(x) {
  x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}

#set default packages for functions
renderDataTable <- DT::renderDataTable

con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/screen.db")

con_expression <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/expression_data_counts_tpm_tmm.db")

con_slamseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/slamseq.db")

con_quantseq <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/quantseq.db")

con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/sgRNAs.db")

con_correlations <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/correlations.db")

con_correlations_tissue <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/correlations_tissue.db")

con_cell_lines <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cell_line_meta_data.db")

con_patient_data <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "databases/cBioportal_mutations_CNAs_fusions.db")

pheno <- con %>%
  tbl("pheno") %>%
  collect()

species <- pheno %>%
  dplyr::select(species) %>%
  distinct %>%
  .$species

features <- con %>%
  tbl("features") %>%
  dplyr::select(guide_id, gene_id, symbol, entrez_id, sequence, sequence_matching, id_entrez_23mer, context, library_id) %>%
  dplyr::filter(gene_id != "AMBIGUOUS") %>%
  dplyr::filter(gene_id != "UNMAPPED") %>%
  dplyr::filter(gene_id != "NOFEATURE") %>%
  dplyr::filter(gene_id != "SAFETARGETING") %>%
  dplyr::filter(gene_id != "NONTARGETING") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct %>%
  collect()

contrasts <- con %>%
  tbl("contrasts") %>%
  dplyr::select(contrast_id, contrast_id_QC, library_id, cellline_name, tissue_name, species, type, reference_type, dynamic_range, auc, treatment, control) %>%
  collect()

tissue_list_screens <- contrasts %>%
  dplyr::filter(type == "dropout") %>%
  dplyr::select(tissue_name) %>%
  distinct() %>%
  arrange(tissue_name) %>%
  .$tissue_name

libraries <- contrasts %>%
  dplyr::select(library_id, cellline_name, tissue_name, species, type, reference_type) %>%
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
  arrange(symbol) %>%
  collect()

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
#expression data
cellline_list_expressionData <- con_expression %>%
  tbl("expression_data_meta_info") %>%
  dplyr::select(sample_id, cell_line_name, tissue_name, species) %>%
  distinct() %>%
  arrange(cell_line_name) %>%
  collect()

gene_list_expressionData <- con_expression %>%
  tbl("expression_data_genes") %>%
  collect()

#slamseq
sample_list_slamseq <- con_slamseq %>%
  tbl("slamseq_sample_metadata") %>%
  dplyr::rename(cell_line_name = cell_line, tissue_name=tissue) %>%
  distinct() %>%
  arrange(cell_line_name) %>%
  collect() %>%
  mutate(sample_id = as.character(sample_id),
         dataset = "SLAM-seq") %>%
  select(dataset, everything())

deseq2_list_slamseq <- con_slamseq %>%
  tbl("slamseq_deseq2_metadata") %>%
  distinct() %>%
  collect() %>%
  mutate(dataset = "SLAM-seq") %>%
  select(dataset, everything())

gene_list_slamseq <- con_slamseq %>%
  tbl("slamseq_genes") %>%
  collect()

#quantseq
sample_list_quantseq <- con_quantseq %>%
  tbl("quantseq_sample_metadata") %>%
  dplyr::rename(cell_line_name = cell_line, tissue_name=tissue) %>%
  distinct() %>%
  arrange(cell_line_name) %>%
  collect() %>%
  mutate(sample_id = as.character(sample_id),
         dataset = "Quant-seq") %>%
  select(dataset, everything())

deseq2_list_quantseq <- con_quantseq %>%
  tbl("quantseq_deseq2_metadata") %>%
  distinct() %>%
  collect() %>%
  mutate(dataset = "Quant-seq") %>%
  select(dataset, everything())

gene_list_quantseq <- con_quantseq %>%
  tbl("quantseq_genes") %>%
  collect()

loadExpressionDataTissueList <- F

#cell line
cellline_list_cellLine <- con_cell_lines %>%
  tbl("cell_line_meta") %>%
  dplyr::select(cell_line_name, tissue_name, cell_line_id, SangerCellModelPassports_cancer_type, SangerCellModelPassports_cancer_type_detail,
                BroadDepMap_OncotreeSubtype, BroadDepMap_OncotreePrimaryDisease) %>%
  distinct() %>%
  arrange(cell_line_name) %>%
  collect

gene_list_cellLine <- con_cell_lines %>%
  tbl("cell_line_genes") %>%
  collect %>%
  arrange(symbol)

tissue_list_cellLine <- cellline_list_cellLine %>%
  dplyr::select(tissue_name) %>%
  distinct() %>%
  arrange(tissue_name) %>%
  .$tissue_name

cn_category_cellLine <- con_cell_lines %>%
  tbl("cell_line_gene_cnv") %>%
  dplyr::select(cn_category) %>%
  distinct() %>%
  collect() %>%
  .$cn_category

#patient
patient_cancer_types <- con_patient_data %>%
  tbl("curated_set_non_redundant_sample_info") %>%
  dplyr::select(cancer_type, TCGA_pan_cancer) %>%
  distinct() %>%
  arrange(cancer_type) %>%
  collect()

patient_cancer_types_choices <- patient_cancer_types$cancer_type

patient_cancer_types_pediatric <- con_patient_data %>%
  tbl("curated_set_non_redundant_sample_info") %>%
  dplyr::select(cancer_type_pediatric) %>%
  distinct() %>%
  arrange(cancer_type_pediatric) %>%
  collect()

patient_studies <- con_patient_data %>%
  tbl("curated_set_non_redundant_sample_info") %>%
  dplyr::select(study_name, TCGA_pan_cancer) %>%
  distinct() %>%
  arrange(study_name) %>%
  collect()

patient_studies_pediatric <- con_patient_data %>%
  tbl("curated_set_non_redundant_sample_info") %>%
  dplyr::filter(pediatric) %>%
  dplyr::select(study_name) %>%
  distinct() %>%
  arrange(study_name) %>%
  collect()

patient_genes_all <- con_patient_data %>%
  tbl("genes") %>%
  dplyr::select(cancer_type, cancer_type_pediatric, study_name, entrez_id, symbol) %>%
  distinct() %>%
  arrange(symbol) %>%
  collect()

patient_genes_all_init <- patient_genes_all %>%
  dplyr::mutate(gene = ifelse(is.na(symbol), paste0("No symbol found (" , entrez_id, " )"), paste0(" ", symbol , " ( ", entrez_id, " )"))) %>%
  dplyr::select(gene) %>%
  distinct %>%
  arrange(gene) %>%
  .$gene

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

default_df<-NULL
#distinguish between internal and external version

if("hs_gw_zuber_v2" %in% (libraries %>% collect %>% .$library_id)){
  view <- "internal"
  #load default df
  default_df <- read_tsv(file="essentiality_data/human_scaledLFC_fdr_adjusted_geneLevel_HQscreens.tsv", col_names = T)
}else{
  view <- "external"
}

displayed_table <- NULL

#close db conections
DBI::dbDisconnect(con)
DBI::dbDisconnect(con_expression)
DBI::dbDisconnect(con_slamseq)
DBI::dbDisconnect(con_quantseq)
DBI::dbDisconnect(con_sgRNAs)
DBI::dbDisconnect(con_correlations)
DBI::dbDisconnect(con_correlations_tissue)
DBI::dbDisconnect(con_cell_lines)
DBI::dbDisconnect(con_patient_data)
