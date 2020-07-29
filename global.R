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

con_facs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "screen_facs.db")


pheno <- con %>%
  tbl("pheno") %>%
  collect()

libraries <- pheno %>%
  select(library_id, species) %>%
  distinct %>%
  collect()

species <- pheno %>%
  select(species) %>%
  distinct %>%
  .$species

features <- con %>%
  tbl("features") %>%
  select(guide_id, gene_id, hgnc_symbol, entrez_id, library_id, sequence) %>%
  filter(gene_id != "AMBIGUOUS") %>%
  filter(gene_id != "UNMAPPED") %>%
  filter(gene_id != "NOFEATURE") %>%
  filter(gene_id != "SAFETARGETING") %>%
  filter(gene_id != "NONTARGETING") %>%
  filter(!is.na(gene_id)) %>%
  distinct %>% collect()

features_facs <- con_facs %>%
  tbl("features") %>%
  select(guide_id, gene_id, entrez_id, symbol, sequence) %>%
  distinct %>% collect()

contrasts <- con %>%
  tbl("contrasts") %>%
  collect()

contrasts_facs <- con_facs %>%
  tbl("contrasts") %>%
  collect()

con_sgRNAs <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "sgRNAs.db")

gene_list_human <- con_sgRNAs %>%
  tbl("sgRNAs_human") %>%
  select(HGNC_approved) %>%
  rename(Symbol=HGNC_approved) %>%
  distinct %>%
  arrange(Symbol) %>% 
  collect %>%
  .$Symbol

gene_list_mouse <- con_sgRNAs %>%
  tbl("sgRNAs_mouse") %>%
  select(Symbol=`MGI name`) %>%
  distinct %>%
  arrange(Symbol) %>%
  collect() %>%
  .$Symbol

#make dictionary
dict_joined <- read_tsv("dict/dict_joined.txt")


# sgRNAs <- con_sgRNAs %>%
#   tbl("sgRNAs_human") %>%
#   rename(Symbol = HGNC_approved, `23mer` = sgRNA_23mer,  Sequence = cloning_sgRNAs_trimmed, Length = len_cloning_sgRNA, maps_to_genome = check) %>%
#   select(-sgRNA_ID) %>%
#   collect()
# 
# sgRNAs_human_genes <- sgRNAs_human %>%
#   select(Symbol) %>%
#   distinct %>%
#   arrange(Symbol) %>%
#   .$Symbol
# 
# sgRNAs_mouse <- con_sgRNAs %>%
#   tbl("sgRNAs_mouse") %>%
#   rename(Symbol = MGI.name, `23mer` = sgRNA_23mer,  Sequence = cloning_sgRNAs_trimmed, Length = len_cloning_sgRNA, maps_to_genome = check) %>%
#   select(-sgRNA_ID) %>%
#   collect()
# 
# sgRNAs_mouse_genes <- sgRNAs_mouse %>%
#   select(Symbol) %>%
#   distinct %>%
#   arrange(Symbol) %>%
#   .$Symbol

