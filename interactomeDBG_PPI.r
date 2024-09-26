suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(reshape)
  library(reshape2)
  library(tidyverse)

  library(RandomWalkRestartMH)
  library(multinet)
  library(igraph)

  library(ggplot2)

  library(ontologyIndex)
  library(ontologySimilarity)
  library(STRINGdb)
  library(hgu95av2.db)
})
print("###############################################")
print("### DBG interactome ####")
print("###############################################")
set.seed(444)
print("###############################################")
print("PPI Network")
print("###############################################")
string_prot_prot <- read_table("data/9606.protein.links.full.v12.0.txt", col_names = TRUE)
string_prot_prot_filt <- string_prot_prot[which(string_prot_prot$experiments > 0 | string_prot_prot$database > 0), ]
string_prot_prot_filt <- string_prot_prot_filt[, c("protein1", "protein2", "combined_score")]
# Transformación ENSP --> GeneSymbol
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 0)
# STRINGdb da problemas al ejecutarlo mediante job, no conecta con su servidor,ejecutarlo directamente en terminal
names1 <- string_db$add_proteins_description(data.frame(STRING_id = string_prot_prot_filt$protein1))
names1 <- names1[, 1:2]
names2 <- string_db$add_proteins_description(data.frame(STRING_id = string_prot_prot_filt$protein2))
names2 <- names2[, 1:2]
names_list <- setNames(as.list(names2$preferred_name), names2$STRING_id)

names2$STRING_id <- string_prot_prot_filt$protein2
names2$preferred_name <- unlist(names_list[string_prot_prot_filt$protein2])

PPI_df <- data.frame(
  ENSP1 = names1$STRING_id,
  PROT1 = names1$preferred_name,
  ENSP2 = names2$STRING_id,
  PROT2 = names2$preferred_name,
  weight = string_prot_prot_filt$combined_score
)
PPI_df$ENSP1 <- sub("^9606\\.", "", PPI_df$ENSP1)
PPI_df$ENSP2 <- sub("^9606\\.", "", PPI_df$ENSP2)
PPI_df_filtered <- PPI_df %>%
  filter(ENSP1 != PROT1 & ENSP2 != PROT2)
# distinct_pairs <- PPI_df %>% distinct(ENSP1, PROT1, ENSP2, PROT2)
aux <- PPI_df_filtered[, c("PROT1", "PROT2", "weight")]
save(aux,file = "data/proteins_matrix")
prot_prot_graph <- graph_from_data_frame(PPI_df_filtered[, c("PROT1", "PROT2", "weight")], directed = FALSE)
PPI_MultiplexObject <- create.multiplex(list(PPI = prot_prot_graph))
PPI_MultiplexObject
save(PPI_MultiplexObject, file = "data/networks/PPI")
