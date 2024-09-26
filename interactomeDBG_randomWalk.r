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
  library(Cairo)

  library(AnnotationDbi)
  library(hgu95av2.db)
})
print("###############################################")
print("### DBG interactome ####")
print("###############################################")
set.seed(444)
print("###############################################")
load("data/networks/PPI") # PPI_MultiplexObject
load("data/networks/GGI") # GGI_MultiplexObject
# load("data/networks/GGI_PPI") # GGI_PPI_network
# load("data/networks/TransitionMatrix") # TransitionMatrix
ls()

## DRUG/DISEASE
# ICD10_disease_drug <- read_csv("data/ICD10_disease_drug.csv")
# gene_disease_rel <- read_csv("data/gene_disease_rel.csv")

# gene_list_original <- unique(matrix_melt$GENE1)
# gene_list_filt <- unique(gene_disease_rel$GeneSymbol)
# genes_lost <- setdiff(gene_list_original, gene_list_filt)

# print(paste0("Genes originales: ", length(gene_list_original), " / Genes mapeados: ", length(gene_list_filt), " / Genes perdidos: ", length(genes_lost)))
# print('###############################################')
# genes_lost[1:30]

get_network <- function(GGI, PPI) {
  GGI_nodes <- GGI_MultiplexObject$Pool_of_Nodes
  PPI_nodes <- PPI_MultiplexObject$Pool_of_Nodes
  binary_matrix <- matrix(0, nrow = length(GGI_nodes), ncol = length(PPI_nodes))
  common_nodes <- intersect(GGI_nodes, PPI_nodes)
  rownames(binary_matrix) <- GGI_nodes
  colnames(binary_matrix) <- PPI_nodes
  for (node in common_nodes) {
    i <- which(GGI_nodes == node)
    j <- which(PPI_nodes == node)
    binary_matrix[i, j] <- 1
  }
  melted_df <- melt(binary_matrix)
  colnames(melted_df) <- c("GGI_node", "PPI_node", "weight")
  GGI_PPI_rel <- as.data.frame(melted_df)

  GGI_PPI_network <- create.multiplexHet(
    GGI_MultiplexObject,
    PPI_MultiplexObject, GGI_PPI_rel
  )

  save(GGI_PPI_network, file = "data/networks/GGI_PPI")
  return(GGI_PPI_network)
}

GGI_PPI_network <- get_network(GGI_MultiplexObject, PPI_MultiplexObject)
TransitionMatrix <- compute.transition.matrix(GGI_PPI_network)
colnames(TransitionMatrix) <- gsub("_1$", "", colnames(TransitionMatrix))
rownames(TransitionMatrix) <- gsub("_1$", "", rownames(TransitionMatrix))
save(TransitionMatrix, file = "data/networks/TransitionMatrix")

GGI_nodes <- GGI_MultiplexObject$Pool_of_Nodes
PPI_nodes <- PPI_MultiplexObject$Pool_of_Nodes
binary_matrix <- matrix(0, nrow = length(GGI_nodes), ncol = length(PPI_nodes))
common_nodes <- intersect(GGI_nodes, PPI_nodes)
rownames(binary_matrix) <- GGI_nodes
colnames(binary_matrix) <- PPI_nodes
for (node in common_nodes) {
  i <- which(GGI_nodes == node)
  j <- which(PPI_nodes == node)
  binary_matrix[i, j] <- 1
}
melted_df <- melt(binary_matrix)
colnames(melted_df) <- c("GGI_node", "PPI_node", "weight")
GGI_PPI_rel <- as.data.frame(melted_df)

seedGenes <- c("IL7R")
seedProt <- c("ARF5")
RWRH_Results <-
  Random.Walk.Restart.MultiplexHet(
    TransitionMatrix,
    GGI_PPI_network, seedGenes, seedProt
  )
RWRH_Results
TopResults <-
  create.multiplexHetNetwork.topResults(RWRH_Results,
    GGI_PPI_network, GGI_PPI_rel,
    k = 10
  )


