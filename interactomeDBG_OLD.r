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
load("data/networks/GGI_PPI") # GGI_PPI_network
load("data/networks/TransitionMatrix") # TransitionMatrix
ls()

GGI_nodes <- GGI_MultiplexObject$Pool_of_Nodes
PPI_nodes <- PPI_MultiplexObject$Pool_of_Nodes

# Crear la matriz binaria y verificar nodos comunes
binary_matrix <- matrix(0, nrow = length(GGI_nodes), ncol = length(PPI_nodes))
common_nodes <- intersect(GGI_nodes, PPI_nodes)
rownames(binary_matrix) <- GGI_nodes
colnames(binary_matrix) <- PPI_nodes

for (node in common_nodes) {
  i <- which(GGI_nodes == node)
  j <- which(PPI_nodes == node)
  binary_matrix[i, j] <- 1
}

# Convertir la matriz a data frame
melted_df <- melt(binary_matrix)
colnames(melted_df) <- c("GGI_node", "PPI_node", "weight")
GGI_PPI_rel <- as.data.frame(melted_df)

# Filtrar filas con weight == 1 para asegurar nodos relevantes
# GGI_PPI_rel <- GGI_PPI_rel[GGI_PPI_rel$weight == 1, ]

GGI_PPI_network <- create.multiplexHet(
    GGI_MultiplexObject,
    PPI_MultiplexObject, GGI_PPI_rel
)
#save(GGI_PPI_network, file = "data/networks/GGI_PPI")
TransitionMatrix <- compute.transition.matrix(GGI_PPI_network)
colnames(TransitionMatrix) <- gsub("_1$", "", colnames(TransitionMatrix))
rownames(TransitionMatrix) <- gsub("_1$", "", rownames(TransitionMatrix))
#save(TransitionMatrix, file = "data/networks/TransitionMatrix")

seedGenes <- c("IL7R")
seedProt <- c("ARF5")

# Verificar si los seedGenes y seedProt están presentes en TransitionMatrix
if (!all(seedGenes %in% rownames(TransitionMatrix))) {
  stop("Some seedGenes are not present in the TransitionMatrix")
}

if (!all(seedProt %in% colnames(TransitionMatrix))) {
  stop("Some seedProt are not present in the TransitionMatrix")
}

# Ejecutar el Random Walk Restart
RWRH_Results <- Random.Walk.Restart.MultiplexHet(
  TransitionMatrix,
  GGI_PPI_network, seedGenes, seedProt
)
RWRH_Results

# Verificar la estructura de GGI_PPI_rel antes de pasarlo a la función
if (ncol(GGI_PPI_rel) != 3) {
  stop("GGI_PPI_rel should have exactly three columns: GGI_node, PPI_node, weight")
}

# Ejecutar create.multiplexHetNetwork.topResults
TopResults <- create.multiplexHetNetwork.topResults(
  RWRH_Results,
  GGI_PPI_network, GGI_PPI_rel,
  k = 10
)

# Verificar y mostrar los resultados
TopResults

## We print that cluster with its interactions.
CairoPNG("profile_001.png", width = 8, height = 6, units = "in", res = 300)
plot(TopResults,
  vertex.label.color = "black",
  vertex.frame.color = "#ffffff",
  vertex.size = 20, edge.curved = .2,
  vertex.color = ifelse(V(TopResults)$name == "IL7R" |
    V(TopResults)$name == "269880", "yellow",
  ifelse(V(TopResults)$name %in%
    GGI_PPI_network$Multiplex1$Pool_of_Nodes, "#00CCFF", "Grey75")
  ),
  edge.color = ifelse(E(TopResults)$type == "PPI", "blue",
    ifelse(E(TopResults)$type == "Disease", "black", "grey50")
  ),
  edge.width = 0.8,
  edge.lty = ifelse(E(TopResults)$type == "bipartiteRelations",
    2, 1
  ),
  vertex.shape = ifelse(V(TopResults)$name %in%
    GGI_PPI_network$Multiplex1$Pool_of_Nodes, "circle", "rectangle")
)
dev.off()
