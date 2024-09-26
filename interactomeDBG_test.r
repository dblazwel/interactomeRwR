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
print("### DBG interactome - Profiles GGI ####")
print("###############################################")
set.seed(444)
print("###############################################")
load("data/networks/PPI") # PPI_MultiplexObject
load("data/networks/GGI") # GGI_MultiplexObject
ls()

GGI_nodes <- GGI_MultiplexObject$Pool_of_Nodes
PPI_nodes <- PPI_MultiplexObject$Pool_of_Nodes
common_nodes <- intersect(GGI_nodes, PPI_nodes)

AdjMatrix_GGI <- compute.adjacency.matrix(GGI_MultiplexObject)
AdjMatrixNorm_GGI <- normalize.multiplex.adjacency(AdjMatrix_GGI)

AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

# Crear una lista para almacenar los resultados de cada nodo
list_RWR <- list()
for (seed in common_nodes) {
    RWR_result <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_MultiplexObject, seed)
    RWR_df <- RWR_result$RWRM_Results
    list_RWR[[seed]] <- RWR_df
}

list_RWR
save(list_RWR, file = "list_RWR_PPI")
