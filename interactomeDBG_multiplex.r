suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(stringr)
  library(reshape)
  library(reshape2)
  library(tidyverse)

  library(ggplot2)
  library(Cairo)

  library(RandomWalkRestartMH)
  library(igraph)
  library(multinet)

  library(ontologyIndex)
  library(ontologySimilarity)
  library(GO.db)
  library(AnnotationDbi)
  library(hgu95av2.db)
})
print("###############################################")
print("### DBG interactome ####")
print("###############################################")


## We create a 2-layers Multiplex object
gene_matrix_melt <- fread("data/genes_matrix.csv")
gene_gene_graph <- graph_from_data_frame(gene_matrix_melt, directed = FALSE)
load("data/proteins_matrix")
prot_prot_graph <- graph_from_data_frame(aux, directed = FALSE)

GGI_PPI_Multiplex <- create.multiplex(list(GGI=gene_gene_graph,PPI=prot_prot_graph))
GGI_PPI_Multiplex

nodes <- GGI_PPI_Multiplex$Pool_of_Nodes

AdjMatrix <- compute.adjacency.matrix(GGI_PPI_Multiplex)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)

list_RWR <- list()
for (seed in nodes) {
## We launch the algorithm with the default parameters (See details on manual)
    RWR_res <- Random.Walk.Restart.Multiplex(AdjMatrixNorm,
                        GGI_PPI_Multiplex,seed)
    RWR_df <- RWR_res$RWRM_Results
    list_RWR[[seed]] <- RWR_df
}

list_RWR
save(list_RWR, file = "list_RWR_GGIPPI")
