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

  library(parallel)
})
print("###############################################")
print("### DBG interactome ####")
print("###############################################")
set.seed(444)
print("###############################################")
getwd()
setwd("/gpfs/projects/bsc08/shared_projects/DB_perfilesDifusion/multilayer/layers") 
getwd()
cocitation_table <- read.table("cocitation_layer.txt", header=TRUE, sep="\t") # 
functional_table <- read.table("functional_layer.txt", header=TRUE, sep="\t") # 
PPI_table <- read.table("ppi_layer.txt", header=TRUE, sep="\t") # 
ls()

num_cores <- floor(detectCores()*0.8)
run_RWR <- function(seed, norm_matrix, multiplex_object) {
    RWR_result <- Random.Walk.Restart.Multiplex(norm_matrix, multiplex_object, seed)
    RWR_df <- RWR_result$RWRM_Results
    return(RWR_df)
}

#### COCITATION ####
cocitation_nodes <- c(cocitation_table$gene1, cocitation_table$gene2)
cocitation_nodes <- unique(cocitation_nodes)
cocitation_graph <- graph_from_data_frame(cocitation_table, directed = FALSE)
GGI_cocitation_MultiplexObject <- create.multiplex(list(GGI=cocitation_graph))

AdjMatrix_GGI_cocitation <- compute.adjacency.matrix(GGI_cocitation_MultiplexObject)
AdjMatrixNorm_GGI_cocitation <- normalize.multiplex.adjacency(AdjMatrix_GGI_cocitation)

list_RWR_GGI_cocitation <- mclapply(cocitation_nodes, function(seed) {
    run_RWR(seed, AdjMatrixNorm_GGI_cocitation, GGI_cocitation_MultiplexObject)
}, mc.cores = num_cores)
names(list_RWR_GGI_cocitation) <- cocitation_nodes
save(list_RWR_GGI_cocitation, file = "list_RWR_GGI_cocitation")

#### FUNCTIONAL ####
functional_nodes <- c(functional_table$gene1, functional_table$gene2)
functional_nodes <- unique(functional_nodes)
functional_graph <- graph_from_data_frame(functional_table, directed = FALSE)
GGI_functional_MultiplexObject <- create.multiplex(list(GGI=functional_graph))

AdjMatrix_GGI_functional <- compute.adjacency.matrix(GGI_functional_MultiplexObject)
AdjMatrixNorm_GGI_functional <- normalize.multiplex.adjacency(AdjMatrix_GGI_functional)

list_RWR_GGI_functional <- mclapply(functional_nodes, function(seed) {
    run_RWR(seed, AdjMatrixNorm_GGI_functional, GGI_functional_MultiplexObject)
}, mc.cores = num_cores)
names(list_RWR_GGI_functional) <- functional_nodes
save(list_RWR_GGI_functional, file = "list_RWR_GGI_functional")

#### PPI ####
PPI_nodes <- c(PPI_table$gene1, PPI_table$gene2)
PPI_nodes <- unique(PPI_nodes)
PPI_graph <- graph_from_data_frame(PPI_table, directed = FALSE)
PPI_MultiplexObject <- create.multiplex(list(PPI=PPI_graph))

AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

list_RWR_PPI <- mclapply(PPI_nodes, function(seed) {
    run_RWR(seed, AdjMatrixNorm_PPI, PPI_MultiplexObject)
}, mc.cores = num_cores)
names(list_RWR_PPI) <- PPI_nodes
save(list_RWR_PPI, file = "list_RWR_PPI")


### We create a 3-layers Multiplex object ###
GGI_PPI_Multiplex <- create.multiplex(list(GGI_cocitation=cocitation_graph,GGI_functional=functional_graph,PPI=PPI_graph))
nodes <- GGI_PPI_Multiplex$Pool_of_Nodes
head(nodes)
AdjMatrix <- compute.adjacency.matrix(GGI_PPI_Multiplex)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)

run_RWR_tau <- function(seed, norm_matrix, multiplex_object, tau) {
    RWR_result <- Random.Walk.Restart.Multiplex(norm_matrix, multiplex_object, seed, tau = tau)
    RWR_df <- RWR_result$RWRM_Results
    return(RWR_df)
}

tau_combinations <- list()

for (i in 0:10) {
    for (j in 0:(10-i)) {
        k <- 10 - i - j 
        tau_combinations[[length(tau_combinations) + 1]] <- c(i, j, k)*3/10 # suma / nº capas = 1
    }
}

setwd("/gpfs/home/bsc/bsc232415/profiles")
getwd()

for (tau in tau_combinations) {
    result_tau <- mclapply(nodes, function(seed) {
        run_RWR_tau(seed, AdjMatrixNorm, GGI_PPI_MultiplexObject, tau)
    }, mc.cores = num_cores)
    names(result_tau) <- nodes
    result_tau$tau <- tau
    tau_string <- paste(tau, collapse = "_") 
    save(result_tau, file = paste0("profile_",tau_string,".RData"))
} 
