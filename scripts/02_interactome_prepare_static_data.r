#!/usr/bin/env Rscript

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

set.seed(444)
setwd("/gpfs/projects/bsc08/shared_projects/DB_perfilesDifusion/multilayer/layers")

# Multilayer no directed & unweighted
cocitation_table <- read.table("cocitation_layer.txt", header = TRUE, sep = "\t") # gene1 gene2 weight
functional_table <- read.table("functional_layer.txt", header = TRUE, sep = "\t") # gene1 gene2
PPI_table <- read.table("ppi_layer.txt", header = TRUE, sep = "\t") # gene1 gene2 weight1 weight2

cocitation_graph <- graph_from_data_frame(cocitation_table[, 1:2], directed = FALSE)
functional_graph <- graph_from_data_frame(functional_table, directed = FALSE)
PPI_graph <- graph_from_data_frame(PPI_table[, 1:2], directed = FALSE)

### We create a 3-layers Multiplex object ###
GGI_PPI_Multiplex <- create.multiplex(list(
    GGI_cocitation = cocitation_graph,
    GGI_functional = functional_graph,
    PPI = PPI_graph
))
AdjMatrix <- compute.adjacency.matrix(GGI_PPI_Multiplex)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)

setwd("/gpfs/home/bsc/bsc232415/profiles")
save(GGI_PPI_Multiplex, AdjMatrixNorm, file = "static_data.RData")
ls()
