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

run_RWR_tau <- function(seed, norm_matrix, multiplex_object, tau) {
  tryCatch(
    {
      RWR_result <- Random.Walk.Restart.Multiplex(norm_matrix, multiplex_object, seed, tau = tau)
      RWR_df <- RWR_result$RWRM_Results
      return(RWR_df)
    },
    error = function(e) {
      message("Error with seed: ", seed, " and tau: ", tau, " - ", e)
      return(NULL)
    }
  )
}

calculate_diffusion_profile <- function(tau1, tau2, tau3) {
  print("### DBG interactome ####")
  print("###############################################")
  cat("### Calculate diffusion profile: ", tau1, tau2, tau3, "\n")
  print("###############################################")

  setwd("/gpfs/home/bsc/bsc232415/profiles")
  load("static_data.RData")

  setwd("/gpfs/scratch/bsc08/bsc232415/profiles")
  num_cores <- max(1, floor(detectCores() * 0.8))
  tau <- c(tau1, tau2, tau3)
  nodes <- GGI_PPI_Multiplex$Pool_of_Nodes

  result_tau <- mclapply(nodes, function(seed) {
    run_RWR_tau(seed, AdjMatrixNorm, GGI_PPI_Multiplex, tau)
  }, mc.cores = num_cores)
  names(result_tau) <- nodes
  result_tau$tau <- tau
  tau_string <- paste(tau * 10 / 3, collapse = "_")
  tau_string

  save(result_tau, file = paste0("profile_", tau_string, ".RData"))
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Error arguments: tau1, tau2, tau3")
}

tau1 <- as.numeric(args[1])
tau2 <- as.numeric(args[2])
tau3 <- as.numeric(args[3])

calculate_diffusion_profile(tau1, tau2, tau3)
