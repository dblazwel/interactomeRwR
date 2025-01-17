suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(parallel)
  library(lsa)
})

set.seed(444)

setwd("/gpfs/scratch/bsc08/bsc232415")
getwd()
load("results/Disease_profiles_GGI_PPI.RData")
load("results/Drug_profiles_GGI_PPI.RData")


dis_gene_matrix[is.na(dis_gene_matrix)] <- 0
drug_gene_matrix[is.na(drug_gene_matrix)] <- 0

dis_gene_list <- as.list(as.data.frame(t(dis_gene_matrix)))
drug_gene_list <- as.list(as.data.frame(t(drug_gene_matrix)))

calculate_pairwise_similarity <- function(profile_list, num_chunks = 224, num_cores = 112) {
  entities <- names(profile_list)
  combinations <- expand.grid(Entity1 = entities, Entity2 = entities, stringsAsFactors = FALSE)
  
  
  # Dividir en chunks
  chunk_size <- ceiling(nrow(combinations) / num_chunks) 
  combinations_split <- split(combinations, (seq_len(nrow(combinations)) - 1) %/% chunk_size)
  
  # Procesar cada chunk en paralelo
  results <- mclapply(combinations_split, function(chunk) {
    result <- lapply(1:nrow(chunk), function(idx) {
      entity1 <- chunk$Entity1[idx]
      entity2 <- chunk$Entity2[idx]
      profile1 <- profile_list[[entity1]]
      profile2 <- profile_list[[entity2]]
      
      if (all(profile1 == 0) || all(profile2 == 0)) {
        return(data.frame(
          Entity1 = entity1,
          Entity2 = entity2,
          Cosine = NA,
          Pearson = NA,
          Pearson_pval = NA,
          Spearman = NA,
          Spearman_pval = NA,
          L1 = NA,
          L2 = NA
        ))
      }
      
      cosine_sim <- cosine(profile1, profile2)
      pearson_test <- tryCatch(cor.test(profile1, profile2, method = "pearson"), error = function(e) NULL)
      spearman_test <- tryCatch(cor.test(profile1, profile2, method = "spearman"), error = function(e) NULL)
      
      pearson_sim <- if (!is.null(pearson_test)) pearson_test$estimate else NA
      pearson_pval <- if (!is.null(pearson_test)) pearson_test$p.value else NA
      spearman_sim <- if (!is.null(spearman_test)) spearman_test$estimate else NA
      spearman_pval <- if (!is.null(spearman_test)) spearman_test$p.value else NA
      l1_sim <- sum(abs(profile1 - profile2)) / length(profile1)
      l2_sim <- sqrt(sum((profile1 - profile2)^2)) / length(profile1)
      
      data.frame(
        Entity1 = entity1,
        Entity2 = entity2,
        Cosine = cosine_sim,
        Pearson = pearson_sim,
        Pearson_pval = pearson_pval,
        Spearman = spearman_sim,
        Spearman_pval = spearman_pval,
        L1 = l1_sim,
        L2 = l2_sim
      )
    })
    rbindlist(result)
  }, mc.cores = num_cores)
  
  # Combinar resultados de los chunks
  rbindlist(results)
}

cat("Calculando matriz enfermedad-enfermedad...\n")
# disease_similarity_matrix <- calculate_pairwise_similarity(dis_gene_list, num_chunks = 25, num_cores = 112)
# save(disease_similarity_matrix, file = "results/disease_similarity_matrix.RData")

cat("Calculando matriz fármaco-fármaco...\n")
drug_similarity_matrix <- calculate_pairwise_similarity(drug_gene_list, num_chunks = 25, num_cores = 112)
save(drug_similarity_matrix, file = "results/drug_similarity_matrix.RData")
