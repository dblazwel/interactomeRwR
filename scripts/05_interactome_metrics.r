suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(parallel)
  library(lsa)
})

set.seed(444)

setwd("/gpfs/scratch/bsc08/bsc232415")
getwd()
load("results/tau/Drug_profiles_profile_5_3_2.RData")
load("results/tau/Disease_profiles_profile_5_3_2.RData")

dis_gene_matrix[is.na(dis_gene_matrix)] <- 0
drug_gene_matrix[is.na(drug_gene_matrix)] <- 0

drug_gene_list <- as.list(as.data.frame(t(drug_gene_matrix)))
dis_gene_list <- as.list(as.data.frame(t(dis_gene_matrix)))

calculate_similarity <- function(disease_name, disease_profile, drug_list) {
  result <- lapply(names(drug_list), function(drug_name) {
    drug_profile <- drug_list[[drug_name]]
    if (all(disease_profile == 0) || all(drug_profile == 0)) {
      return(data.frame(
        Disease = disease_name,
        Drug = drug_name,
        Cosine = NA,
        Pearson = NA,
        Pearson_pval = NA,
        Spearman = NA,
        Spearman_pval = NA,
        L1 = NA,
        L2 = NA
      ))
    }

    cosine_sim <- cosine(disease_profile, drug_profile)
    pearson_test <- tryCatch(cor.test(disease_profile, drug_profile, method = "pearson"), error = function(e) NULL)
    spearman_test <- tryCatch(cor.test(disease_profile, drug_profile, method = "spearman"), error = function(e) NULL)

    pearson_sim <- if (!is.null(pearson_test)) pearson_test$estimate else NA
    pearson_pval <- if (!is.null(pearson_test)) pearson_test$p.value else NA
    spearman_sim <- if (!is.null(spearman_test)) spearman_test$estimate else NA
    spearman_pval <- if (!is.null(spearman_test)) spearman_test$p.value else NA
    l1_sim <- sum(abs(disease_profile - drug_profile)) / length(disease_profile)
    l2_sim <- sqrt(sum((disease_profile - drug_profile)^2)) / length(disease_profile)
    data.frame(
      Disease = disease_name,
      Drug = drug_name,
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
}

num_cores <- max(1, detectCores() - 2) # Usa menos núcleos
chunk_size <- ceiling(length(dis_gene_list) / 25) 
disease_chunks <- split(names(dis_gene_list), (seq_len(length(dis_gene_list)) - 1) %/% chunk_size)
length(disease_chunks) # 223

results <- mclapply(disease_chunks, function(disease_chunk) {
  chunk_results <- lapply(disease_chunk, function(disease_name) {
    disease_profile <- dis_gene_list[[disease_name]]
    calculate_similarity(disease_name, disease_profile, drug_gene_list)
  })
   rbindlist(chunk_results)
}, mc.cores = num_cores)

similarity_matrix_5_3_2 <- rbindlist(results)
save(similarity_matrix_5_3_2, file = "results/tau/matrix/similarity_matrix_5_3_2.RData")
