suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(pROC)
})

setwd("/gpfs/scratch/bsc08/bsc232415")
getwd()
load("results/similarity_matrix_cocitation.RData")
#load("results/tau/matrix/similarity_matrix_0_3_7.RData")

similarity_matrix <- similarity_matrix_cocitation

#########################################################################################
### Get the drugs and diseases with less than 1,000 genes associated ###
druggen <- fread("associations/Drug_gene_associations.txt", stringsAsFactors = F, sep = "\t")
druggen$GeneSymbol <- toupper(druggen$GeneSymbol)
setkey(druggen, "ChemicalName")

disgen <- fread("associations/Disease_gene_associations.txt", stringsAsFactors = F, sep = "\t")
disgen$GeneSymbol <- toupper(disgen$GeneSymbol)
setkey(disgen, "DiseaseID")

dd <- fread("associations/Drug_disease_associations.txt", stringsAsFactors = F, sep = "\t")
valid_drugs <- intersect(dd$ChemicalName, druggen$ChemicalName)
valid_diseases <- intersect(dd$DiseaseID, disgen$DiseaseID)
filtered_dd <- dd %>%
    filter(ChemicalName %in% valid_drugs & DiseaseID %in% valid_diseases)

# **1. Tabla de AUC para Enfermedades**
# Seleccionar enfermedades con al menos un fármaco relacionado
results <- list()
for (disease in valid_diseases) {
    tryCatch(
        {
            predictors_cos <- similarity_matrix$Cosine[similarity_matrix$Disease == disease]
            predictors_pear <- similarity_matrix$Pearson[similarity_matrix$Disease == disease]
            predictors_spe <- similarity_matrix$Spearman[similarity_matrix$Disease == disease]
            predictors_l1 <- similarity_matrix$L1[similarity_matrix$Disease == disease]
            predictors_l2 <- similarity_matrix$L2[similarity_matrix$Disease == disease]

            associated_drugs <- filtered_dd$ChemicalName[filtered_dd$DiseaseID == disease]
            binary_response <- ifelse(similarity_matrix$Drug[similarity_matrix$Disease == disease] %in% associated_drugs, 1, 0)

            roc_object_cos <- roc(
                response = binary_response,
                predictor = predictors_cos
            )
            roc_object_pear <- roc(
                response = binary_response,
                predictor = predictors_pear
            )
            roc_object_spe <- roc(
                response = binary_response,
                predictor = predictors_spe
            )
            roc_object_l1 <- roc(
                response = binary_response,
                predictor = predictors_l1
            )
            roc_object_l2 <- roc(
                response = binary_response,
                predictor = predictors_l2
            )
            results[[disease]] <- data.frame(
                Disease = disease,
                AUC_Cosine = as.numeric(auc(roc_object_cos)),
                AUC_Pearson = as.numeric(auc(roc_object_pear)),
                AUC_Spearman = as.numeric(auc(roc_object_spe)),
                AUC_L1 = as.numeric(auc(roc_object_l1)),
                AUC_L2 = as.numeric(auc(roc_object_l2))
            )
        },
        error = function(e) {
            results[[disease]] <- data.frame(
                Disease = disease,
                AUC_Cosine = 0,
                AUC_Pearson = 0,
                AUC_Spearman = 0,
                AUC_L1 = 0,
                AUC_L2 = 0
            )
        }
    )
}
auc_diseases <- do.call(rbind, results)
save(auc_diseases, file = "results/AUC/AUC_diseases_cocitation_all_metrics.RData")

# **2. Tabla de AUC para Fármacos**
# Seleccionar fármacos con al menos una enfermedad asociada
results <- list()
for (drug in valid_drugs) {
    tryCatch(
        {
            predictors_cos <- similarity_matrix$Cosine[similarity_matrix$Disease == disease]
            predictors_pear <- similarity_matrix$Pearson[similarity_matrix$Disease == disease]
            predictors_spe <- similarity_matrix$Spearman[similarity_matrix$Disease == disease]
            predictors_l1 <- similarity_matrix$L1[similarity_matrix$Disease == disease]
            predictors_l2 <- similarity_matrix$L2[similarity_matrix$Disease == disease]

            associated_diseases <- filtered_dd$DiseaseID[filtered_dd$ChemicalName == drug]
            binary_response <- ifelse(similarity_matrix$Disease[similarity_matrix$Drug == drug] %in% associated_diseases, 1, 0)
            roc_object_cos <- roc(
                response = binary_response,
                predictor = predictors_cos
            )
            roc_object_pear <- roc(
                response = binary_response,
                predictor = predictors_pear
            )
            roc_object_spe <- roc(
                response = binary_response,
                predictor = predictors_spe
            )
            roc_object_l1 <- roc(
                response = binary_response,
                predictor = predictors_l1
            )
            roc_object_l2 <- roc(
                response = binary_response,
                predictor = predictors_l2
            )
            results[[disease]] <- data.frame(
                Disease = disease,
                AUC_Cosine = as.numeric(auc(roc_object_cos)),
                AUC_Pearson = as.numeric(auc(roc_object_pear)),
                AUC_Spearman = as.numeric(auc(roc_object_spe)),
                AUC_L1 = as.numeric(auc(roc_object_l1)),
                AUC_L2 = as.numeric(auc(roc_object_l2))
            )
        },
        error = function(e) {
            results[[disease]] <- data.frame(
                Disease = disease,
                AUC_Cosine = 0,
                AUC_Pearson = 0,
                AUC_Spearman = 0,
                AUC_L1 = 0,
                AUC_L2 = 0
            )
        }
    )
}
auc_drugs <- do.call(rbind, results)
save(auc_drugs, file = "results/AUC/AUC_drugs_cocitation_all_metrics.RData")
