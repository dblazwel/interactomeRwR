suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(parallel)
})

set.seed(444)
setwd("/gpfs/scratch/bsc08/bsc232415")

  #### Get the associations
  tt<-fread("data/CTD_chem_gene_ixns.tsv",stringsAsFactors = F,sep="\t")
  tt<-tt[which(tt$Organism=="Homo sapiens"),]
  if(dim(tt)[2]==11){tt<-tt[,c(1,4,10,11)]}
  nrow(tt) # 984.006
  write.table(tt,"associations/Drug_gene_associations.txt",quote=F,sep="\t",row.names = F)
  
  tt<-fread("data/CTD_curated_genes_diseases.tsv",stringsAsFactors = F,sep="\t")
  nrow(tt) # 34.047
  write.table(tt,"associations/Disease_gene_associations.txt",quote=F,sep="\t",row.names = F)
  
  tt<-fread("data/CTD_chemicals_diseases.tsv",stringsAsFactors = F,sep="\t")
  tt<-tt[which(tt$DirectEvidence=="therapeutic"),c(1:6,10)]
  nrow(tt) # 36.721
  write.table(tt,"associations/Drug_disease_associations.txt",quote=F,sep="\t",row.names = F)

  #########################################################################################
  ### Get the drugs and diseases with less than 1,000 genes associated ###
  druggen<-fread("associations/Drug_gene_associations.txt",stringsAsFactors = F,sep="\t")
  druggen$GeneSymbol<-toupper(druggen$GeneSymbol)
  setkey(druggen,"ChemicalName")

  disgen<-fread("associations/Disease_gene_associations.txt",stringsAsFactors = F,sep="\t")
  disgen$GeneSymbol<-toupper(disgen$GeneSymbol)
  setkey(disgen,"DiseaseID")
  
  dd<-fread("associations/Drug_disease_associations.txt",stringsAsFactors = F,sep="\t")
  
  ## number of drugs with associated genes and viceversa ##
  length(intersect(dd$ChemicalName,druggen$ChemicalName)) # 3.821
  length(intersect(dd$DiseaseID,disgen$DiseaseID)) # 1.563

  drugs<-table(druggen$ChemicalName)
  drugs<-drugs[which(drugs<=1000)]   
  length(drugs) # 9.890
  drug_names <- names(drugs)

  diseases<-table(disgen$DiseaseID)
  diseases<-diseases[which(diseases<=1000)]   
  length(diseases) # 5.852
  disease_names <- names(diseases)

# Calculas los perfiles de difusion para los 9,891 fármacos y 5,853 enfermedades, 

load("data/list_RWR_PPI")
profiles_list <- list_RWR_PPI
genes <- names(profiles_list)


### Drug-gene associations
drug_gene_matrix <- matrix(NA,
    nrow = length(drugs), ncol = length(genes),
    dimnames = list(drug_names, genes)
)

for (drug in drug_names) {
    print(drug)
    related_genes <- druggen %>%
        filter(ChemicalName == drug) %>%
        pull(GeneSymbol)

    difussion_profiles <- list()
    for (gene in related_genes) {
        if (gene %in% genes) {
            difussion_profiles[[gene]] <- profiles_list[[gene]]
        }
    }

    if (length(difussion_profiles) == 1) {
        print("DBG entra por 1")
        unified_profile <- difussion_profiles[[1]]
        matching_nodes <- unified_profile$NodeNames %in% colnames(drug_gene_matrix)
        filtered_profile <- unified_profile[matching_nodes, ]
        if (nrow(filtered_profile) > 0) {
            drug_gene_matrix[drug, filtered_profile$NodeNames] <- filtered_profile$Score
        } else {
            print(paste("No se encontraron nodos coincidentes para el fármaco:", drug))
        }
    }

    if (length(difussion_profiles) > 1) {
        print("DBG entra por >1")
        combined_profiles <- list()
        for (profile in difussion_profiles) {
            combined_profiles[[length(combined_profiles) + 1]] <- profile
        }
        unified_df <- Reduce(function(df1, df2) {
            merge(df1, df2, by = "NodeNames", all = TRUE, suffixes = c("", ""))
        }, combined_profiles)
        unified_df$MedianScore <- apply(unified_df[, -1], 1, median, na.rm = TRUE)
        unified_profile <- unified_df[, c("NodeNames", "MedianScore")]
        matching_nodes <- unified_profile$NodeNames %in% colnames(drug_gene_matrix)
        filtered_profile <- unified_profile[matching_nodes, ]
        if (nrow(filtered_profile) > 0) {
            drug_gene_matrix[drug, filtered_profile$NodeNames] <- filtered_profile$MedianScore
        } else {
            print(paste("No se encontraron nodos coincidentes para el fármaco:", drug))
        }
    }
}

head(drug_gene_matrix)
save(drug_gene_matrix, file = "results/Drug_profiles_PPI.RData")

### Disease-gene associations
dis_gene_matrix <- matrix(NA,
    nrow = length(diseases), ncol = length(genes),
    dimnames = list(disease_names, genes)
)

for (disease in disease_names) {
    print(disease)
    related_genes <- disgen %>%
        filter(DiseaseID == disease) %>%
        pull(GeneSymbol)

    difussion_profiles <- list()
    for (gene in related_genes) {
        if (gene %in% genes) {
            difussion_profiles[[gene]] <- profiles_list[[gene]]
        }
    }

    if (length(difussion_profiles) == 1) {
        print("DBG entra por 1")
        unified_profile <- difussion_profiles[[1]]
        matching_nodes <- unified_profile$NodeNames %in% colnames(dis_gene_matrix)
        filtered_profile <- unified_profile[matching_nodes, ]
        if (nrow(filtered_profile) > 0) {
            dis_gene_matrix[disease, filtered_profile$NodeNames] <- filtered_profile$Score
        } else {
            print(paste("No se encontraron nodos coincidentes para la enfermedad:", disease))
        }
    }

    if (length(difussion_profiles) > 1) {
        print("DBG entra por >1")
        combined_profiles <- list()
        for (profile in difussion_profiles) {
            combined_profiles[[length(combined_profiles) + 1]] <- profile
        }
        unified_df <- Reduce(function(df1, df2) {
            merge(df1, df2, by = "NodeNames", all = TRUE, suffixes = c("", ""))
        }, combined_profiles)
        unified_df$MedianScore <- apply(unified_df[, -1], 1, median, na.rm = TRUE)
        unified_profile <- unified_df[, c("NodeNames", "MedianScore")]
        matching_nodes <- unified_profile$NodeNames %in% colnames(dis_gene_matrix)
        filtered_profile <- unified_profile[matching_nodes, ]
        if (nrow(filtered_profile) > 0) {
            dis_gene_matrix[disease, filtered_profile$NodeNames] <- filtered_profile$MedianScore
        } else {
            print(paste("No se encontraron nodos coincidentes para la enfermedad:", disease))
        }
    }
}

head(dis_gene_matrix)
save(dis_gene_matrix, file = "results/Disease_profiles_PPI.RData")