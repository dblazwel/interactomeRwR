library(dplyr)
library(readr)
library(stringr)
library(reshape)
library(reshape2)
library(tidyverse)
library(RandomWalkRestartMH)
library(igraph)
library(ggplot2)
library(ontologyIndex)
#library(ontologySimilarity)
library(GO.db)

print('###############################################')
print('### DBG interactome ####')
print('###############################################')

set.seed(444)

## DRUG/DISEASE
ICD10_disease_drug <- read_csv("data/ICD10_disease_drug.csv")
ICD10_dictionary <- read_csv("data/ICD10_dictionary.csv")
CTD_gene_disease <- read.table("data/CTD_genes_diseases_direct.tsv", 
                 header = TRUE,   # La primera línea contiene los nombres de las columnas
                 sep = "\t",      # Separador de columnas es tabulación
                 quote = "",      # No hay comillas que rodeen los campos
                 comment.char = "")# No hay caracteres de comentario
print('ICD10_disease_drug')
head(ICD10_disease_drug)
print('###############################################')
print('ICD10_dictionary')
head(ICD10_dictionary)
print('###############################################')
print('CTD_gene_disease')
head(CTD_gene_disease)
print('###############################################')

gene_list_original <- levels(matrix_melt$GENE1)
gene_list_CTD <- CTD_gene_disease$GeneSymbol
disease_list_CTD <- CTD_gene_disease$DiseaseName

filtered_diseases <- ICD10_dictionary %>%  filter(str_detect(tolower(disease), paste(tolower(disease_list_CTD), collapse = "|")))
filtered_diseases
write_csv(filtered_diseases, "ICD10_dictionary_filtered.csv")
print(paste0("Enfermedades en el diccionario: ", dim(ICD10_dictionary), " / Enfermedades mapeadas: ", dim(filtered_diseases)))
print('###############################################')
matched_rows <- list()
# Recorrer cada enfermedad filtrada y encontrar las filas correspondientes en CTD_gene_disease
for (i in seq_along(filtered_diseases$disease)) {
  current_disease <- filtered_diseases$disease[i]
  ICD10_index <- filtered_diseases$icd10[i]
  # Buscar las filas correspondientes en CTD_gene_disease
  matches <- CTD_gene_disease %>% filter(str_detect(tolower(DiseaseName), tolower(current_disease)))
  if (nrow(matches) > 0) {
    matched_rows[[i]] <- data.frame(GeneSymbol = matches$GeneSymbol,
                                  GeneID = matches$GeneID,
                                  DiseaseName = matches$DiseaseName,
                                  GeneralTerm = current_disease,
                                  ICD10 = rep(ICD10_index, nrow(matches)))
  }
}
gene_disease_rel <- do.call(rbind, matched_rows)
write_csv(gene_disease_rel, "data/gene_disease_rel.csv")

