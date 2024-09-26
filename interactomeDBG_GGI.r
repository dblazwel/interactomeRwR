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

set.seed(444)
print("GGI Network")
data("go")
go_terms <- go$id
get_go_ontology <- function(go_term) {
  term <- GOTERM[[go_term]]
  if (!is.null(term)) {
    return(Ontology(term))
  } else {
    return(NA)
  }
}
go_ontology <- sapply(go_terms, get_go_ontology)
go_data <- data.frame(go_id = go_terms, ontology = go_ontology)
go_bp <- go_data %>% filter(ontology == "BP") # Nos quedamos solo con BP
go_id_list <- as.list(go_bp$go_id)
print(paste0("Num de terminos total: ", length(go_terms), " / Num de terminos BP: ", length(go_id_list)))
save(go_id_list, file = "GO_id_list")
print("###############################################")
load("data/GO_id_list")
all_similarity_dfs <- list()
for (i in 1:length(go_id_list)) {
  # Calcular la matriz de similitud
  similarity_matrix <- get_sim_grid(
    ontology = go,
    term_sim_method = "resnik",
    term_sets = go_id_list[i],
    term_sets2 = go_id_list,
    combine = "average"
  )
  # Convertir la matriz a un formato data.frame y almacenarlo en la lista
  similarity_df <- as.data.frame(as.matrix(similarity_matrix))
  all_similarity_dfs[[i]] <- similarity_df
}
# Combinar todos los data.frames en uno solo
combined_similarity_df <- do.call(rbind, all_similarity_dfs)
rownames(combined_similarity_df) <- go_id_list # [25234:27597]
colnames(combined_similarity_df) <- go_id_list
head(combined_similarity_df)
# write.csv(x = combined_similarity_df, file = "similarity_matrix_FINAL.csv", row.names = TRUE, eol = "\r\n")
# sim_data <- fread("similarity_matrix_FINAL.csv")
sim_data <- combined_similarity_df
matriz_similitud <- as.matrix(sim_data[, -1])
rownames(matriz_similitud) <- sim_data$V1
print("Numero de GGIs inicial: ")
dim(matriz_similitud)
valores_superior <- upper.tri(matriz_similitud, diag = FALSE) # Nos quedamos solo con los valores de la diagonal superior
valores_superior <- valores_superior * matriz_similitud
Q99 <- quantile(valores_superior, 0.99)
Q99
gene_relation_Matrix <- as.matrix(valores_superior > Q99)
gene_relation_Matrix <- gene_relation_Matrix * valores_superior

matrix_melt <- melt(gene_relation_Matrix, varnames = c("GENE1", "GENE2"), value.name = "IC")
matrix_melt <- matrix_melt[matrix_melt$IC != 0, ] # Filter out those GGIs with no IC

print("Numero de GGIs final: ")
dim(matrix_melt)
write_csv(matrix_melt, "matrixMelted.csv")
head(matrix_melt)

CairoPNG("distribucion_ic.png", width = 800, height = 600)
hist(matrix_melt$IC, main = "Distribución de IC", xlab = "IC", ylab = "Frecuencia", col = "orange", border = "black")
dev.off()

# matrix_melt <- fread("data/matrixMelted.csv")

# Transformación ENSG --> GeneSymbol
# hgu95av2.db --> ENSEMBL / SYMBOL
keys1 <- as.character(matrix_melt$GENE1)
symbols1 <- mapIds(hgu95av2.db, keys = keys1, column = "SYMBOL", keytype = "GO", multiVals = "list")
matrix_melt$SYMBOL1 <- symbols1[as.character(matrix_melt$GENE1)]

keys2 <- as.character(matrix_melt$GENE2)
symbols2 <- mapIds(hgu95av2.db, keys = keys2, column = "SYMBOL", keytype = "GO", multiVals = "list")
matrix_melt$SYMBOL2 <- symbols2[as.character(matrix_melt$GENE2)]

# Filter out NAs & Combine all possible symbols 
matrix_melt_filt <- matrix_melt[!is.na(SYMBOL1) & !is.na(SYMBOL2),]

matrix_melt_FINAL <- as.data.table(matrix_melt_filt[1:100,])
matrix_melt_FINAL <- matrix_melt_FINAL[, {
  # Explosión de listas de SYMBOL1 y SYMBOL2 en combinaciones
  comb <- CJ(unlist(SYMBOL1), unlist(SYMBOL2), unique = TRUE, sorted = FALSE)
  comb[, IC := IC] 
  comb
}, by = .(GENE1, GENE2)]
matrix_melt_FINAL <- matrix_melt_FINAL[, c("V1", "V2", "IC")]
colnames(matrix_melt_FINAL) <- c("GENE1", "GENE2", "IC")

#matrix_melt_filt <- matrix_melt[complete.cases(matrix_melt[, c("SYMBOL1", "SYMBOL2")]), c("SYMBOL1", "SYMBOL2", "IC")]
print("Nº GGIs mapeadas a Gene Symbol: ")
dim(matrix_melt_FINAL)
head(matrix_melt_FINAL)
write_csv(matrix_melt_FINAL,"data/genes_matrix.csv")
gene_gene_graph <- graph_from_data_frame(matrix_melt_FINAL, directed = FALSE)
GGI_MultiplexObject <- create.multiplex(list(GGI=gene_gene_graph))
# AdjMatrix_GGI <- compute.adjacency.matrix(GGI_MultiplexObject)
GGI_MultiplexObject
save(GGI_MultiplexObject, file = "data/networks/GGI")