library(edgeR)
library(GSVA)

source("functions.R")

load("cluster_assignments.Rdata")

#Keep only paired samples and remove lowly expressed genes (cpm<1 in all tumour or normal samples)
expression <- list()
for(tumour in tumours_normal){
  expression[[tumour]] <- read_expression_paired(tumour)
  print(tumour)
}

##Calculate ssGSEA scores of tumour modules in normal samples to compare how these genes
#have changed expression
ssGSEA_normal_exp_tumour_modules <- list()
for(tumour in tumours_normal){
  temp_normal <- expression[[tumour]][["normal"]]
  temp_tumour <- expression[[tumour]][["tumour"]]
  all_genes <- unique(c(rownames(temp_tumour), rownames(temp_normal)))
  colnames(temp_normal) <- paste("N", colnames(temp_normal), sep="_")
  colnames(temp_tumour) <- paste("T", colnames(temp_tumour), sep="_")
  
  temp_all <- data.frame(Genes=all_genes,
                         temp_tumour[match(all_genes, rownames(temp_tumour)),],
                         temp_normal[match(all_genes, rownames(temp_normal)),])
  rownames(temp_all) <- temp_all$Genes
  temp_all$Genes <- NULL
  
  tumour_samples <- colnames(temp_tumour)
  normal_samples <- colnames(temp_normal)
  
  temp_ssGSEA <- gsva(as.matrix(temp_all), cluster_assignments[["tumour"]][[tumour]], method="ssgsea", kcdf="Gaussian", ssgsea.norm=TRUE)
  temp_ssGSEA_tumour <- temp_ssGSEA[,colnames(temp_ssGSEA) %in% tumour_samples]
  temp_ssGSEA_normal <- temp_ssGSEA[,colnames(temp_ssGSEA) %in% normal_samples]
  
  ssGSEA_normal_exp_tumour_modules[[tumour]][["normal"]] <- temp_ssGSEA_normal
  ssGSEA_normal_exp_tumour_modules[[tumour]][["tumour"]] <- temp_ssGSEA_tumour
  print(tumour)
}

#save(ssGSEA_normal_exp_tumour_modules, file="ssGSEA_normal_exp_tumour_modules.Rdata")
