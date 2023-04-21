library(WGCNA)
library(edgeR)
library(ggplot2)
library(reshape2)
source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

load("soft_threshold.Rdata")
load("cluster_assignments.Rdata")
load("expression_WGCNA_paired.Rdata")


#Adjacency matrix and TOM

for(tissue_type in c("normal", "tumour")){
  for(tumour in tumours){
    if(tissue_type == "tumour" || (tumour %in% tumours_normal & tissue_type == "normal")){
      sub_networks <- list()
      adjacency_list <- adjacency(expression_WGCNA[[tumour]][[tissue_type]], power = soft_threshold[[tissue_type]][[tumour]])
      TOM = TOMsimilarity(adjacency_list)
      
      rownames(TOM) <- rownames(adjacency_list)
      colnames(TOM) <- colnames(adjacency_list)
      
      local_modules <- cluster_assignments[[tissue_type]][[tumour]]
      for(mod in names(local_modules)){
        if(mod != "grey"){
          genes_in_mod <- local_modules[[mod]]
          sub_networks[[tissue_type]][[tumour]][[mod]] <- TOM[genes_in_mod, genes_in_mod]
        }
      }
      save(sub_networks, file=paste("Subnetworks_", tumour, "_", tissue_type, ".Rdata", sep=""))
      
    }
  }
}

