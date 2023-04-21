library(reshape2)
library(ggraph)
library(igraph)

source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

subnet_normal <- list()
subnet_tumour <- list()
for(tumour in tumours){
  if(tumour %in% tumours_normal){
    load(paste("Subnetworks_", tumour, "_normal.Rdata", sep=""))
    subnet_normal[[tumour]] <- sub_networks$normal[[tumour]]
  }
  
  load(paste("Subnetworks_", tumour, "_tumour.Rdata", sep=""))
  subnet_tumour[[tumour]] <- sub_networks$tumour[[tumour]]
}

##Calculate degree of tumour and normal modules. Normal ones seem to have stronger
#hubs than tumour ones

degree_modules <- vector()
for(tumour in tumours){
  degree_modules_tumour <- calculate_degree_in_modules(subnet_tumour[[tumour]], tumour, "Tumour")
  degree_modules_normal <- calculate_degree_in_modules(subnet_normal[[tumour]], tumour, "Normal")
  
  degree_modules <- rbind(degree_modules,
                          degree_modules_tumour,
                          degree_modules_normal)
  
}
degree_modules$Tissue_type <- factor(degree_modules$Tissue_type, levels=c("Normal", "Tumour"))
#save(degree_modules, file="degree_modules.Rdata")