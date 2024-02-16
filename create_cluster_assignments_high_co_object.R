##Keep only the strongest links
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

tumours_other <- tumours[!(tumours %in% tumours_normal)]

subnet_normal <- list()
subnet_tumour <- list()
cluster_assignments_high <- list()
for(tumour in tumours_other){
  if(tumour %in% tumours_normal){
    load(paste("Subnetworks_", tumour, "_normal.Rdata", sep=""))
    subnet_normal[[tumour]] <- sub_networks$normal[[tumour]]
    
    for(module in names(subnet_normal[[tumour]])){
      temp <- subnet_normal[[tumour]][[module]]
      temp[temp == 1] <- NA
      local_cutoff <- median(temp, na.rm=TRUE)
      temp[temp < local_cutoff] <- NA
      n_connections <- apply(temp, 1, function(x){ return(sum(!is.na(x)))})
      #at least 50% of connections
      n_connections <- n_connections[n_connections >= nrow(subnet_normal[[tumour]][[module]])/2]
      cluster_assignments_high[["normal"]][[tumour]][[module]] <- names(n_connections)
    }
  }
  
  load(paste("Subnetworks_", tumour, "_tumour.Rdata", sep=""))
  subnet_tumour[[tumour]] <- sub_networks$tumour[[tumour]]
  
  for(module in names(subnet_tumour[[tumour]])){
    temp <- subnet_tumour[[tumour]][[module]]
    temp[temp == 1] <- NA
    local_cutoff <- median(temp, na.rm=TRUE)
    temp[temp < local_cutoff] <- NA
    n_connections <- apply(temp, 1, function(x){ return(sum(!is.na(x)))})
    #at least 50% of connections
    n_connections <- n_connections[n_connections >= nrow(subnet_tumour[[tumour]][[module]])/2]
    cluster_assignments_high[["tumour"]][[tumour]][[module]] <- names(n_connections)
  }
  print(tumour)
}


save(cluster_assignments_high, file="Objects/cluster_assignments_high.Rdata")
