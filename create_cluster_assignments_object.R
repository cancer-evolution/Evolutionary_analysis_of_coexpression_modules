library(WGCNA)
library(edgeR)
library(ggplot2)
library(reshape2)
source("functions.R")

##read in data
tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]


load("expression_WGCNA_paired.Rdata")

###Module detection

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft <- list()
for(tissue_type in c("tumour")){
  for(tumour in tumours){
    sft[[tissue_type]][[tumour]] = pickSoftThreshold(expression_WGCNA[[tumour]][[tissue_type]], powerVector = powers, verbose = 5)
  }
}

##Samples with normal samples
tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

for(tissue_type in c("normal")){
  for(tumour in tumours_normal){
    sft[[tissue_type]][[tumour]] = pickSoftThreshold(expression_WGCNA[[tumour]][[tissue_type]], powerVector = powers, verbose = 5)
  }
}

#save(sft, file="sft_paired.Rdata")

load("sft_paired.Rdata")
pdf("Choosing_soft_thresholds_paired.pdf")
par(mfrow=c(2,2))
for(tissue_type in c("tumour")){
  for(tumour in tumours){
    local_sft <- sft[[tissue_type]][[tumour]]
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence", tumour, tissue_type, sep="\n"));
    text(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         labels=powers,col="red");
    abline(h=0.90,col="red")
    
    # Mean connectivity as a function of the soft-thresholding power
    plot(local_sft$fitIndices[,1], local_sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity", tumour, tissue_type, sep="\n"))
    text(local_sft$fitIndices[,1], local_sft$fitIndices[,5], labels=powers, col="red")
  }
}
for(tissue_type in c("normal")){
  for(tumour in tumours_normal){
    local_sft <- sft[[tissue_type]][[tumour]]
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence", tumour, tissue_type, sep="\n"));
    text(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         labels=powers,col="red");
    abline(h=0.90,col="red")
    
    # Mean connectivity as a function of the soft-thresholding power
    plot(local_sft$fitIndices[,1], local_sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity", tumour, tissue_type, sep="\n"))
    text(local_sft$fitIndices[,1], local_sft$fitIndices[,5], labels=powers, col="red")
  }
}
dev.off()

##Choosing thresholds for each tumours type (based on the results with all genes)
soft_threshold <- list()
soft_threshold[["normal"]][["BLCA"]] <- 8
soft_threshold[["normal"]][["BRCA"]] <- 7
soft_threshold[["normal"]][["COAD"]] <- 9
soft_threshold[["normal"]][["ESCA"]] <- 16
soft_threshold[["normal"]][["HNSC"]] <- 9
soft_threshold[["normal"]][["KICH"]] <- 10
soft_threshold[["normal"]][["KIRC"]] <- 19
soft_threshold[["normal"]][["KIRP"]] <- 12
soft_threshold[["normal"]][["LIHC"]] <- 5
soft_threshold[["normal"]][["LUAD"]] <- 4
soft_threshold[["normal"]][["LUSC"]] <- 5
soft_threshold[["normal"]][["PRAD"]] <- 14
soft_threshold[["normal"]][["READ"]] <- 20
soft_threshold[["normal"]][["STAD"]] <- 9
soft_threshold[["normal"]][["THCA"]] <- 5
soft_threshold[["normal"]][["UCEC"]] <- 16


soft_threshold[["tumour"]][["ACC"]] <- 7
soft_threshold[["tumour"]][["BLCA"]] <- 6
soft_threshold[["tumour"]][["BRCA"]] <- 6
soft_threshold[["tumour"]][["CESC"]] <- 4
soft_threshold[["tumour"]][["CHOL"]] <- 5
soft_threshold[["tumour"]][["COAD"]] <- 6
soft_threshold[["tumour"]][["ESCA"]] <- 12
soft_threshold[["tumour"]][["GBM"]] <- 4
soft_threshold[["tumour"]][["HNSC"]] <- 6
soft_threshold[["tumour"]][["KICH"]] <- 10
soft_threshold[["tumour"]][["KIRC"]] <- 8
soft_threshold[["tumour"]][["KIRP"]] <- 9
soft_threshold[["tumour"]][["LGG"]] <- 5
soft_threshold[["tumour"]][["LIHC"]] <- 4
soft_threshold[["tumour"]][["LUAD"]] <- 7
soft_threshold[["tumour"]][["LUSC"]] <- 5
soft_threshold[["tumour"]][["MESO"]] <- 4
soft_threshold[["tumour"]][["OV"]] <- 4
soft_threshold[["tumour"]][["PAAD"]] <- 7
soft_threshold[["tumour"]][["PCPG"]] <- 7
soft_threshold[["tumour"]][["PRAD"]] <- 7
soft_threshold[["tumour"]][["READ"]] <- 20
soft_threshold[["tumour"]][["SARC"]] <- 5
soft_threshold[["tumour"]][["SKCM"]] <- 12
soft_threshold[["tumour"]][["STAD"]] <- 8
soft_threshold[["tumour"]][["TGCT"]] <- 12
soft_threshold[["tumour"]][["THCA"]] <- 9
soft_threshold[["tumour"]][["THYM"]] <- 6
soft_threshold[["tumour"]][["UCEC"]] <- 18
soft_threshold[["tumour"]][["UCS"]] <- 4
soft_threshold[["tumour"]][["UVM"]] <- 14

save(soft_threshold, file="soft_threshold.Rdata")

#Adjacency matrix and TOM
minModuleSize = 30;

cluster_assignments <- list()

for(tissue_type in c("tumour")){
  for(tumour in tumours){
    adjacency_list <- adjacency(expression_WGCNA[[tumour]][[tissue_type]], power = soft_threshold[[tissue_type]][[tumour]])
    TOM = TOMsimilarity(adjacency_list)
    dissTOM = 1-TOM
    
    rownames(TOM) <- rownames(adjacency_list)
    colnames(TOM) <- colnames(adjacency_list)
    geneTree = hclust(as.dist(dissTOM), method = "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, 
                                distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    dynamicColors = labels2colors(dynamicMods)
    
     # Calculate eigengenes
    MEList = moduleEigengenes(expression_WGCNA[[tumour]][[tissue_type]], colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");

    MEDissThres = 0.25  ##Merge modules with a correlation of at least 0.75

    # Call an automatic merging function
    merge = mergeCloseModules(expression_WGCNA[[tumour]][[tissue_type]], dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs;
    cluster_names <- unique(mergedColors)
    
    for(cluster in cluster_names){
      indices <- mergedColors == cluster
      cluster_assignments[[tissue_type]][[tumour]][[cluster]] <- colnames(expression_WGCNA[[tumour]][[tissue_type]])[indices]
    }
    print(tumour)
  }
}

for(tissue_type in c("normal")){
  for(tumour in tumours_normal){
    adjacency_list <- adjacency(expression_WGCNA[[tumour]][[tissue_type]], power = soft_threshold[[tissue_type]][[tumour]])
    TOM = TOMsimilarity(adjacency_list)
    dissTOM = 1-TOM
    
    rownames(TOM) <- rownames(adjacency_list)
    colnames(TOM) <- colnames(adjacency_list)
    geneTree = hclust(as.dist(dissTOM), method = "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, 
                                distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    dynamicColors = labels2colors(dynamicMods)
    

    # Calculate eigengenes
    MEList = moduleEigengenes(expression_WGCNA[[tumour]][[tissue_type]], colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");
    
    MEDissThres = 0.25  ##Merge modules with a correlation of at least 0.75

    # Call an automatic merging function
    merge = mergeCloseModules(expression_WGCNA[[tumour]][[tissue_type]], dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs;
    
    cluster_names <- unique(mergedColors)
    
    for(cluster in cluster_names){
      indices <- mergedColors == cluster
      cluster_assignments[[tissue_type]][[tumour]][[cluster]] <- colnames(expression_WGCNA[[tumour]][[tissue_type]])[indices]
    }
    print(tumour)
  }
}

#save(cluster_assignments, file="cluster_assignments.Rdata")
