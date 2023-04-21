library(WGCNA)
library(ggplot2)

expression_pheo <- read.delim("merged_expression_admixture_removed.tsv")
genes <- expression_pheo$X

annotations_pheo <- read.delim("merged_expression_admixture_removed_metadata.tsv")
annotations_pheo$Malignancy <- as.character(annotations_pheo$Malignancy)
annotations_pheo$Malignancy[annotations_pheo$Genotype_Cluster == "Normal"] <- "Normal"

annotations_pheo_SDH <- annotations_pheo[annotations_pheo$Genotype_Cluster %in% c("SDHx", "SDHx (H&N)"),]
table(annotations_pheo_SDH$Malignancy)

normal <- as.character(annotations_pheo$Sample[annotations_pheo$Malignancy == "Normal"])
benign <- as.character(annotations_pheo_SDH$Sample[annotations_pheo_SDH$Malignancy == "Benign"])
malignant <- as.character(annotations_pheo_SDH$Sample[annotations_pheo_SDH$Malignancy == "Malignant"])
normal <- normal[!is.na(normal)]
benign <- benign[!is.na(benign)]
malignant <- malignant[!is.na(malignant)]

expression <- list()
expression[["normal"]] <- expression_pheo[,colnames(expression_pheo) %in% normal]
rownames(expression[["normal"]]) <- genes
expression[["benign"]] <- expression_pheo[,colnames(expression_pheo) %in% benign]
rownames(expression[["benign"]]) <- genes
expression[["malignant"]] <- expression_pheo[,colnames(expression_pheo) %in% malignant]
rownames(expression[["malignant"]]) <- genes

expression[["normal"]] <- t(expression[["normal"]])
expression[["benign"]] <- t(expression[["benign"]])
expression[["malignant"]] <- t(expression[["malignant"]])

for(type in c("normal", "benign", "malignant")){
  gsg = goodSamplesGenes(expression[[type]], verbose = 3);
  gsg$allOK
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(expression[[type]])[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(expression[[type]])[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    expression[[type]] = expression[[type]][gsg$goodSamples, gsg$goodGenes]
  }
}

for(type in c("normal", "benign", "malignant")){
  sampleTree = hclust(dist(expression[[type]]), method = "average");
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = paste("Sample clustering to detect outliers", type, sep="\n"), sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)  
}


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft <- list()
for(type in c("normal", "benign", "malignant")){
    sft[[type]] = pickSoftThreshold(expression[[type]], powerVector = powers, verbose = 5,
                                    dataIsExpr = T)
}

#save(sft, file="sft_paired_pheo_SDH.Rdata")
load("sft_paired_pheo_SDH.Rdata")


par(mfrow=c(2,2))
for(type in c("normal", "benign", "malignant")){
    local_sft <- sft[[type]]
    # Scale-free topology fit index as a function of the soft-thresholding power
  plot(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence", type, sep="\n"));
  text(local_sft$fitIndices[,1], -sign(local_sft$fitIndices[,3])*local_sft$fitIndices[,2],
         labels=powers,col="red");
  abline(h=0.90,col="red")
    
    # Mean connectivity as a function of the soft-thresholding power
  plot(local_sft$fitIndices[,1], local_sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity", type, sep="\n"))
  text(local_sft$fitIndices[,1], local_sft$fitIndices[,5], labels=powers, col="red")
}

##Choosing thresholds for each tumours type (based on the results with all genes)
soft_threshold <- list()
soft_threshold[["normal"]] <- 5
soft_threshold[["benign"]] <- 7
soft_threshold[["malignant"]] <- 8


#Adjacency matrix and TOM
minModuleSize = 30;

modules_pheo <- list()

for(type in c("normal", "benign", "malignant")){
    adjacency_list <- adjacency(expression[[type]], power = soft_threshold[[type]])
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

    
    cluster_names <- unique(dynamicColors)
    sub_networks <- list()
    for(cluster in cluster_names){
      indices <- dynamicColors == cluster
      modules_pheo[[type]][[cluster]] <- colnames(expression[[type]])[indices]
      local_modules <- modules_pheo[[type]]
      for(mod in names(local_modules)){
        if(mod != "grey"){
          genes_in_mod <- local_modules[[mod]]
          sub_networks[[type]][[mod]] <- TOM[genes_in_mod, genes_in_mod]
        }
      }
      save(sub_networks, file=paste("Subnetworks_pheo_", type, "_SDH.Rdata", sep=""))
      
    }
    print(type)
}

save(modules_pheo, file="modules_pheo_SDH.Rdata")


##Calculate age and novelty, and centrality of genes

library(ggplot2)
library(gridExtra)
library(Biobase)
library(GO.db)
library(igraph)
library(reshape2)

source("functions.R")

load("modules_pheo_SDH.Rdata")

gene_ages <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_in_clusters <- unique(unname(unlist(modules_pheo)))

gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 1:3),1]
MC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 4:16),1]


percentage_UC_cluster <- vector()
for(type in c("normal", "benign", "malignant")){
  all_genes <- unique(unname(unlist(modules_pheo[[type]])))
  local_exp_per_UC <- sum(all_genes %in% UC_genes)/c(sum(all_genes %in% UC_genes) + sum(all_genes %in% MC_genes))*100
  cluster_names <- names(modules_pheo[[type]])
  for(cluster_name in cluster_names){
    genes_in_cluster <- modules_pheo[[type]][[cluster_name]]
    UC_genes_cluster <- sum(genes_in_cluster %in% UC_genes)
    MC_genes_cluster <- sum(genes_in_cluster %in% MC_genes)
    p_UC <- UC_genes_cluster/(UC_genes_cluster + MC_genes_cluster)*100
    percentage_UC_cluster <- rbind(percentage_UC_cluster,
                                     c(type, cluster_name, length(genes_in_cluster), p_UC,
                                       p_UC-local_exp_per_UC))
  }
}

percentage_UC_cluster <- as.data.frame(percentage_UC_cluster)
colnames(percentage_UC_cluster) <- c("type", "cluster_name", "number_genes_cluster", "percentage_UC", "diff_p_UC")
percentage_UC_cluster$percentage_UC <- as.numeric(as.character(percentage_UC_cluster$percentage_UC))
percentage_UC_cluster$diff_p_UC <- as.numeric(as.character(percentage_UC_cluster$diff_p_UC))

#Enrichment tests (taking into account that each sample type has a diff number of total UC/MC genes)
age_enrichment <- list()

for(type in c("normal", "benign", "malignant")){
  cluster_names <- names(modules_pheo[[type]])
  all_genes <- unique(unname(unlist(modules_pheo[[type]])))
  UC_genes_total <- all_genes[all_genes %in% UC_genes]
  MC_genes_total <- all_genes[all_genes %in% MC_genes]
      
  for(cluster_name in cluster_names){
      genes_in_cluster <- modules_pheo[[type]][[cluster_name]]
      UC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% UC_genes]
        
      MC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% MC_genes]
        
      df <- cbind(c(length(UC_genes_cluster), length(MC_genes_cluster)),
                    c(length(UC_genes_total), length(MC_genes_total)))
      UC_p <- fisher.test(df, alternative="greater")$p.value
      MC_p <- fisher.test(df, alternative="less")$p.value
        
      age_enrichment[[type]] <- rbind(age_enrichment[[type]],
                                     c(type, cluster_name, UC_p, MC_p))
  }
  colnames(age_enrichment[[type]]) <- c("type", "cluster", "p_UC_enrichment", "p_MC_enrichment")
  age_enrichment[[type]] <- as.data.frame(age_enrichment[[type]])
  age_enrichment[[type]]$p_UC_enrichment <- as.numeric(as.character(age_enrichment[[type]]$p_UC_enrichment))
  age_enrichment[[type]]$p_MC_enrichment <- as.numeric(as.character(age_enrichment[[type]]$p_MC_enrichment))
  age_enrichment[[type]]$p_UC_adj <- p.adjust(age_enrichment[[type]]$p_UC_enrichment, method="BH")
  age_enrichment[[type]]$p_MC_adj <- p.adjust(age_enrichment[[type]]$p_MC_enrichment, method="BH")
  age_enrichment[[type]]$Module_age <- ifelse(age_enrichment[[type]]$p_UC_adj < 0.05, "UC", 
                                               ifelse(age_enrichment[[type]]$p_MC_adj < 0.05, "MC",
                                                      "Mixed"))
}

save(age_enrichment, file="age_enrichment_pheo_SDH.Rdata")
#load("age_enrichment_pheo_SDH.Rdata")

table(age_enrichment$normal$Module_age)

table(age_enrichment$benign$Module_age)

table(age_enrichment$malignant$Module_age)

number_of_shared_genes <- vector()
for(type1 in c("normal", "benign", "malignant")){
  for(type2 in c("normal", "benign", "malignant")){
    if(type1 != type2){
      clusters_of_local_sample_type <- modules_pheo[[type1]]
      clusters_of_other_sample_type <- modules_pheo[[type2]]
      
      for(local_name1 in names(clusters_of_local_sample_type)){
        genes_in_local_cluster <- clusters_of_local_sample_type[[local_name1]]
        for(local_name2 in names(clusters_of_other_sample_type)){
          genes_in_other_cluster <- clusters_of_other_sample_type[[local_name2]]
          shared_genes <- sum(genes_in_local_cluster %in% genes_in_other_cluster)/length(genes_in_local_cluster)*100
          number_of_shared_genes <- rbind(number_of_shared_genes,
                                          c(local_name1, local_name2, type1, 
                                            type2, shared_genes))
          print(local_name2)
        }
      }
    }
  }
}

number_of_shared_genes <- as.data.frame(number_of_shared_genes)
colnames(number_of_shared_genes) <- c("Module1", "Module2", "Type1", "Type2", "Shared_genes")
number_of_shared_genes$Shared_genes <- as.numeric(as.character(number_of_shared_genes$Shared_genes))

#save(number_of_shared_genes, file="number_of_shared_genes_pheo_SDH.Rdata")

##Novelty
source("functions.R")
load("number_of_shared_genes_pheo_SDH.Rdata")
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module1 != "grey",]
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module2 != "grey",]

preservation_pheo <- calculate_preservation_pheo(number_of_shared_genes)

#save(preservation_pheo, file="preservation_pheo_SDH.Rdata")
load("preservation_pheo_SDH.Rdata")

cluster_sizes <- vector()
for(type in c("normal", "benign", "malignant")){
  temp <- modules_pheo[[type]]
  mod <- sapply(temp, length)
  mod <- cbind(Module=names(mod), Size=as.vector(mod), Type=type)
  cluster_sizes <- rbind(cluster_sizes, mod)
}

cluster_sizes <- as.data.frame(cluster_sizes)
cluster_sizes$Size <- as.numeric(as.character(cluster_sizes$Size))
cluster_sizes$Label <- paste(cluster_sizes$Type, cluster_sizes$Module)

preservation_pheo$Label <- paste(preservation_pheo$Type, preservation_pheo$Cluster)

preservation_pheo$Size <- cluster_sizes[match(preservation_pheo$Label, cluster_sizes$Label),"Size"]

preservation_pheo$Per_50 <- as.numeric(as.character(preservation_pheo$Per_50))

save(preservation_pheo, file="preservation_pheo_SDH.Rdata")



load("age_enrichment_pheo_SDH.Rdata")



all_preservation_pheo <- preservation_pheo[,c("Type", "Cluster", "Per_50", "Label", "Size")]

#Divide preservation score by module size?
all_preservation_pheo$Preservation_ratio <- all_preservation_pheo$Per_50/all_preservation_pheo$Size


##Age
age_enrichment_df <- vector()
for(type in c("normal", "benign", "malignant")){
  local_age_enrichment <- age_enrichment[[type]]
  age_enrichment_df <- rbind(age_enrichment_df, local_age_enrichment)
}
age_enrichment_df$Label <- paste(age_enrichment_df$type, age_enrichment_df$cluster)

all_preservation_pheo$Age <- age_enrichment_df$Module_age[match(all_preservation_pheo$Label, 
                                                                age_enrichment_df$Label)]
all_preservation_pheo$Age <- factor(all_preservation_pheo$Age, levels=c("UC", "Mixed", "MC"))

all_preservation_pheo <- define_categories_preservation_pheo(all_preservation_pheo)

all_preservation_pheo$Age <- factor(all_preservation_pheo$Age,
                                    levels=c("UC", "Mixed", "MC"))

save(all_preservation_pheo, file="all_preservation_pheo_SDH.Rdata")
