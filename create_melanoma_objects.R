##Data accessed from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98394

library(WGCNA)
library(ggplot2)

expression_melanoma <- read.delim("GSE98394_expression.txt")

samples <- colnames(expression_melanoma)


metadata <- read.delim("GSE98394_series_matrix.txt",
                       skip=51)
metadata <- metadata[-c(1,2,3,4,5,6,8,9,21:53),]
metadata <- t(metadata)
colnames(metadata) <- c("Source", "Tissue", "Tumour_thickness", "TNM_T", "TNM_N", "TNM_M", "Stage", "Ulceration", "Time_death_followup", "Survival_status", "Treatment", "Biopsy")
metadata <- metadata[-1,]
metadata <- as.data.frame(metadata)
metadata$Source <- NULL
metadata$Tissue <- gsub("tissue: ", "", metadata$Tissue)
metadata$Tumour_thickness <- as.character(metadata$Tumour_thickness)
metadata$Tumour_thickness <- substr(metadata$Tumour_thickness, 44, nchar(metadata$Tumour_thickness))
metadata$Tumour_thickness[metadata$Tumour_thickness == "n.a."] <- NA
metadata$TNM_T <- as.character(metadata$TNM_T)
metadata$TNM_T <- substr(metadata$TNM_T, 50, nchar(metadata$TNM_T))
metadata$TNM_T[metadata$TNM_T == "n.a."] <- NA

metadata$TNM_N <- as.character(metadata$TNM_N)
metadata$TNM_N <- substr(metadata$TNM_N, 51, nchar(metadata$TNM_N))
metadata$TNM_N[grep("n.a.", metadata$TNM_N)] <- NA

metadata$TNM_M <- as.character(metadata$TNM_M)
metadata$TNM_M <- substr(metadata$TNM_M, 55, nchar(metadata$TNM_M))
metadata$TNM_M[grep("n.a.", metadata$TNM_M)] <- NA

metadata$Stage <- as.character(metadata$Stage)
metadata$Stage <- substr(metadata$Stage, 36, nchar(metadata$Stage))
metadata$Stage[grep("n.a.", metadata$Stage)] <- NA

metadata$Ulceration <- as.character(metadata$Ulceration)
metadata$Ulceration <- substr(metadata$Ulceration, 13,  nchar(metadata$Ulceration))
metadata$Ulceration[grep("n.a.", metadata$Ulceration)] <- NA

metadata$Time_death_followup <- as.character(metadata$Time_death_followup)
metadata$Time_death_followup <- substr(metadata$Time_death_followup, 38,  nchar(metadata$Time_death_followup))
metadata$Time_death_followup[grep("n.a.", metadata$Time_death_followup)] <- NA

metadata$Survival_status <- as.character(metadata$Survival_status)
metadata$Survival_status <- substr(metadata$Survival_status, 36, nchar(metadata$Survival_status))
metadata$Survival_status[grep("n.a.", metadata$Survival_status)] <- NA

metadata$Treatment <- NULL
metadata$Biopsy <- NULL

rownames(metadata) <- substr(rownames(metadata), nchar(rownames(metadata))-5, nchar(rownames(metadata))-1)


sample_map <- data.frame(Samples = samples,
                         Type = "Primary")
nevi <- c("JC050","JC052","JC053","JC054","JC055","JC056","JC058","JC059","JC060","JC066","JC067",
          "JC069","JC072","JC077","JC247","JC248","JC262","JC271","JC272","JC273","JC274","JC281",
          "JC282","JC286","JC301","JC302","JC303")

sample_map$Type <- as.character(sample_map$Type)
sample_map$Type[sample_map$Samples %in% nevi] <- "Nevus"

sample_map$Thickness <- metadata$Tumour_thickness[match(sample_map$Samples, rownames(metadata))]
sample_map$TNM_T <- metadata$TNM_T[match(sample_map$Samples, rownames(metadata))]
sample_map$TNM_N <- metadata$TNM_N[match(sample_map$Samples, rownames(metadata))]
sample_map$Stage <- metadata$Stage[match(sample_map$Samples, rownames(metadata))]
sample_map$Ulceration <- metadata$Ulceration[match(sample_map$Samples, rownames(metadata))]
sample_map$Time_death_followup <- metadata$Time_death_followup[match(sample_map$Samples, rownames(metadata))]
sample_map$Survival_status <- metadata$Survival_status[match(sample_map$Samples, rownames(metadata))]


nevi <- sample_map$Samples[sample_map$Type == "Nevus"]
primary <- sample_map$Samples[sample_map$Type == "Primary"]
##27 Nevi, 51 primary



library(biomaRt)
mart <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
ensg <- getBM(mart = mart, attributes=c('ensembl_gene_id', 'external_gene_name'), filter='ensembl_gene_id',
              values=expression_melanoma$gene)
#save( ensg, file='/home/atrigos/Paper_4/Objects/ensg.Rdata')
load('/home/atrigos/Paper_4/Objects/ensg.Rdata')
expression_melanoma$GeneID <- ensg$external_gene_name[match(expression_melanoma$gene,
                                                            ensg$ensembl_gene_id)]
expression_melanoma <- expression_melanoma[!is.na(expression_melanoma$GeneID),]

rownames(expression_melanoma) <- expression_melanoma$GeneID

expression_melanoma$gene <- NULL
expression_melanoma$GeneID <- NULL

expression <- list()
expression[["nevus"]] <- expression_melanoma[,colnames(expression_melanoma) %in% nevi]
expression[["primary"]] <- expression_melanoma[,colnames(expression_melanoma) %in% primary]

expression[["nevus"]] <- t(expression[["nevus"]])
expression[["primary"]] <- t(expression[["primary"]])


for(type in c("nevus", "primary")){
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

for(type in c("nevus", "primary")){
  sampleTree = hclust(dist(expression[[type]]), method = "average");
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = paste("Sample clustering to detect outliers", type, sep="\n"), sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)  
}

primary_remove <- sampleTree$labels[sampleTree$height > 200]
expression$primary <- expression$primary[!(rownames(expression$primary) %in% primary_remove),]


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft <- list()
for(type in c("nevus", "primary")){
    sft[[type]] = pickSoftThreshold(expression[[type]], powerVector = powers, verbose = 5,
                                    dataIsExpr = T)
}

#save(sft, file="sft_paired_melanoma.Rdata")
load("sft_paired_melanoma.Rdata")


par(mfrow=c(2,2))
for(type in c("nevus", "primary")){
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
soft_threshold[["nevus"]] <- 8
soft_threshold[["primary"]] <- 10


#Adjacency matrix and TOM
minModuleSize = 30;

modules_melanoma <- list()


save(expression, file="expression_melanoma.Rdata")



for(type in c("nevus", "primary")){
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
      modules_melanoma[[type]][[cluster]] <- colnames(expression[[type]])[indices]
      local_modules <- modules_melanoma[[type]]
      for(mod in names(local_modules)){
        if(mod != "grey"){
          genes_in_mod <- local_modules[[mod]]
          sub_networks[[type]][[mod]] <- TOM[genes_in_mod, genes_in_mod]
        }
      }
      save(sub_networks, file=paste("Subnetworks_melanoma_", type, ".Rdata", sep=""))
      
    }
    print(type)
}

#save(modules_melanoma, file="modules_melanoma.Rdata")
load("modules_melanoma.Rdata")

##Calculate age and novelty, and centrality of genes
library(ggplot2)
library(gridExtra)
library(Biobase)
library(GO.db)
library(igraph)
library(reshape2)

source("functions.R")

load("modules_melanoma.Rdata")

gene_ages <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_in_clusters <- unique(unname(unlist(modules_melanoma)))

gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 1:3),1]
MC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 4:16),1]


percentage_UC_cluster <- vector()
for(type in c("nevus", "primary", "malignant")){
  all_genes <- unique(unname(unlist(modules_melanoma[[type]])))
  local_exp_per_UC <- sum(all_genes %in% UC_genes)/c(sum(all_genes %in% UC_genes) + sum(all_genes %in% MC_genes))*100
  cluster_names <- names(modules_melanoma[[type]])
  for(cluster_name in cluster_names){
    genes_in_cluster <- modules_melanoma[[type]][[cluster_name]]
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

for(type in c("nevus", "primary")){
  cluster_names <- names(modules_melanoma[[type]])
  all_genes <- unique(unname(unlist(modules_melanoma[[type]])))
  UC_genes_total <- all_genes[all_genes %in% UC_genes]
  MC_genes_total <- all_genes[all_genes %in% MC_genes]
      
  for(cluster_name in cluster_names){
      genes_in_cluster <- modules_melanoma[[type]][[cluster_name]]
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

#save(age_enrichment, file="age_enrichment_melanoma.Rdata")
load("age_enrichment_melanoma.Rdata")

table(age_enrichment$nevus$Module_age)

table(age_enrichment$primary$Module_age)


number_of_shared_genes <- vector()
for(type1 in c("nevus", "primary")){
  for(type2 in c("nevus", "primary")){
    if(type1 != type2){
      clusters_of_local_sample_type <- modules_melanoma[[type1]]
      clusters_of_other_sample_type <- modules_melanoma[[type2]]
      
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

#save(number_of_shared_genes, file="number_of_shared_genes_melanoma.Rdata")

##Novelty
source("functions.R")
load("number_of_shared_genes_melanoma.Rdata")
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module1 != "grey",]
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module2 != "grey",]

preservation_melanoma <- calculate_preservation_melanoma(number_of_shared_genes)

cluster_sizes <- vector()
for(type in c("nevus", "primary")){
  temp <- modules_melanoma[[type]]
  mod <- sapply(temp, length)
  mod <- cbind(Module=names(mod), Size=as.vector(mod), Type=type)
  cluster_sizes <- rbind(cluster_sizes, mod)
}

cluster_sizes <- as.data.frame(cluster_sizes)
cluster_sizes$Size <- as.numeric(as.character(cluster_sizes$Size))
cluster_sizes$Label <- paste(cluster_sizes$Type, cluster_sizes$Module)

preservation_melanoma$Label <- paste(preservation_melanoma$Type, preservation_melanoma$Cluster)

preservation_melanoma$Size <- cluster_sizes[match(preservation_melanoma$Label, cluster_sizes$Label),"Size"]

preservation_melanoma$Per_50 <- as.numeric(as.character(preservation_melanoma$Per_50))

#save(preservation_melanoma, file="preservation_melanoma.Rdata")
