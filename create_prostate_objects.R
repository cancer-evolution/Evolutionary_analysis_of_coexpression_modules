library(WGCNA)
library(ggplot2)

path_to_clinical <- "Clinical_information_2019/"
TCGA_clinical <- read.delim(paste(path_to_clinical, "PRAD.clin.merged.txt", sep=""))

TCGA_clinical <- TCGA_clinical[TCGA_clinical[,1] %in% c("patient.bcr_patient_barcode", 
                                                        "patient.stage_event.gleason_grading.gleason_score",
                                                        "patient.stage_event.gleason_grading.primary_pattern",
                                                        "patient.stage_event.gleason_grading.secondary_pattern"),]
TCGA_clinical <- t(TCGA_clinical)
colnames(TCGA_clinical) <- c("Patient_barcode", "Gleason", "Pattern_1", "Pattern_2")
TCGA_clinical <- TCGA_clinical[-1,]
TCGA_clinical <- as.data.frame(TCGA_clinical)
TCGA_clinical$Patient <- substr(TCGA_clinical$Patient_barcode, 9, 12)


TCGA_clinical$Grade_group <- ifelse(TCGA_clinical$Pattern_1 == 3 & TCGA_clinical$Pattern_2 == 3, "GG1",
                                    ifelse(TCGA_clinical$Pattern_1 == 2 & TCGA_clinical$Pattern_2 == 4, "GG1",
                                      ifelse(TCGA_clinical$Pattern_1 == 3 & TCGA_clinical$Pattern_2 == 4, "GG2",
                                           ifelse(TCGA_clinical$Pattern_1 == 4 & TCGA_clinical$Pattern_2 == 3, "GG3",
                                                  ifelse(TCGA_clinical$Pattern_1 == 4 & TCGA_clinical$Pattern_2 == 4, "GG4",
                                                         ifelse(TCGA_clinical$Pattern_1 == 3 & TCGA_clinical$Pattern_2 == 5, "GG4",
                                                                ifelse(TCGA_clinical$Pattern_1 == 5 & TCGA_clinical$Pattern_2 == 3, "GG4", "GG5")))))))



TCGA_clinical$Grade_class <- ifelse(TCGA_clinical$Grade_group %in% c("GG4", "GG5"), "High_grade", "Low_grade")

high_grade_patients <- toupper(TCGA_clinical$Patient[TCGA_clinical$Grade_class == "High_grade"])
low_grade_patients <- toupper(TCGA_clinical$Patient[TCGA_clinical$Grade_class == "Low_grade"])

source("functions.R")

tumour <- "PRAD"
library(limma)
library(edgeR)
expression_all <- read_expression_paired("PRAD")

expression_tumour <- expression_all$tumour

expression_high <- expression_tumour[,toupper(colnames(expression_tumour)) %in% high_grade_patients]
expression_low <- expression_tumour[,toupper(colnames(expression_tumour)) %in% low_grade_patients]

expression <- list()
expression[["normal"]] <- expression_all$normal
expression[["low_grade"]] <- expression_low
expression[["high_grade"]] <- expression_high


expression[["normal"]] <- t(expression[["normal"]])
expression[["low_grade"]] <- t(expression[["low_grade"]])
expression[["high_grade"]] <- t(expression[["high_grade"]])


for(type in c("normal", "low_grade", "high_grade")){
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

for(type in c("normal", "low_grade", "high_grade")){
  sampleTree = hclust(dist(expression[[type]]), method = "average");
  plot(sampleTree, main = paste("Sample clustering to detect outliers", type, sep="\n"), sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)  
}


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft <- list()
for(type in c("normal", "low_grade", "high_grade")){
  sft[[type]] = pickSoftThreshold(expression[[type]], powerVector = powers, verbose = 5,
                                  dataIsExpr = T)
}

save(sft, file="sft_paired_PRAD_2.Rdata")
load("sft_paired_PRAD_2.Rdata")


par(mfrow=c(2,2))
for(type in c("normal", "low_grade", "high_grade")){
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
soft_threshold[["normal"]] <- 12
soft_threshold[["low_grade"]] <- 7
soft_threshold[["high_grade"]] <- 16


#Adjacency matrix and TOM
minModuleSize = 30;

modules_PRAD <- list()

for(type in c("normal", "low_grade", "high_grade")){
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
    modules_PRAD[[type]][[cluster]] <- colnames(expression[[type]])[indices]
  }
  print(type)
}

save(modules_PRAD, file="modules_PRAD_2.Rdata")


##Calculate age and novelty, and centrality of genes
library(ggplot2)
library(gridExtra)
library(Biobase)
library(GO.db)
library(igraph)
library(reshape2)

source("functions.R")

load("modules_PRAD_2.Rdata")

##Number of genes
mean(sapply(modules_PRAD$normal, length))
mean(sapply(modules_PRAD$low_grade, length))
mean(sapply(modules_PRAD$high_grade, length))


gene_ages <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_in_clusters <- unique(unname(unlist(modules_PRAD)))

gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 1:3),1]
MC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 4:16),1]


percentage_UC_cluster <- vector()
for(type in c("normal", "low_grade", "high_grade")){
  all_genes <- unique(unname(unlist(modules_PRAD[[type]])))
  local_exp_per_UC <- sum(all_genes %in% UC_genes)/c(sum(all_genes %in% UC_genes) + sum(all_genes %in% MC_genes))*100
  cluster_names <- names(modules_PRAD[[type]])
  for(cluster_name in cluster_names){
    genes_in_cluster <- modules_PRAD[[type]][[cluster_name]]
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

for(type in c("normal", "low_grade", "high_grade")){
  cluster_names <- names(modules_PRAD[[type]])
  all_genes <- unique(unname(unlist(modules_PRAD[[type]])))
  UC_genes_total <- all_genes[all_genes %in% UC_genes]
  MC_genes_total <- all_genes[all_genes %in% MC_genes]
  
  for(cluster_name in cluster_names){
    genes_in_cluster <- modules_PRAD[[type]][[cluster_name]]
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

save(age_enrichment, file="age_enrichment_PRAD_2.Rdata")
#load("age_enrichment_PRAD_2.Rdata")

table(age_enrichment$normal$Module_age)

table(age_enrichment$low_grade$Module_age)

table(age_enrichment$high_grade$Module_age)

par(mfrow=c(2,2))
hist(sapply(modules_PRAD$normal, length), breaks=100, main="Normal")
hist(sapply(modules_PRAD$low_grade, length), breaks=100, main="low_grade")
hist(sapply(modules_PRAD$high_grade, length), breaks=100, main="high_grade")


##Number of genes
module_size_ages <- vector()
for(tissue in c("normal", "low_grade", "high_grade")){
  temp <- modules_PRAD[[tissue]]
  age_temp <- age_enrichment[[tissue]]
  age_temp_UC <- age_temp[age_temp$Module_age == "UC","cluster"]
  age_temp_MC <- age_temp[age_temp$Module_age == "MC","cluster"]
  age_temp_Mixed <- age_temp[age_temp$Module_age == "Mixed","cluster"]
  
  temp_UC <- temp[names(temp) %in% age_temp_UC]
  temp_MC <- temp[names(temp) %in% age_temp_MC]
  temp_Mixed <- temp[names(temp) %in% age_temp_Mixed]

  print(tissue)
  print(paste("UC:", mean(sapply(temp_UC, length))))
  
  print(paste("Mixed:", mean(sapply(temp_Mixed, length))))
  print(paste("MC:", mean(sapply(temp_MC, length))))

  module_size_ages <- rbind(module_size_ages,
                            c(tissue, mean(sapply(temp_UC, length)),
                              mean(sapply(temp_Mixed, length)),
                              mean(sapply(temp_MC, length))))
    
}
colnames(module_size_ages) <- c("Tissue", "UC_modules", "Mixed_modules", "MC_modules")
module_size_ages <- as.data.frame(module_size_ages)
module_size_ages$UC_modules <- as.numeric(as.character(module_size_ages$UC_modules))
module_size_ages$Mixed_modules <- as.numeric(as.character(module_size_ages$Mixed_modules))
module_size_ages$MC_modules <- as.numeric(as.character(module_size_ages$MC_modules))
module_size_ages_melt <- melt(module_size_ages)
module_size_ages_melt$variable <- gsub("_modules", "", module_size_ages_melt$variable)

module_size_ages_melt$Tissue <- factor(module_size_ages_melt$Tissue,
                                         levels=c("normal", "low_grade", "high_grade"))

module_size_ages_melt$variable <- factor(module_size_ages_melt$variable,
                                         levels=c("UC", "Mixed", "MC"))

ggplot(module_size_ages_melt, aes(x=variable, y=value))+
  geom_bar(stat='identity', aes(fill=variable))+
  xlab("Module age")+
  ylab("Mean module size")+
  facet_grid(.~Tissue)


mean(sapply(modules_PRAD$normal, length))
mean(sapply(modules_PRAD$low_grade, length))
mean(sapply(modules_PRAD$high_grade, length))



number_of_shared_genes <- vector()
for(type1 in c("normal", "low_grade", "high_grade")){
  for(type2 in c("normal", "low_grade", "high_grade")){
    if(type1 != type2){
      clusters_of_local_sample_type <- modules_PRAD[[type1]]
      clusters_of_other_sample_type <- modules_PRAD[[type2]]
      
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

save(number_of_shared_genes, file="number_of_shared_genes_PRAD_2.Rdata")

##Novelty
source("functions.R")
load("number_of_shared_genes_PRAD_2.Rdata")
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module1 != "grey",]
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Module2 != "grey",]

preservation <- calculate_preservation_PRAD(number_of_shared_genes)

save(preservation, file="preservation_PRAD_2.Rdata")
load("preservation_PRAD_2.Rdata")

cluster_sizes <- vector()
for(type in c("normal", "low_grade", "high_grade")){
  temp <- modules_PRAD[[type]]
  mod <- sapply(temp, length)
  mod <- cbind(Module=names(mod), Size=as.vector(mod), Type=type)
  cluster_sizes <- rbind(cluster_sizes, mod)
}

cluster_sizes <- as.data.frame(cluster_sizes)
cluster_sizes$Size <- as.numeric(as.character(cluster_sizes$Size))
cluster_sizes$Label <- paste(cluster_sizes$Type, cluster_sizes$Module)

preservation$Label <- paste(preservation$Type, preservation$Cluster)

preservation$Size <- cluster_sizes[match(preservation$Label, cluster_sizes$Label),"Size"]

preservation$Per_50 <- as.numeric(as.character(preservation$Per_50))

save(preservation, file="preservation_PRAD_2.Rdata")



load("age_enrichment_PRAD_2.Rdata")



all_preservation <- preservation[,c("Type", "Cluster", "Per_50", "Label", "Size")]
all_preservation$Preservation_ratio <- all_preservation$Per_50/all_preservation$Size


##Age
age_enrichment_df <- vector()
for(type in c("normal", "low_grade", "high_grade")){
  local_age_enrichment <- age_enrichment[[type]]
  age_enrichment_df <- rbind(age_enrichment_df, local_age_enrichment)
}
age_enrichment_df$Label <- paste(age_enrichment_df$type, age_enrichment_df$cluster)

all_preservation$Age <- age_enrichment_df$Module_age[match(all_preservation$Label, 
                                                                age_enrichment_df$Label)]
all_preservation$Age <- factor(all_preservation$Age, levels=c("UC", "Mixed", "MC"))
all_preservation$Type <- as.character(all_preservation$Type)

all_preservation$Type[all_preservation$Type == "high_grade"] <- "low_to_high_grade"
all_preservation$Type[all_preservation$Type == "low_grade"] <- "normal_to_low_grade"
all_preservation$Type <- factor(all_preservation$Type, 
                                     levels=c("normal_to_low_grade", "low_to_high_grade"))

normal_to_low <- all_preservation[all_preservation$Type == "normal_to_low_grade",]
low_to_high <- all_preservation[all_preservation$Type == "low_to_high_grade",]
wilcox.test(normal_to_low$Preservation_ratio, low_to_high$Preservation_ratio)

ggplot(all_preservation, aes(x=Type, y=Preservation_ratio))+
  geom_boxplot()+
  xlab("Transition")+ylab("Novelty")


#Greater ratio means less preservation
library(ggplot2)
ggplot(all_preservation, aes(x=Age, y=Preservation_ratio))+
  geom_boxplot(aes(fill=Age))+
  facet_grid(.~Type)+
  ylab("Novelty")+
  theme_bw()

wilcox.test(normal_to_low$Preservation_ratio[normal_to_low$Age == "UC"],
            normal_to_low$Preservation_ratio[normal_to_low$Age != "UC"])
wilcox.test(normal_to_low$Preservation_ratio[normal_to_low$Age == "MC"],
            normal_to_low$Preservation_ratio[normal_to_low$Age != "MC"])
wilcox.test(normal_to_low$Preservation_ratio[normal_to_low$Age == "Mixed"],
            normal_to_low$Preservation_ratio[normal_to_low$Age != "Mixed"])

wilcox.test(low_to_high$Preservation_ratio[low_to_high$Age == "UC"],
            low_to_high$Preservation_ratio[low_to_high$Age != "UC"], alternative="greater")
wilcox.test(low_to_high$Preservation_ratio[low_to_high$Age == "MC"],
            low_to_high$Preservation_ratio[low_to_high$Age != "MC"], alternative="greater")
wilcox.test(low_to_high$Preservation_ratio[low_to_high$Age == "Mixed"],
            low_to_high$Preservation_ratio[low_to_high$Age == "UC"], alternative="greater")

wilcox.test(low_to_high$Preservation_ratio[low_to_high$Age == "Mixed"],
            low_to_high$Preservation_ratio[low_to_high$Age == "MC"], alternative="greater")