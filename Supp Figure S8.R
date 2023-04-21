#Script for Figure S8.

library(reshape2)
library(ggraph)
library(igraph)
library(gridExtra)
library(ggrepel)

source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")


##Gene ages
genes_phy <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_phy_categorical <- genes_phy
genes_phy_categorical$Phylostrata <- ifelse(genes_phy_categorical$Phylostrata %in% 1:3, "UC",
                                            ifelse(genes_phy_categorical$Phylostrata %in% 4:9, "EM",
                                                   ifelse(genes_phy_categorical$Phylostrata %in% 10:16, "MM", NA)))


###See differences in strength of correlation of UC, MC and mixed modules, and by preservation
UC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "UC", "GeneID"])
MC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata != "UC", "GeneID"])


load("age_enrichment.Rdata")

#Subnetworks here correspond to the networks of each of the modules
number_subnetworks <- vector()
for(tumour in tumours){
  for(tissue_type in c("normal", "tumour")){
    if((tissue_type == "tumour") | (tissue_type == "normal" & tumour %in% tumours_normal)){
      load(paste("Subnetworks_", tumour, "_", tissue_type,
                 ".Rdata", sep=""))
      local_subnet <- sub_networks[[tissue_type]][[tumour]]
      
      local_ages <- age_enrichment[[tumour]]
      local_ages <- local_ages[local_ages$tissue_type == tissue_type,]
      
      number <- calculate_number_of_connections(local_subnet, tumour, tissue_type)
      
      number$Module_age <- local_ages[match(number$Module, local_ages$cluster), "Module_age"]
      number_subnetworks <- rbind(number_subnetworks, number)
    }
  }
  print(tumour)
}

#save(number_subnetworks, file="number_subnetworks.Rdata")

load("number_subnetworks.Rdata")
number_subnetworks_melt <- melt(number_subnetworks)

number_subnetworks_melt$Module_age <- factor(number_subnetworks_melt$Module_age,
                                             levels=c("UC", "Mixed", "MC"))
colnames(number_subnetworks_melt)[5] <- "Connection_type"
colnames(number_subnetworks_melt)[6] <- "Percentage"

number_subnetworks_melt$Connection_type <- gsub("_connections", "", number_subnetworks_melt$Connection_type)
number_subnetworks_melt$Connection_type <- factor(number_subnetworks_melt$Connection_type,
                                                   levels=c("UC_UC", "UC_MC", "MC_MC"))

pdf("Figure_S8.pdf", height=3, width=7)
g <- ggplot(number_subnetworks_melt, aes(x=Connection_type, y=Percentage))+
  geom_boxplot(aes(fill=Connection_type))+
  facet_grid(.~Module_age)+
  ylab("Percentage of connections")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()
