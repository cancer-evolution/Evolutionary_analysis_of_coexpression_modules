##Script for Figures S13 and S14

library(reshape2)
library(ggraph)
library(igraph)
library(gridExtra)
library(ggrepel)
library(readr)

source("functions.R")

load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")

load("cluster_assignments.Rdata")
load("all_preservation_t_to_n2.Rdata")
load("all_preservation_n_to_t2.Rdata")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]
tumours <- tumours[tumours != "MESO"]

genes_phy <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_phy_categorical <- genes_phy
genes_phy_categorical$Phylostrata <- ifelse(genes_phy_categorical$Phylostrata %in% 1:3, "UC",
                                            ifelse(genes_phy_categorical$Phylostrata %in% 4:9, "EM",
                                                   ifelse(genes_phy_categorical$Phylostrata %in% 10:16, "MM", NA)))

UC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "UC", "GeneID"])
MC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata != "UC", "GeneID"])

n_genes_phy <- table(genes_phy$Phylostrata)
n_genes_phy <- data.frame(Phy = names(n_genes_phy), Number = as.vector(n_genes_phy))

n_genes_phy$Phy <- factor(n_genes_phy$Phy, levels=1:16)

n_genes_phy$Age <- ifelse(n_genes_phy$Phy %in% 1:3, "UC",
                          ifelse(n_genes_phy$Phy %in% 4:9, "EM",
                                 ifelse(n_genes_phy$Phy %in% 10:16, "MM", NA)))

load("CNVs_above_fraction_0.25.Rdata")

CNVs_df <- load_CNVs(only_focal="FOCAL")

mutations_df <- load_mutations()
mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
mutations_df <- mutations_df[mutations_df$Syn_ratio >1, ]

subnet_tumour <- list()
for(tumour in tumours){
  load(paste("Subnetworks_", tumour, "_tumour.Rdata", sep=""))
  subnet_tumour[[tumour]] <- sub_networks$tumour[[tumour]]
}

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")
subnet_normal <- list()
for(tumour in tumours_normal){
  load(paste("Subnetworks_", tumour, "_normal.Rdata", sep=""))
  subnet_normal[[tumour]] <- sub_networks$normal[[tumour]]
}


load("degree_modules.Rdata")
load("age_enrichment.Rdata")

degree_nodes_all2 <- vector()
attributes_age_all2 <- vector()
attributes_CNV_all2 <- vector()
attributes_mut_all2 <- vector()

age_degree_p_all2 <- vector()
CNV_degree_p_all2 <- vector()
mut_degree_p_all2 <- vector()

for(tumour in tumours){
  load(paste("degree_nodes_all_",
             tumour, ".Rdata", sep=""))
  if(!(tumour %in% tumours_normal)){
    degree_nodes_all$Module_normal <- NA
    degree_nodes_all$Degree_normal <- NA
    degree_nodes_all$Degree_norm_normal <- NA
    degree_nodes_all$Degree_rank_normal <- NA
    degree_nodes_all$Degree_rank_norm_normal <- NA
    
  }
  
  degree_nodes_all2 <- rbind(degree_nodes_all2, degree_nodes_all)
}

degree_nodes_all <- degree_nodes_all2

degree_nodes_all$Degree_rank_norm_normal <- 1-degree_nodes_all$Degree_rank_norm_normal
degree_nodes_all$Degree_rank_norm_tumour <- 1-degree_nodes_all$Degree_rank_norm_tumour

degree_nodes_all$Tumour <- factor(degree_nodes_all$Tumour, levels=tumours)

degree_nodes_all$Alt <- factor(degree_nodes_all$Alt,
                               levels=c("None", "Amp", "Del", "Miss", "LoF"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

degree_nodes_all2 <- degree_nodes_all[,c("Tumour", "Alt", "Degree_rank_norm_normal", "Degree_rank_norm_tumour")]

degree_nodes_all2 <- melt(degree_nodes_all2)
colnames(degree_nodes_all2)[3:4] <- c("Tissue", "Centrality")
degree_nodes_all2$Tissue <- as.character(degree_nodes_all2$Tissue)
degree_nodes_all2$Tissue[degree_nodes_all2$Tissue == "Degree_rank_norm_normal"] <- "Normal"
degree_nodes_all2$Tissue[degree_nodes_all2$Tissue == "Degree_rank_norm_tumour"] <- "Tumour"

degree_nodes_all2 <- subset(degree_nodes_all2, Alt %in% c("Amp", "Del", "Miss", "LoF", "None"))


degree_nodes_all$Gene_age <- genes_phy$Phylostrata[match(degree_nodes_all$Gene, genes_phy$GeneID)]
degree_nodes_all$Gene_age[degree_nodes_all$Gene_age %in% c(1,2,3)] <- "UC"
degree_nodes_all$Gene_age[degree_nodes_all$Gene_age %in% c(4,5,6,7,8,9)] <- "EM"
degree_nodes_all$Gene_age[degree_nodes_all$Gene_age %in% c(10:16)] <- "MM"

temp <- subset(degree_nodes_all, Alt %in% c("Amp", "None")
               & !(Tumour %in% c("ACC", "COAD", "KIRP")))##No amplifications


##S13
pdf("Figure_S13.pdf", height=2, width=17)
g <- ggplot(temp, aes(x=Degree_rank_norm_tumour))+
  geom_density(aes(fill=Alt))+
  facet_grid(.~Tumour)+
  xlab("Normalized degree")+
  ggtitle("Tumour")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()

temp <- subset(degree_nodes_all, Alt %in% c("Del", "None")
               & !(Tumour %in% c("UVM")))

##S14
pdf("Figure_S14.pdf", height=2, width=17)
g <- ggplot(temp, aes(x=Degree_rank_norm_tumour))+
  geom_density(aes(fill=Alt))+
  facet_grid(.~Tumour)+
  xlab("Normalized degree")+
  ggtitle("Tumour")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()
