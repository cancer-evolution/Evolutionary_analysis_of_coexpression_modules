library(reshape2)
library(readr)

source("functions.R")

load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")
load("all_preservation_t_to_n2.Rdata")
load("cluster_assignments.Rdata")
load("CNVs_above_fraction_0.25.Rdata")
load("age_enrichment.Rdata")

tumours <- list.files(path = "/home/atrigos/TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "MESO"]  #Excluded mesothelioma because no CNVs
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]


genes_phy <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_phy_categorical <- genes_phy
genes_phy_categorical$Phylostrata <- ifelse(genes_phy_categorical$Phylostrata %in% 1:3, "UC",
                                            ifelse(genes_phy_categorical$Phylostrata %in% 4:9, "EM",
                                                   ifelse(genes_phy_categorical$Phylostrata %in% 10:16, "MM", NA)))
n_genes_phy <- table(genes_phy$Phylostrata)
n_genes_phy <- data.frame(Phy = names(n_genes_phy), Number = as.vector(n_genes_phy))

n_genes_phy$Phy <- factor(n_genes_phy$Phy, levels=1:16)

n_genes_phy$Age <- ifelse(n_genes_phy$Phy %in% 1:3, "UC",
                          ifelse(n_genes_phy$Phy %in% 4:9, "EM",
                                 ifelse(n_genes_phy$Phy %in% 10:16, "MM", NA)))

all_preservation_t_to_n <- all_preservation_t_to_n2

CNVs_df <- load_CNVs(only_focal="FOCAL")
amp_genes <- CNVs_df[CNVs_df$Patients_amp >= 0.1,c("Genes", "Tumour")]
del_genes <- CNVs_df[CNVs_df$Patients_del >= 0.1,c("Genes", "Tumour")]

tumours <- c(tumours, "MESO")
mutations_df <- load_mutations()

mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
mutations_df <- mutations_df[mutations_df$Syn_ratio >1, ]

per_somatic_modules <- vector()
for(tumour in tumours){
  local_tumour_modules <- cluster_assignments$tumour[[tumour]]
  local_normal_modules <- cluster_assignments$normal[[tumour]]
  
  local_point <- mutations_df[mutations_df$Tumour == tumour,]
  
  local_amp <- amp_genes[amp_genes$Tumour == tumour,]
  local_del <- del_genes[del_genes$Tumour == tumour,]
  
  for(mod in names(local_tumour_modules)){
    local_genes <- local_tumour_modules[[mod]]
    per_miss <- (sum(local_point[local_point$Variant_type == "Missense", "Hugo_Symbol"] %in% local_genes)/length(local_genes))*100
    per_lof <- (sum(local_point[local_point$Variant_type == "LoF", "Hugo_Symbol"] %in% local_genes)/length(local_genes))*100
    per_amp <- (sum(local_amp$Genes %in% local_genes)/length(local_genes))*100
    per_del <- (sum(local_del$Genes %in% local_genes)/length(local_genes))*100
    per_somatic_modules <- rbind(per_somatic_modules,
                                 c(tumour, mod, per_miss, per_lof, per_amp, per_del, paste(tumour, mod)))
  }
  
}
per_somatic_modules <- as.data.frame(per_somatic_modules)
colnames(per_somatic_modules) <- c("Tumour", "Module", "Per_miss", "Per_LoF", "Per_amp", "Per_del", "Label")

module_information <- per_somatic_modules
all_preservation_t_to_n$Label <- paste(all_preservation_t_to_n$Tumour, all_preservation_t_to_n$Cluster)
module_information$Novelty <- all_preservation_t_to_n$Preservation_ratio[match(module_information$Label, all_preservation_t_to_n$Label)]
module_information$Novelty_category <- all_preservation_t_to_n$Category[match(module_information$Label, all_preservation_t_to_n$Label)]

age_enrichment_df <- vector()
for(tumour in names(age_enrichment)){
  age_enrichment_df <- rbind(age_enrichment_df, age_enrichment[[tumour]])
}
age_enrichment_df$Label <- paste(age_enrichment_df$tumour, age_enrichment_df$cluster)

age_enrichment_df_tumour <- age_enrichment_df[age_enrichment_df$tissue_type == "tumour",]

module_information$Age <- age_enrichment_df_tumour$Module_age[match(module_information$Label,
                                                                    age_enrichment_df_tumour$Label)]
module_information$Age <- factor(module_information$Age, levels=c("UC", "Mixed", "MC"))

module_information$Per_miss <- per_somatic_modules[match(module_information$Label, per_somatic_modules$Label), "Per_miss"]
module_information$Per_lof <- per_somatic_modules[match(module_information$Label, per_somatic_modules$Label), "Per_LoF"]
module_information$Per_amp <- per_somatic_modules[match(module_information$Label, per_somatic_modules$Label), "Per_amp"]
module_information$Per_del <- per_somatic_modules[match(module_information$Label, per_somatic_modules$Label), "Per_del"]

module_information$Per_miss <- as.numeric(as.character(module_information$Per_miss))
module_information$Per_lof <- as.numeric(as.character(module_information$Per_lof))
module_information$Per_amp <- as.numeric(as.character(module_information$Per_amp))
module_information$Per_del <- as.numeric(as.character(module_information$Per_del))

module_information$Novelty_category <- as.character(module_information$Novelty_category)
module_information$Novelty_category[module_information$Novelty_category == "Low_score"] <- "Low"
module_information$Novelty_category[module_information$Novelty_category == "Median_score"] <- "Moderate"
module_information$Novelty_category[module_information$Novelty_category == "High_score"] <- "High"
module_information$Novelty_category[module_information$Novelty_category == "Inf_score"] <- "High"
module_information$Novelty_category <- factor(module_information$Novelty_category,
                                              levels=c("Low", "Moderate", "High"))

save(module_information, file="somatic_muta_in_modules.Rdata")
