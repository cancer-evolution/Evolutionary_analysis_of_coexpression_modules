##Script for Figure 3, panels A-F, and Figures S10-S12.

library(ggplot2)
library(ggsci)
library(clinfun)

##Panel A
source("functions.R")

load("cluster_assignments.Rdata")
load("all_preservation_t_to_n2.Rdata")
all_preservation_t_to_n <- all_preservation_t_to_n2

##2022 version (downloaded Dec 22, 2022)
cancer_genes <- read.csv("Census_allWed Dec 21 23_46_30 2022.csv")
cancer_genes$Gene.Symbol <- toupper(cancer_genes$Gene.Symbol)

#Using the cancer genes defined for each tumour types

#Checking which regular expressions to use
all_tumours <- unique(trim(unlist(strsplit(as.character(cancer_genes$Tumour.Types.Somatic.), ","))))
all_tumours[grep("adrenocortical", all_tumours)]
all_tumours[grep("^bladder", all_tumours)]
all_tumours[grep("breast", all_tumours)]
all_tumours[grep("cervical", all_tumours)]
all_tumours[grep("cholangio", all_tumours)]
all_tumours[grep("colon", all_tumours)]
all_tumours[grep("colorectal", all_tumours)]
all_tumours[grep("esopha", all_tumours)]
all_tumours[grep("^glioblastoma", all_tumours)]
all_tumours[grep("head", all_tumours)]
all_tumours[grep("kidney", all_tumours)]
all_tumours[grep("hepatocellular", all_tumours)]
all_tumours[grep("lung", all_tumours)]
all_tumours[grep("meso", all_tumours)]
all_tumours[grep("ovarian", all_tumours)]
all_tumours[grep("pancrea", all_tumours)]
all_tumours[grep("pheo", all_tumours)]
all_tumours[grep("paragan", all_tumours)]
all_tumours[grep("prostate", all_tumours)]
all_tumours[grep("rectal", all_tumours)]
all_tumours[grep("sarcoma", all_tumours)]
all_tumours[grep("melano", all_tumours)] ##minus "uveal melanoma"
all_tumours[c(grep("stomach", all_tumours),
              grep("gastric", all_tumours))]
all_tumours[grep("testicular", all_tumours)]
all_tumours[grep("thyroid", all_tumours)]
all_tumours[c(grep("^endometrial$", all_tumours),
              grep("endometrial stromal tumour", all_tumours),
              grep("endometrial carcinoma", all_tumours))]
all_tumours[grep("uterine serous carcinoma", all_tumours)]
all_tumours[grep("uveal melano", all_tumours)]


cancer_genes$Tumour.Types.Somatic. <- as.character(cancer_genes$Tumour.Types.Somatic.)

cancer_genes_per_tumour <- list()
cancer_genes_per_tumour[["ACC"]] <- cancer_genes[grep("adrenocortical", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["BLCA"]] <- cancer_genes[grep("^bladder", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["BRCA"]] <- cancer_genes[grep("breast", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["CESC"]] <- cancer_genes[grep("cervical", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["CHOL"]] <- cancer_genes[grep("cholangio", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["COAD"]] <- cancer_genes[c(grep("colon", cancer_genes$Tumour.Types.Somatic.),
                                                    grep("colorectal", cancer_genes$Tumour.Types.Somatic.)), "Gene.Symbol"]
cancer_genes_per_tumour[["ESCA"]] <- cancer_genes[grep("esopha", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["GBM"]] <- cancer_genes[grep("^glioblastoma", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["HNSC"]] <- cancer_genes[grep("head", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["KICH"]] <- cancer_genes[grep("kidney", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["KIRC"]] <- cancer_genes[grep("kidney", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["KIRP"]] <- cancer_genes[grep("kidney", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]

cancer_genes_per_tumour[["LIHC"]] <- cancer_genes[grep("hepatocellular", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["LUAD"]] <- cancer_genes[grep("lung", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["LUSC"]] <- cancer_genes[grep("lung", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]

cancer_genes_per_tumour[["MESO"]] <- cancer_genes[grep("meso", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["OV"]] <- cancer_genes[grep("ovarian", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["PAAD"]] <- cancer_genes[grep("pancrea", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["PCPG"]] <- cancer_genes[c(grep("pheo", cancer_genes$Tumour.Types.Somatic.),
                                                    grep("paragan", cancer_genes$Tumour.Types.Somatic.)), "Gene.Symbol"]

cancer_genes_per_tumour[["PRAD"]] <- cancer_genes[grep("prostate", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["READ"]] <- cancer_genes[grep("rectal", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["SARC"]] <- cancer_genes[grep("sarcoma", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["SKCM"]] <- cancer_genes[grep("melano", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"] ##minus "uveal melanoma"
cancer_genes_per_tumour[["STAD"]] <- cancer_genes[c(grep("stomach", cancer_genes$Tumour.Types.Somatic.),
                                                    grep("gastric", cancer_genes$Tumour.Types.Somatic.)), "Gene.Symbol"]
cancer_genes_per_tumour[["TGCT"]] <- cancer_genes[grep("testicular", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["THCA"]] <- cancer_genes[grep("thyroid", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["UCEC"]] <- cancer_genes[c(grep("^endometrial$", cancer_genes$Tumour.Types.Somatic.),
                                                    grep("endometrial stromal tumour", cancer_genes$Tumour.Types.Somatic.),
                                                    grep("endometrial carcinoma", cancer_genes$Tumour.Types.Somatic.)), "Gene.Symbol"]
cancer_genes_per_tumour[["UCS"]] <- cancer_genes[grep("uterine serous carcinoma", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_genes_per_tumour[["UVM"]] <- cancer_genes[grep("uveal melano", cancer_genes$Tumour.Types.Somatic.), "Gene.Symbol"]


tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")
all_preservation_t_to_n_cc <- vector()

for(tumour in tumours_normal){
  all_preservation_t_to_n_cc <- rbind(all_preservation_t_to_n_cc,
                                      add_number_cc_genes(all_preservation_t_to_n2, tumour, cancer_genes_per_tumour))
}
all_preservation_t_to_n_cc <- as.data.frame(all_preservation_t_to_n_cc)

all_preservation_t_to_n_cc$Tumour <- factor(all_preservation_t_to_n_cc$Tumour,
                                            levels=tumours_normal)

all_preservation_t_to_n_cc_all <- all_preservation_t_to_n_cc

all_preservation_t_to_n_cc_all <- all_preservation_t_to_n_cc_all[all_preservation_t_to_n_cc_all$N_cancer_census != 0,]

all_preservation_t_to_n_cc_all$Category <- factor(all_preservation_t_to_n_cc_all$Category,
                                                  levels=c("High_score", "Median_score", "Low_score"))

##Panel A
pdf("Figure_3A.pdf", height=3.2, width=3.5)
g <- ggplot(all_preservation_t_to_n_cc_all, aes(x=Category, y=Per_cancer_census))+
  geom_boxplot(aes(fill=Category))+
  scale_fill_jco()+
  ylab("Percentage\ncancer census genes")+
  xlab("Module novelty")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
print(g)
dev.off()

wilcox.test(all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Category == "High_score"],
all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Category == "Median_score"],
alternative="greater")
wilcox.test(all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Category == "High_score"],
            all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Category == "Low_score"],
            alternative="greater")


#Figure S10
pdf("Figure_S10.pdf", height=3.85, width=5)
g <- ggplot(all_preservation_t_to_n_cc_all, aes(y=Preservation_ratio, x=Per_cancer_census))+
  geom_point(aes(colour=Category))+
  scale_colour_jco()+
  geom_smooth(se=FALSE, colour="black", method="lm", size=0.25)+
  #facet_grid(.~Tumour)+
  xlab("Percentage of\nknown cancer drivers\nin module")+
  ylab("Novelty score of module")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

cor.test(all_preservation_t_to_n_cc_all$Per_cancer_census, all_preservation_t_to_n_cc_all$Preservation_ratio, method="sp")$p.value


##Panel B
pdf("Figure_3B.pdf", height=3.2, width=3.5)
g <- ggplot(all_preservation_t_to_n_cc_all, aes(x=Age, y=Per_cancer_census))+
  geom_boxplot(aes(fill=Age))+
  ylab("Percentage\ncancer census genes")+
  xlab("Module age")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
print(g)
dev.off()

wilcox.test(all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Age == "Mixed"],
            all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Age == "UC"], alternative="greater")

wilcox.test(all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Age == "Mixed"],
            all_preservation_t_to_n_cc_all$Per_cancer_census[all_preservation_t_to_n_cc_all$Age == "MC"], alternative="greater")


##Panel D & E
library(ggplot2)
library(reshape2)
library("ggsci")
library(readr)
library(gridExtra)

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
tumours <- tumours[tumours != "MESO"]  #Had to exclude mesothelioma because no CNVs
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")


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

age_enrichment_df <- vector()
for(tumour in names(age_enrichment)){
  age_enrichment_df <- rbind(age_enrichment_df, age_enrichment[[tumour]])
}
age_enrichment_df$Label <- paste(age_enrichment_df$tumour, age_enrichment_df$cluster)

age_enrichment_df_tumour <- age_enrichment_df[age_enrichment_df$tissue_type == "tumour",]


load("somatic_muta_in_modules.Rdata")
module_information$Module_age <- module_information$Age
module_information <- module_information[module_information$Module != "grey",]


per_somatic_modules_normal <- vector()
for(tumour in tumours){
  local_normal_modules <- cluster_assignments$normal[[tumour]]
  local_point <- mutations_df[mutations_df$Tumour == tumour,]
  local_amp <- amp_genes[amp_genes$Tumour == tumour,]
  local_del <- del_genes[del_genes$Tumour == tumour,]
  
  for(mod in names(local_normal_modules)){
    local_genes <- local_normal_modules[[mod]]
    per_miss <- (sum(local_point[local_point$Variant_type == "Missense", "Hugo_Symbol"] %in% local_genes)/length(local_genes))*100
    per_lof <- (sum(local_point[local_point$Variant_type == "LoF", "Hugo_Symbol"] %in% local_genes)/length(local_genes))*100
    per_amp <- (sum(local_amp$Genes %in% local_genes)/length(local_genes))*100
    per_del <- (sum(local_del$Genes %in% local_genes)/length(local_genes))*100
    per_somatic_modules_normal <- rbind(per_somatic_modules_normal,
                                        c(tumour, mod, per_miss, per_lof, per_amp, per_del))
  }
}
per_somatic_modules_normal <- as.data.frame(per_somatic_modules_normal)
colnames(per_somatic_modules_normal) <- c("Tumour", "Module", "Per_miss", "Per_lof", "Per_amp", "Per_del")
per_somatic_modules_normal$Per_miss <- as.numeric(as.character(per_somatic_modules_normal$Per_miss))
per_somatic_modules_normal$Per_lof <- as.numeric(as.character(per_somatic_modules_normal$Per_lof))
per_somatic_modules_normal$Per_amp <- as.numeric(as.character(per_somatic_modules_normal$Per_amp))
per_somatic_modules_normal$Per_del <- as.numeric(as.character(per_somatic_modules_normal$Per_del))

per_somatic_modules_normal$Novelty <- NA
per_somatic_modules_normal$Novelty_category <- NA

age_enrichment_df_normal <- age_enrichment_df[age_enrichment_df$tissue_type == "normal",]
per_somatic_modules_normal$Age <- age_enrichment_df_normal$Module_age[match(per_somatic_modules_normal$Module,
                                                                            age_enrichment_df_normal$cluster)]
per_somatic_modules_normal$Age <- factor(per_somatic_modules_normal$Age, levels=c("UC", "Mixed", "MC"))
per_somatic_modules_normal$Tissue <- "Normal"

per_somatic_modules_normal <- per_somatic_modules_normal[per_somatic_modules_normal$Module != "grey",]

module_information_tumor <- module_information[,c("Tumour", "Module", "Per_miss", "Per_lof", "Per_amp", "Per_del", "Novelty", "Novelty_category", "Age")]
module_information_tumor$Tissue <- "Tumour"

module_information_tumor$Novelty_category <- factor(module_information_tumor$Novelty_category,
                                                    levels=c("High", "Moderate", "Low"))
module_information_tumor$Age <- factor(module_information_tumor$Age,
                                       levels=c("UC", "Mixed", "MC"))

module_information_both <- rbind(per_somatic_modules_normal, module_information_tumor)
module_information_both$Label <- paste(module_information_both$Tumour, module_information_both$Tissue, module_information_both$Module, sep="_")

label_order_amp <- module_information_both$Label[order(module_information_both$Per_amp, decreasing=TRUE)]
label_order_del <- module_information_both$Label[order(module_information_both$Per_del, decreasing=TRUE)]
label_order_miss <- module_information_both$Label[order(module_information_both$Per_miss, decreasing=TRUE)]
label_order_lof <- module_information_both$Label[order(module_information_both$Per_lof, decreasing=TRUE)]


module_information_both$Label <- factor(module_information_both$Label,
                                        levels=label_order_amp)
module_information_both$Novelty_category <- factor(module_information_both$Novelty_category,
                                                   levels=c("High", "Moderate", "Low"))
module_information_both_amp <- subset(module_information_both, Per_amp > 0)
module_information_both_amp$Label <- factor(module_information_both_amp$Label,
                                            levels=label_order_amp)
module_information_both_amp$Amp_10 <- ifelse(module_information_both_amp$Per_amp >= 10, "More_10", "Less_10")

table(module_information_both_amp[,c("Age", "Tissue", "Amp_10")])

table(module_information_both_amp[,c("Age", "Tumour", "Amp_10")])

count_more_10_amp_age <- vector()
for(tumour in tumours_normal){
  temp <- module_information_both_amp[module_information_both_amp$Tumour == tumour,]
  if(nrow(temp) > 0){
    temp_normal <- temp[temp$Tissue == "Normal",]
    temp_tumour <- temp[temp$Tissue == "Tumour",]
    print(tumour)
    
    tumour_count <- table(temp_tumour[,c("Age", "Amp_10")])
    tumour_count <- melt(tumour_count)
    tumour_count <- tumour_count[tumour_count$Amp_10 == "More_10",]
    if(nrow(tumour_count) > 0){
      tumour_count$Percentage <- (tumour_count$value/sum(tumour_count$value))*100
      tumour_count$Tissue <- "Tumour"
      tumour_count$Tumour <- tumour
      count <- tumour_count
    }
    
    if(nrow(temp_normal) > 0){
        normal_count <- table(temp_normal[,c("Age", "Amp_10")])
        normal_count <- melt(normal_count)
        normal_count <- normal_count[normal_count$Amp_10 == "More_10",]
        if(nrow(normal_count) > 0){
          normal_count$Percentage <- (normal_count$value/sum(normal_count$value))*100
          normal_count$Tissue <- "Normal"
          normal_count$Tumour <- tumour
          count <- rbind(count, normal_count)
        }
    }
    
    if(!is.null(nrow(count))){
      count_more_10_amp_age <- rbind(count_more_10_amp_age, count)
      count <- NA
    }
  }
}

#Panel D
pdf("Figure_3D_amp.pdf", height=2, width=6)
g <- ggplot(count_more_10_amp_age, aes(x=Tumour, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Age))+
  ylab("Percentage\nof modules")+
  ggtitle("Amp")+
  facet_grid(Tissue~., scales="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

count_more_10_amp_novelty <- vector()
for(tumour in tumours_normal){
  temp <- module_information_both_amp[module_information_both_amp$Tumour == tumour,]
  if(nrow(temp) > 0){
    temp_tumour <- temp[temp$Tissue == "Tumour",]
    
    tumour_count <- table(temp_tumour[,c("Novelty_category", "Amp_10")])
    tumour_count <- melt(tumour_count)
    tumour_count <- tumour_count[tumour_count$Amp_10 == "More_10",]
    if(nrow(tumour_count) > 0){
      tumour_count$Percentage <- (tumour_count$value/sum(tumour_count$value))*100
      tumour_count$Tissue <- "Tumour"
      tumour_count$Tumour <- tumour
      count <- tumour_count
      if(nrow(count) > 0){
        count_more_10_amp_novelty <- rbind(count_more_10_amp_novelty, count)
      }
    }
  }
}

#Panel E
pdf("Figure_3E_amp.pdf", height=2, width=6)
g <- ggplot(count_more_10_amp_novelty, aes(x=Tumour, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Novelty_category))+
  ylab("Percentage of modules")+
  scale_fill_jco()+
  ggtitle("Amp")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

module_information_both_del <- subset(module_information_both, Per_del > 0)
module_information_both_del$Label <- factor(module_information_both_del$Label,
                                            levels=label_order_del)
module_information_both_del$del_10 <- ifelse(module_information_both_del$Per_del >= 10, "More_10", "Less_10")

count_more_10_del_age <- vector()
for(tumour in tumours_normal){
  temp <- module_information_both_del[module_information_both_del$Tumour == tumour,]
  if(nrow(temp) > 0){
    temp_normal <- temp[temp$Tissue == "Normal",]
    temp_tumour <- temp[temp$Tissue == "Tumour",]
    print(tumour)
    
    tumour_count <- table(temp_tumour[,c("Age", "del_10")])
    tumour_count <- melt(tumour_count)
    tumour_count <- tumour_count[tumour_count$del_10 == "More_10",]
    if(nrow(tumour_count) > 0){
      tumour_count$Percentage <- (tumour_count$value/sum(tumour_count$value))*100
      tumour_count$Tissue <- "Tumour"
      tumour_count$Tumour <- tumour
      count <- tumour_count
    }
    
    if(nrow(temp_normal) > 0){
      normal_count <- table(temp_normal[,c("Age", "del_10")])
      normal_count <- melt(normal_count)
      normal_count <- normal_count[normal_count$del_10 == "More_10",]
      if(nrow(normal_count) > 0){
        normal_count$Percentage <- (normal_count$value/sum(normal_count$value))*100
        normal_count$Tissue <- "Normal"
        normal_count$Tumour <- tumour
        count <- rbind(count, normal_count)
      }
    }
    
    if(!is.null(nrow(count))){
      count_more_10_del_age <- rbind(count_more_10_del_age, count)
      count <- NA
    }
  }
}


#Panel D
pdf("Figure_3D_del.pdf", height=2, width=6)
g <- ggplot(count_more_10_del_age, aes(x=Tumour, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Age))+
  ylab("Percentage\nof modules")+
  ggtitle("Del")+
  facet_grid(Tissue~., scales="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

count_more_10_del_novelty <- vector()
for(tumour in tumours_normal){
  temp <- module_information_both_del[module_information_both_del$Tumour == tumour,]
  if(nrow(temp) > 0){
    temp_tumour <- temp[temp$Tissue == "Tumour",]
    
    tumour_count <- table(temp_tumour[,c("Novelty_category", "del_10")])
    tumour_count <- melt(tumour_count)
    tumour_count <- tumour_count[tumour_count$del_10 == "More_10",]
    if(nrow(tumour_count) > 0){
      tumour_count$Percentage <- (tumour_count$value/sum(tumour_count$value))*100
      tumour_count$Tissue <- "Tumour"
      tumour_count$Tumour <- tumour
      count <- tumour_count
      if(nrow(count) > 0){
        count_more_10_del_novelty <- rbind(count_more_10_del_novelty, count)
      }
    }
  }
}

count_more_10_del_novelty$Novelty_category <- factor(count_more_10_del_novelty$Novelty_category,
                                                     levels=c("High","Moderate","Low"))



#Panel E
pdf("Figure_3E_del.pdf", height=2, width=6)
g <- ggplot(count_more_10_del_novelty, aes(x=Tumour, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Novelty_category))+
  ylab("Percentage of modules")+
  scale_fill_jco()+
  ggtitle("Del")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

table(module_information_both_del[,c("Novelty_category", "Tissue", "del_10")])
table(module_information_both_del[,c("Age", "Tissue", "del_10")])


#Panel C
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

tumours <- list.files(path = "/home/atrigos/TCGA_all_data")
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
tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

# subnet_tumour <- list()
# for(tumour in tumours){
#   load(paste("Subnetworks_", tumour, "_tumour.Rdata", sep=""))
#   subnet_tumour[[tumour]] <- sub_networks$tumour[[tumour]]
# }
# 
# subnet_normal <- list()
# for(tumour in tumours_normal){
#   load(paste("Subnetworks_", tumour, "_normal.Rdata", sep=""))
#   subnet_normal[[tumour]] <- sub_networks$normal[[tumour]]
# }
# 
# load("degree_modules.Rdata")
# load("age_enrichment.Rdata")
# 
# for(tumour in tumours){
#   degree_nodes_all <- vector()
#   attributes_age_all <- vector()
#   attributes_CNV_all <- vector()
#   attributes_mut_all <- vector()
#   
#   age_degree_p_all <- vector()
#   CNV_degree_p_all <- vector()
#   mut_degree_p_all <- vector()
#   load(paste("Subnetworks_",
#              tumour, "_tumour.Rdata", sep=""))
#   subnet_tumour <- sub_networks$tumour[[tumour]]
#   
#   
#   if(tumour %in% tumours_normal){
#     load(paste("Subnetworks_",
#                tumour, "_normal.Rdata", sep=""))
#     subnet_normal <- sub_networks$normal[[tumour]]
#   }
#   
#   #Somatic mutations
#   local_mut <- mutations_df[mutations_df$Tumour == tumour,]
#   genes_lof <- as.character(local_mut[local_mut$Variant_type == "LoF","Hugo_Symbol"])
#   genes_miss <- as.character(local_mut[local_mut$Variant_type == "Missense","Hugo_Symbol"])
#   
#   #CNVs
#   local_CNVs <- CNVs_df[CNVs_df$Tumour == tumour,]
#   genes_amp <- as.character(local_CNVs[local_CNVs$Patients_amp > 0.1, "Genes"])
#   genes_del <- as.character(local_CNVs[local_CNVs$Patients_del > 0.1, "Genes"])
#   
#   strength_tumour <- calculate_strength_of_subnetworks(subnet_tumour, tumour, "tumour")
#   
#   if(tumour %in% tumours_normal){
#     tumour_preservation <- all_preservation_t_to_n2[all_preservation_t_to_n2$Tumour == tumour,]
#     normal_preservation <- all_preservation_n_to_t2[all_preservation_n_to_t2$Tumour == tumour,]
#     
#     tumour_preservation$Label <- paste(tumour_preservation$Tumour, "Tumour",
#                                        tumour_preservation$Cluster, sep="_")
#     normal_preservation$Label <- paste(normal_preservation$Tumour, "Tumour",
#                                        normal_preservation$Cluster, sep="_")
#     
#     all_preservation <- rbind(tumour_preservation, normal_preservation)
#     strength_tumour$Novelty <- tumour_preservation[match(strength_tumour$Module, tumour_preservation$Cluster), "Category"]
#     
#   }else{
#     strength_tumour$Novelty <- NA
#   }
#   local_age <- age_enrichment[[tumour]]
#   local_age$Module_ID <- paste(local_age$tumour, local_age$tissue_type, local_age$cluster, sep="_")
#   strength_tumour$Module_ID <- paste(strength_tumour$Tumour, strength_tumour$Tissue_type, strength_tumour$Module, sep="_")
#   
#   strength_tumour$Module_age <- local_age[match(strength_tumour$Module_ID, local_age$Module_ID), "Module_age"]
#   strength_tumour$Module_ID <- NULL
#   
#   if(tumour %in% tumours_normal){
#     strength_normal <- calculate_strength_of_subnetworks(subnet_normal, tumour, "normal")
#     strength_normal$Module_age <- normal_preservation[match(strength_normal$Module, normal_preservation$Cluster), "Age"]
#     strength_normal$Novelty <- normal_preservation[match(strength_normal$Module, normal_preservation$Cluster), "Category"]
#     strength_subnetworks <- data.frame(rbind(strength_tumour, strength_normal))
#   }else{
#     strength_subnetworks <- data.frame(rbind(strength_tumour))
#   }
#   
#   
#   strength_subnetworks$Average_strength <- as.numeric(as.character(strength_subnetworks$Average_strength))
#   
#   strength_subnetworks$Connection_type <- factor(strength_subnetworks$Connection_type,
#                                                  levels=c("All", "UC", "Mixed", "MC"))
#   
#   degree_all <- strength_subnetworks[strength_subnetworks$Connection_type == "All",c("Tumour", "Tissue_type",
#                                                                                      "Module", "Strength",
#                                                                                      "Module_age", "Novelty")]
#   genes_in_tumour_modules <- sapply(cluster_assignments$tumour[[tumour]], length)
#   
#   degree_tumour <- degree_all[degree_all$Tissue_type == "tumour",]
#   degree_tumour$Number_genes <- genes_in_tumour_modules[match(degree_tumour$Module, names(genes_in_tumour_modules))]
#   degree_tumour$Overall_degree <- degree_tumour$Strength/degree_tumour$Number_genes
#   gene_attributes_tumour <- get_metrics_subnetworks(subnet_tumour, "tumour", tumour)
#   
#   
#   if(tumour %in% tumours_normal){
#     genes_in_normal_modules <- sapply(cluster_assignments$normal[[tumour]], length)
#     degree_normal <- degree_all[degree_all$Tissue_type == "normal",]
#     degree_normal$Number_genes <- genes_in_normal_modules[match(degree_normal$Module, names(genes_in_normal_modules))]
#     degree_normal$Overall_degree <- degree_normal$Strength/degree_normal$Number_genes
#     overall_degree <- rbind(degree_tumour, degree_normal)
#     gene_attributes_normal <- get_metrics_subnetworks(subnet_normal, "normal", tumour)
#     gene_attributes <- rbind(gene_attributes_tumour, gene_attributes_normal)
#     
#   }else{
#     overall_degree <- rbind(degree_tumour)
#     gene_attributes <- rbind(gene_attributes_tumour)
#     
#   }
#   
#   gene_attributes$CNV[is.na(gene_attributes$CNV)] <- "No_CNV"
#   gene_attributes$Mut[is.na(gene_attributes$Mut)] <- "No_Mut"
#   
#   attributes_age_mean <- aggregate(Degree ~ Age+Module+Tissue_type, gene_attributes, mean)
#   attributes_CNV_mean <- aggregate(Degree ~ CNV+Module+Tissue_type, gene_attributes, mean)
#   attributes_mut_mean <- aggregate(Degree ~ Mut+Module+Tissue_type, gene_attributes, mean)
#   
#   attributes_age_mean$Degree_reference <- overall_degree[match(interaction(attributes_age_mean$Tissue_type, attributes_age_mean$Module),
#                                                                interaction(overall_degree$Tissue_type, overall_degree$Module)), "Overall_degree"]
#   attributes_age_mean$Degree_normalized <- attributes_age_mean$Degree/attributes_age_mean$Degree_reference
#   
#   attributes_CNV_mean$Degree_reference <- overall_degree[match(interaction(attributes_CNV_mean$Tissue_type, attributes_CNV_mean$Module),
#                                                                interaction(overall_degree$Tissue_type, overall_degree$Module)), "Overall_degree"]
#   attributes_CNV_mean$Degree_normalized <- attributes_CNV_mean$Degree/attributes_CNV_mean$Degree_reference
#   
#   attributes_mut_mean$Degree_reference <- overall_degree[match(interaction(attributes_mut_mean$Tissue_type, attributes_mut_mean$Module),
#                                                                interaction(overall_degree$Tissue_type, overall_degree$Module)), "Overall_degree"]
#   attributes_mut_mean$Degree_normalized <- attributes_mut_mean$Degree/attributes_mut_mean$Degree_reference
#   
#   attributes_age_mean$Tumour <- tumour
#   attributes_CNV_mean$Tumour <- tumour
#   attributes_mut_mean$Tumour <- tumour
#   
#   attributes_age_all <- rbind(attributes_age_all, attributes_age_mean)
#   attributes_CNV_all <- rbind(attributes_CNV_all, attributes_CNV_mean)
#   attributes_mut_all <- rbind(attributes_mut_all, attributes_mut_mean)
#   
#   age_degree_results <- vector()
#   for(tissue_type in c("normal", "tumour")){
#     if(tissue_type == "tumour" || (tissue_type == "normal" & tumour %in% tumours_normal)){
#       local_age_mean <- attributes_age_mean[attributes_age_mean$Tissue_type == tissue_type,]
#       
#       local_age_mean_UC <- local_age_mean[local_age_mean$Age == "UC",]
#       local_age_mean_EM <- local_age_mean[local_age_mean$Age == "EM",]
#       local_age_mean_MM <- local_age_mean[local_age_mean$Age == "MM",]
#       
#       p1 <- wilcox.test(local_age_mean_UC$Degree_normalized, local_age_mean_EM$Degree_normalized, alternative="greater")$p.val
#       p2 <- wilcox.test(local_age_mean_EM$Degree_normalized, local_age_mean_MM$Degree_normalized, alternative="greater")$p.val
#       
#       age_degree_results <- rbind(age_degree_results, c(tissue_type, p1, p2))
#     }
#   }
#   colnames(age_degree_results) <- c("Tissue_type", "UC_greater_EM", "EM_greater_MM")
#   
#   age_degree_results <- as.data.frame(age_degree_results)
#   age_degree_results$Tumour <- tumour
#   age_degree_p_all <- rbind(age_degree_p_all, age_degree_results)
#   
#   
#   CNV_degree_results <- vector()
#   for(tissue_type in c("normal", "tumour")){
#     if(tissue_type == "tumour" || (tissue_type == "normal" & tumour %in% tumours_normal)){
#       local_CNV_mean <- attributes_CNV_mean[attributes_CNV_mean$Tissue_type == tissue_type,]
#       
#       local_CNV_mean_Amp <- local_CNV_mean[local_CNV_mean$CNV == "Amp",]
#       local_CNV_mean_Del <- local_CNV_mean[local_CNV_mean$CNV == "Del",]
#       local_CNV_mean_None <- local_CNV_mean[local_CNV_mean$CNV == "No_CNV",]
#       
#       if(nrow(local_CNV_mean_Amp) > 0){
#         p1 <- wilcox.test(local_CNV_mean_Amp$Degree_normalized, local_CNV_mean_None$Degree_normalized, alternative="less")$p.val
#         p2 <- wilcox.test(local_CNV_mean_Amp$Degree_normalized, local_CNV_mean_None$Degree_normalized, alternative="greater")$p.val
#       }else{
#         p1 <- NA
#         p2 <- NA
#       }
#       if(nrow(local_CNV_mean_Del) > 0){
#         p3 <- wilcox.test(local_CNV_mean_Del$Degree_normalized, local_CNV_mean_None$Degree_normalized, alternative="less")$p.val
#       }else{
#         p3 <- NA
#       }
#       CNV_degree_results <- rbind(CNV_degree_results, c(tissue_type, p1, p2, p3))
#     }
#   }
#   colnames(CNV_degree_results) <- c("Tissue_type", "Amp_less_None", "Amp_greater_None", "Del_less_None")
#   
#   CNV_degree_results <- as.data.frame(CNV_degree_results)
#   CNV_degree_results$Tumour <- tumour
#   CNV_degree_p_all <- rbind(CNV_degree_p_all, CNV_degree_results)
#   
#   mut_degree_results <- vector()
#   for(tissue_type in c("normal", "tumour")){
#     if(tissue_type == "tumour" || (tissue_type == "normal" & tumour %in% tumours_normal)){
#       local_mut_mean <- attributes_mut_mean[attributes_mut_mean$Tissue_type == tissue_type,]
#       
#       local_mut_mean_LoF <- local_mut_mean[local_mut_mean$Mut == "LoF",]
#       local_mut_mean_Miss <- local_mut_mean[local_mut_mean$Mut == "Miss",]
#       local_mut_mean_None <- local_mut_mean[local_mut_mean$Mut == "No_Mut",]
#       
#       if(nrow(local_mut_mean_LoF) > 0){
#         p1 <- wilcox.test(local_mut_mean_LoF$Degree_normalized, local_mut_mean_None$Degree_normalized, alternative="less")$p.val
#       }else{
#         p1 <- NA
#       }
#       
#       if(nrow(local_mut_mean_Miss) > 0){
#         p2 <- wilcox.test(local_mut_mean_Miss$Degree_normalized, local_mut_mean_None$Degree_normalized, alternative="less")$p.val
#       }else{
#         p2 <- NA
#       }
#       mut_degree_results <- rbind(mut_degree_results, c(tissue_type, p1, p2))
#     }
#   }
#   colnames(mut_degree_results) <- c("Tissue_type", "LoF_less_None", "Miss_less_None")
#   
#   
#   mut_degree_results <- as.data.frame(mut_degree_results)
#   mut_degree_results$Tumour <- tumour
#   mut_degree_p_all <- rbind(mut_degree_p_all, mut_degree_results)
#   
#   degree_nodes_tumour <- calculate_degree_rank(subnet_tumour, "Tumour", tumour)
#   colnames(degree_nodes_tumour)[2] <- "Degree_tumour"
#   colnames(degree_nodes_tumour)[3] <- "Degree_norm_tumour"
#   colnames(degree_nodes_tumour)[4] <- "Degree_rank_tumour"
#   colnames(degree_nodes_tumour)[5] <- "Degree_rank_norm_tumour"
#   colnames(degree_nodes_tumour)[6] <- "Module_tumour"
#   
#   
#   if(tumour %in% tumours_normal){
#     degree_nodes_normal <- calculate_degree_rank(subnet_normal, "Normal", tumour)
#     colnames(degree_nodes_normal)[2] <- "Degree_normal"
#     colnames(degree_nodes_normal)[3] <- "Degree_norm_normal"
#     colnames(degree_nodes_normal)[4] <- "Degree_rank_normal"
#     colnames(degree_nodes_normal)[5] <- "Degree_rank_norm_normal"
#     colnames(degree_nodes_normal)[6] <- "Module_normal"
#     degree_nodes <- data.frame(degree_nodes_tumour[,c(1,2,3,4,5,6)],
#                                degree_nodes_normal[match(degree_nodes_tumour$Gene,
#                                                          degree_nodes_normal$Gene),c(2,3,4,5,6)])
#     
#   }else{
#     degree_nodes <- data.frame(degree_nodes_tumour[,c(1,2,3,4,5,6)])
#     
#   }
#   
#   
#   
#   
#   degree_nodes$CNV <- ""
#   degree_nodes$CNV[degree_nodes$Gene %in% genes_amp] <- "Amp"
#   degree_nodes$CNV[degree_nodes$Gene %in% genes_del] <- "Del"
#   
#   degree_nodes$Mut <- ""
#   degree_nodes$Mut[degree_nodes$Gene %in% genes_miss] <- "Miss"
#   degree_nodes$Mut[degree_nodes$Gene %in% genes_lof] <- "LoF"
#   degree_nodes$Alt <- paste(degree_nodes$CNV, degree_nodes$Mut, sep="")
#   degree_nodes$Alt <- trimws(degree_nodes$Alt, which = "both")
#   degree_nodes$Alt[degree_nodes$Alt == ""] <- NA
#   
#   degree_nodes_temp <- degree_nodes
#   degree_nodes_temp$Alt[is.na(degree_nodes_temp$Alt)] <- "None"
#   
#   degree_nodes_temp$Tumour <- tumour
#   degree_nodes_all <- rbind(degree_nodes_all, degree_nodes_temp)
#   
#   print(tumour)
#   
#   save(degree_nodes_all, file=paste("degree_nodes_all_",
#                                     tumour, ".Rdata", sep=""))
#   save(attributes_age_all, file=paste("attributes_age_all_",
#                                       tumour, ".Rdata", sep=""))
#   save(attributes_CNV_all, file=paste("attributes_CNV_all_",
#                                       tumour, ".Rdata", sep=""))
#   save(attributes_mut_all, file=paste("attributes_mut_all_",
#                                       tumour, ".Rdata", sep=""))
#   save(age_degree_p_all, file=paste("age_degree_p_all_",
#                                     tumour, ".Rdata", sep=""))
#   save(CNV_degree_p_all, file=paste("CNV_degree_p_all_",
#                                     tumour, ".Rdata", sep=""))
#   save(mut_degree_p_all, file=paste("mut_degree_p_all_",
#                                     tumour, ".Rdata", sep=""))
# }

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
  
  
  load(paste("attributes_age_all_",
             tumour, ".Rdata", sep=""))
  attributes_age_all2 <- rbind(attributes_age_all2, attributes_age_all)
  
  load(paste("attributes_CNV_all_",
             tumour, ".Rdata", sep=""))
  attributes_CNV_all2 <- rbind(attributes_CNV_all2, attributes_CNV_all)
  
  load(paste("attributes_mut_all_",
             tumour, ".Rdata", sep=""))
  attributes_mut_all2 <- rbind(attributes_mut_all2, attributes_mut_all)
  
  load(paste("age_degree_p_all_",
             tumour, ".Rdata", sep=""))
  age_degree_p_all2 <- rbind(age_degree_p_all2, age_degree_p_all)
  
  load(paste("CNV_degree_p_all_",
             tumour, ".Rdata", sep=""))
  CNV_degree_p_all2 <- rbind(CNV_degree_p_all2, CNV_degree_p_all)
  
  load(paste("mut_degree_p_all_",
             tumour, ".Rdata", sep=""))
  mut_degree_p_all2 <- rbind(mut_degree_p_all2, mut_degree_p_all)
}

degree_nodes_all <- degree_nodes_all2
attributes_age_all <- attributes_age_all2
attributes_CNV_all <- attributes_CNV_all2
attributes_mut_all <- attributes_mut_all2

age_degree_p_all <- age_degree_p_all2
CNV_degree_p_all <- CNV_degree_p_all2
mut_degree_p_all <- mut_degree_p_all2

attributes_CNV_all$Tumour <- factor(attributes_CNV_all$Tumour,
                                    levels=tumours)

attributes_CNV_all$CNV[attributes_CNV_all$CNV == "No_CNV"] <- "No CNA"
attributes_CNV_all$Tissue_type[attributes_CNV_all$Tissue_type == "normal"] <- "Normal"
attributes_CNV_all$Tissue_type[attributes_CNV_all$Tissue_type == "tumour"] <- "Tumour"

degree_nodes_all$Degree_rank_norm_normal <- 1-degree_nodes_all$Degree_rank_norm_normal
degree_nodes_all$Degree_rank_norm_tumour <- 1-degree_nodes_all$Degree_rank_norm_tumour

degree_nodes_all$Tumour <- factor(degree_nodes_all$Tumour, levels=tumours)

degree_nodes_all$Alt <- factor(degree_nodes_all$Alt,
                               levels=c("None", "Amp", "Del", "Miss", "LoF"))


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

degree_nodes_amp <- subset(degree_nodes_all, Alt %in% c("Amp")
                           & !(Tumour %in% c("ACC", "COAD", "KIRP")))##No amplifications

degree_nodes_amp <- degree_nodes_amp[,c("Tumour", "Alt", "Degree_rank_norm_tumour", "Degree_rank_norm_normal","Gene")]
degree_nodes_amp <- melt(degree_nodes_amp)

degree_nodes_amp$variable <- as.character(degree_nodes_amp$variable)
degree_nodes_amp$variable[degree_nodes_amp$variable == "Degree_rank_norm_tumour"] <- "Tumour"
degree_nodes_amp$variable[degree_nodes_amp$variable == "Degree_rank_norm_normal"] <- "Normal"

degree_nodes_amp_median <- aggregate(value ~ Tumour+Alt+variable, degree_nodes_amp, median)

##Panel C
pdf("Figure_3C.pdf", height=3.85, width=5)
g <- ggplot(degree_nodes_amp_median, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(colour=Tumour), width=0.1)+
  ylab("Median normalized degree")+
  xlab("Tissue")+
  ggtitle("Amp")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

wilcox.test(degree_nodes_amp_median$value[degree_nodes_amp_median$variable == "Tumour"],
            degree_nodes_amp_median$value[degree_nodes_amp_median$variable == "Normal"], alternative="greater")



degree_nodes_del <- subset(degree_nodes_all, Alt %in% c("Del")
                           & !(Tumour %in% c("UVM")))

degree_nodes_del <- degree_nodes_del[,c("Tumour", "Alt", "Degree_rank_norm_tumour", "Degree_rank_norm_normal","Gene")]
degree_nodes_del <- melt(degree_nodes_del)

degree_nodes_del$variable <- as.character(degree_nodes_del$variable)
degree_nodes_del$variable[degree_nodes_del$variable == "Degree_rank_norm_tumour"] <- "Tumour"
degree_nodes_del$variable[degree_nodes_del$variable == "Degree_rank_norm_normal"] <- "Normal"

degree_nodes_del_median <- aggregate(value ~ Tumour+Alt+variable, degree_nodes_del, median)

pdf("Figure_S11.pdf", height=4, width=5.5)
g <- ggplot(degree_nodes_del_median, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(colour=Tumour), width=0.1)+
  ylab("Median normalized degree")+
  xlab("Tissue")+
  ggtitle("Del")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

degree_nodes_all$Module_name <- paste(degree_nodes_all$Tumour, degree_nodes_all$Module_tumour, sep="_")

all_preservation_t_to_n2$Module_name <- paste(all_preservation_t_to_n2$Tumour, all_preservation_t_to_n2$Cluster, sep="_")

degree_nodes_all$Module_age <- all_preservation_t_to_n2[match(degree_nodes_all$Module_name, all_preservation_t_to_n2$Module_name), "Age"]
degree_nodes_all$Novelty <- all_preservation_t_to_n2[match(degree_nodes_all$Module_name, all_preservation_t_to_n2$Module_name), "Category"]
degree_nodes_all$Diff <- degree_nodes_all$Degree_norm_tumour-degree_nodes_all$Degree_norm_normal

key_genes <- degree_nodes_all
key_genes$Change_centrality <- key_genes$Degree_rank_norm_tumour-key_genes$Degree_rank_norm_normal

key_genes$Change_centrality_abs <- abs(key_genes$Change_centrality)

key_genes$Module_name <- paste(key_genes$Tumour, key_genes$Module_tumour, sep="_")

all_preservation_t_to_n2$Module_name <- paste(all_preservation_t_to_n2$Tumour, all_preservation_t_to_n2$Cluster, sep="_")

key_genes$Module_age <- all_preservation_t_to_n2[match(key_genes$Module_name, all_preservation_t_to_n2$Module_name), "Age"]
key_genes$Novelty <- all_preservation_t_to_n2[match(key_genes$Module_name, all_preservation_t_to_n2$Module_name), "Category"]

key_genes$Alt <- paste(key_genes$CNV, key_genes$Mut, sep=",")


cancer_genes <- read.csv("Census_allWed Dec 21 23_46_30 2022.csv")
cancer_genes$Gene.Symbol <- toupper(cancer_genes$Gene.Symbol)


cancer_genes$Tumour.Types.Somatic. <- as.character(cancer_genes$Tumour.Types.Somatic.)

##Figure S12
for(tumour in tumours_normal){
  temp <- key_genes[key_genes$Tumour == tumour,]
  
  temp$Change_centrality <- as.numeric(as.character(temp$Change_centrality))
  temp$Change_centrality_abs <- abs(temp$Change_centrality)
  temp <- temp[temp$Gene %in% cancer_genes_per_tumour[[tumour]],]
  
  temp$Alt <- as.character(temp$Alt)
  temp$Alt[is.na(temp$Alt)] <- paste(temp$CNV[is.na(temp$Alt)], temp$Mut[is.na(temp$Alt)], sep=",")
  
  temp <- subset(temp, Module_age == "Mixed")
  temp$Alt[temp$Alt == ","] <- "None"
  temp$Alt[temp$Alt == ",LoF"] <- "LoF"
  temp$Alt[temp$Alt == ",Miss"] <- "Miss"
  temp$Alt[temp$Alt == "Amp,"] <- "Amp"
  
  pdf(paste0("Figure_S12_", tumour, ".pdf"), height=4, width=5)
  g <- ggplot(temp, aes(x=Degree_rank_norm_normal, y=Degree_rank_norm_tumour))+
    geom_point(aes(colour=Alt), size=3)+
    geom_abline(intercept = 0, slope=1)+
    ggtitle(tumour)+
    geom_text_repel(aes(label=Gene), size=3)+
    xlab("Centrality in normal modules")+
    ylab("Centrality in tumour modules")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(g)
  dev.off()
}

fisher.test(cbind(c(2,337-2),c(51,1166-51)), alt="l")


fisher.test(cbind(c(1,337-1),c(39,1166-39)), alt="l")
