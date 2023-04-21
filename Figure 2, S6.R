##Script for all panels of Figure 2 and Figure S6.

##Panels B and C
library(ggplot2)
library(gridExtra)
library(ggsci)
library(reshape2)
library(ggraph)
library(igraph)
library(ggrepel)


tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

load("all_preservation_t_to_n2.Rdata")

##Figure 2C
pdf("Figure_2C.pdf", height=3, width=4.5)
g <- ggplot(all_preservation_t_to_n2, aes(x=Age, y =Preservation_ratio))+
  geom_boxplot(aes(fill=Age))+
  ylab("Novelty")+
  xlab("Module age")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(all_preservation_t_to_n2$Preservation_ratio[all_preservation_t_to_n2$Age == "Mixed"],
            all_preservation_t_to_n2$Preservation_ratio[all_preservation_t_to_n2$Age == "UC"],
            alternative="greater")

wilcox.test(all_preservation_t_to_n2$Preservation_ratio[all_preservation_t_to_n2$Age == "Mixed"],
            all_preservation_t_to_n2$Preservation_ratio[all_preservation_t_to_n2$Age == "MC"],
            alternative="greater")

##Supplementary Figure S6
pdf("Figure_S6.pdf", height=2, width=10)
g <- ggplot(all_preservation_t_to_n2, aes(x=Age, y =log10(Preservation_ratio)))+
  geom_boxplot(aes(fill=Age))+
  facet_grid(.~Tumour)+
  ylab("Novelty (log10)")+
  xlab("Module age")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())
print(g)
dev.off()

all_preservation_t_to_n2$Category <- as.character(all_preservation_t_to_n2$Category)
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Low_score"] <- "Low novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Median_score"] <- "Moderate novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "High_score"] <- "High novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Inf_score"] <- "High novelty"
all_preservation_t_to_n2$Category <- factor(all_preservation_t_to_n2$Category,
                                            levels=c("Low novelty", "Moderate novelty", "High novelty"))

fraction_age_novelty2 <- vector()
for(tumour in tumours_normal){
  temp1 <- all_preservation_t_to_n2[all_preservation_t_to_n2$Tumour == tumour,]
  for (age in c("UC", "Mixed", "MC")){
    temp2 <- temp1[temp1$Age == age,]
    values <- as.vector(table(temp2$Category))
    values <- values/sum(values)
    fraction_age_novelty2 <- rbind(fraction_age_novelty2,
                                  c(tumour, age,values))
    
  }
}


fraction_age_novelty2 <- as.data.frame(fraction_age_novelty2)
colnames(fraction_age_novelty2) <- c("Tumour", "Age", "Low novelty", "Moderate novelty", "High novelty")
fraction_age_novelty2[,3] <- as.numeric(as.character(fraction_age_novelty2[,3]))
fraction_age_novelty2[,4] <- as.numeric(as.character(fraction_age_novelty2[,4]))
fraction_age_novelty2[,5] <- as.numeric(as.character(fraction_age_novelty2[,5]))

fraction_age_novelty_melt2 <- melt(fraction_age_novelty2)
colnames(fraction_age_novelty_melt2)[3:4] <- c("Novelty", "Fraction")

fraction_age_novelty_melt2$Novelty <- factor(fraction_age_novelty_melt2$Novelty,
                                            levels=c("Low novelty",
                                                     "Moderate novelty",
                                                     "High novelty"))

fraction_age_novelty_melt2$Age <- factor(fraction_age_novelty_melt2$Age,
                                             levels=c("UC", "Mixed", "MC"))

fraction_age_novelty_melt2$Tumour <- factor(fraction_age_novelty_melt2$Tumour,
                                           levels=tumours)


temp_mean <- aggregate(Fraction ~ Age+Novelty, fraction_age_novelty_melt2, mean)

temp_mean$Novelty <- factor(temp_mean$Novelty, levels = c("High novelty", "Moderate novelty", "Low novelty"))

##Figure 2B
pdf("Figure_2B.pdf", height=3, width=4.25)
g <- ggplot(temp_mean, aes(x=Novelty, y=Fraction))+
  geom_point(aes(colour=Age), size=3)+
  geom_line(aes(group=Age, colour=Age), size=1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()



###Panel D
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

strength_subnetworks <- vector()
for(tumour in tumours){
  for(tissue_type in c("normal", "tumour")){
    if((tissue_type == "tumour") | (tissue_type == "normal" & tumour %in% tumours_normal)){
      load(paste("Subnetworks_", tumour, "_", tissue_type,
                 ".Rdata", sep=""))
      local_subnet <- sub_networks[[tissue_type]][[tumour]]
      
      local_ages <- age_enrichment[[tumour]]
      local_ages <- local_ages[local_ages$tissue_type == tissue_type,]
      
      strength <- calculate_strength_of_subnetworks(local_subnet, tumour, tissue_type)
      strength$Module_age <- local_ages[match(strength$Module, local_ages$cluster), "Module_age"]
      strength_subnetworks <- rbind(strength_subnetworks, strength)
    }
  }
  print(tumour)
}

strength_subnetworks <- as.data.frame(strength_subnetworks)
strength_subnetworks$Average_strength <- as.numeric(as.character(strength_subnetworks$Average_strength))

strength_subnetworks$Connection_type <- factor(strength_subnetworks$Connection_type,
                                               levels=c("All", "UC", "Mixed", "MC"))

strength_subnetworks$Module_age <- factor(strength_subnetworks$Module_age,
                                          levels=c("UC", "Mixed", "MC"))

strength_subnetworks_tumour <- subset(strength_subnetworks, Connection_type != "All" & Tissue_type == "tumour")
strength_subnetworks_tumour$Connection_type2 <- ifelse(strength_subnetworks_tumour$Connection_type == "UC", "UC_UC",
                                                       ifelse(strength_subnetworks_tumour$Connection_type == "MC", "MC_MC",
                                                              ifelse(strength_subnetworks_tumour$Connection_type == "Mixed", "UC_MC", NA)))
strength_subnetworks_tumour$Connection_type2 <- factor(strength_subnetworks_tumour$Connection_type2,
                                                       levels=c("UC_UC", "UC_MC", "MC_MC"))


#Figure 2D
pdf("Figure_2D.pdf", height=3, width=5)
g <- ggplot(strength_subnetworks_tumour, 
       aes(x=Connection_type2, y=Average_strength))+
  geom_boxplot(aes(fill=Connection_type))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(strength_subnetworks_tumour$Average_strength[strength_subnetworks_tumour$Connection_type2 == "UC_MC"],
            strength_subnetworks_tumour$Average_strength[strength_subnetworks_tumour$Connection_type2 == "UC_UC"],
            alternative="less")

wilcox.test(strength_subnetworks_tumour$Average_strength[strength_subnetworks_tumour$Connection_type2 == "UC_MC"],
            strength_subnetworks_tumour$Average_strength[strength_subnetworks_tumour$Connection_type2 == "MC_MC"],
            alternative="less")


##Panel E
load("ssGSEA_normal_exp_tumour_modules.Rdata")
load("age_enrichment.Rdata")
load("all_preservation_t_to_n2.Rdata")


##The raw ssGSEA scores

ssGSEA_scores <- vector()
for(tumour in names(ssGSEA_normal_exp_tumour_modules)){
  tumour_mod <- ssGSEA_normal_exp_tumour_modules[[tumour]]$tumour
  colnames(tumour_mod) <- gsub("T_", "P_", colnames(tumour_mod))
  tumour_mod <- melt(tumour_mod)
  
  tumour_novelty <- all_preservation_t_to_n2[all_preservation_t_to_n2$Tumour == tumour,]
  tumour_mod$Novelty <- tumour_novelty$Category[match(tumour_mod$Var1, tumour_novelty$Cluster)]
  colnames(tumour_mod) <- c("Module", "Patient", "ssGSEA", "Novelty")
  ssGSEA_scores <- rbind(ssGSEA_scores, tumour_mod)
}

ssGSEA_scores <- ssGSEA_scores[!is.na(ssGSEA_scores$Novelty),] 
high_mean <- mean(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "High_score"], na.rm=TRUE)
medium_mean <- mean(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Median_score"], na.rm=TRUE)
low_mean <- mean(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Low_score"], na.rm=TRUE)


SE_low <- sd(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Low_score"])/sqrt(length(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Low_score"]))
SE_medium <- sd(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Median_score"])/sqrt(length(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "Median_score"]))
SE_high <- sd(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "High_score"])/sqrt(length(ssGSEA_scores$ssGSEA[ssGSEA_scores$Novelty == "High_score"]))

df <- data.frame(Novelty = c("Low", "Medium", "High"),
                 ssGSEA = c(low_mean, medium_mean, high_mean),
                 SE = c(SE_low, SE_medium, SE_high))
df$Novelty <- factor(df$Novelty, 
                     levels=c("High", "Medium", "Low"))

#Figure 2E
pdf("Figure_2E.pdf", height=3, width=4)
g <- ggplot(df, aes(x=Novelty, y=ssGSEA)) +
  geom_bar(stat="identity", aes(fill=Novelty)) +
  geom_hline(yintercept=0)+
  scale_fill_jco()+
  geom_errorbar(aes(x=Novelty, ymin=ssGSEA-SE, ymax=ssGSEA+SE),
                width=0.2,linewidth=1)+
  ylab("ssGSEA scores")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()
