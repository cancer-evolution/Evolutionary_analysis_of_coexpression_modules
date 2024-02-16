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

load("Objects/all_preservation_t_to_n2_high.Rdata")

##Figure 2C
pdf("Figure_2C_high.pdf", height=3, width=4.5)
g <- ggplot(all_preservation_t_to_n2_high, aes(x=Age, y =Preservation_ratio))+
  geom_boxplot(aes(fill=Age))+
  ylab("Novelty")+
  xlab("Module age")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()
