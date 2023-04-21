#Script for Figure S7.

library(ggplot2)
library(gridExtra)
library("ggsci")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")


load("all_preservation_t_to_n2.Rdata")

all_preservation_t_to_n2$Category <- as.character(all_preservation_t_to_n2$Category)
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Low_score"] <- "Low novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Median_score"] <- "Moderate novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "High_score"] <- "High novelty"
all_preservation_t_to_n2$Category[all_preservation_t_to_n2$Category == "Inf_score"] <- "High novelty"
all_preservation_t_to_n2$Category <- factor(all_preservation_t_to_n2$Category,
                                            levels=c("Low novelty", "Moderate novelty", "High novelty"))



##Fraction of modules by novelty
fraction_novelty <- vector()
for(tumour in tumours_normal){
  temp1 <- all_preservation_t_to_n2[all_preservation_t_to_n2$Tumour == tumour,]
  for (novelty in c("High novelty", "Moderate novelty", "Low novelty")){
    temp2 <- temp1[temp1$Category == novelty,]
    values <- as.vector(table(temp2$Age))
    values <- values/sum(values)
    fraction_novelty <- rbind(fraction_novelty,
                              c(tumour, novelty,values))
    
  }
}


fraction_novelty <- as.data.frame(fraction_novelty)
colnames(fraction_novelty) <- c("Tumour", "Novelty", "UC", "Mixed", "MC")
fraction_novelty[,3] <- as.numeric(as.character(fraction_novelty[,3]))
fraction_novelty[,4] <- as.numeric(as.character(fraction_novelty[,4]))
fraction_novelty[,5] <- as.numeric(as.character(fraction_novelty[,5]))

library(reshape2)
fraction_novelty_melt <- melt(fraction_novelty)
colnames(fraction_novelty_melt)[3:4] <- c("Age", "Fraction")

fraction_novelty_melt$Novelty <- factor(fraction_novelty_melt$Novelty,
                                        levels=c("Low novelty", "Moderate novelty", "High novelty"))

fraction_novelty_melt$Age <- factor(fraction_novelty_melt$Age,
                                    levels=c("UC", "Mixed", "MC"))

fraction_novelty_melt$Tumour <- factor(fraction_novelty_melt$Tumour,
                                       levels=tumours)

pdf("Figure_S7.pdf", height=4, width=12)
g <- ggplot(fraction_novelty_melt, aes(x=Tumour, y = Fraction))+
  geom_bar(stat="identity", aes(fill=Age))+
  ylab("Fraction of modules")+
  facet_grid(.~Novelty)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g)
dev.off()