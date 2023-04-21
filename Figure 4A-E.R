##Script for Figure 4, panels A-E.

library(ggplot2)

##Panels A-C
source("functions.R")
load("age_enrichment_PRAD_2.Rdata")
load("modules_PRAD_2.Rdata")
load("preservation_PRAD_2.Rdata")

all_preservation_PRAD <- preservation[,c("Type", "Cluster", "Per_50", "Label", "Size")]

#Divide preservation score by module size
all_preservation_PRAD$Preservation_ratio <- all_preservation_PRAD$Per_50/all_preservation_PRAD$Size

##Age
age_enrichment_df <- vector()
for(type in c("normal", "low_grade", "high_grade")){
  local_age_enrichment <- age_enrichment[[type]]
  age_enrichment_df <- rbind(age_enrichment_df, local_age_enrichment)
}
age_enrichment_df$Label <- paste(age_enrichment_df$type, age_enrichment_df$cluster)

all_preservation_PRAD$Age <- age_enrichment_df$Module_age[match(all_preservation_PRAD$Label, 
                                                                age_enrichment_df$Label)]
all_preservation_PRAD$Age <- factor(all_preservation_PRAD$Age, levels=c("UC", "Mixed", "MC"))
all_preservation_PRAD$Type <- as.character(all_preservation_PRAD$Type)

all_preservation_PRAD$Type[all_preservation_PRAD$Type == "high_grade"] <- "low_to_high_grade"
all_preservation_PRAD$Type[all_preservation_PRAD$Type == "low_grade"] <- "normal_to_low_grade"
all_preservation_PRAD$Type <- factor(all_preservation_PRAD$Type, 
                                     levels=c("normal_to_low_grade", "low_to_high_grade"))


library(ggplot2)
#Panel A
pdf("Figure_4A.pdf", height=2.75, width=5)
g <- ggplot(all_preservation_PRAD, aes(x=Age, y=Preservation_ratio))+
  geom_boxplot(aes(fill=Age))+
  facet_grid(.~Type)+
  ylab("Novelty")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()


##Panel B

all_preservation_PRAD_o <- all_preservation_PRAD

all_preservation_PRAD <- define_categories_preservation_pheo(all_preservation_PRAD)

all_preservation_PRAD$Age <- factor(all_preservation_PRAD$Age,
                                    levels=c("UC", "Mixed", "MC"))
preservation_cat_PRAD <- as.data.frame(table(all_preservation_PRAD[,c("Type", "Age", "Category")]))


preservation_cat_PRAD2 <- preservation_cat_PRAD
preservation_cat_PRAD2 <- preservation_cat_PRAD2[preservation_cat_PRAD2$Category != "Inf_score",]

library("ggsci")
preservation_cat_PRAD2$Category <- factor(preservation_cat_PRAD2$Category,
                                          levels=c("High_score", "Medium_score", "Low_score"))

pdf("Figure_4B.pdf", height=2.75, width=5)
g <- ggplot(preservation_cat_PRAD2, aes(x=Category, y=Freq))+
  geom_bar(stat='identity', aes(fill=Category))+
  scale_fill_jco()+
  xlab("Module novelty")+ ylab("Number of modules")+
  facet_grid(Type ~ Age, scales="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()


#Panel C
total_mod_low_grade <- sum(preservation_cat_PRAD[preservation_cat_PRAD$Type == "normal_to_low_grade", "Freq"])
total_mod_high_grade <- sum(preservation_cat_PRAD[preservation_cat_PRAD$Type == "low_to_high_grade", "Freq"])

preservation_cat_PRAD$Total <- NA
preservation_cat_PRAD$Total[preservation_cat_PRAD$Type == "normal_to_low_grade"] <- total_mod_low_grade
preservation_cat_PRAD$Total[preservation_cat_PRAD$Type == "low_to_high_grade"] <- total_mod_high_grade
preservation_cat_PRAD$Percentage <- (preservation_cat_PRAD$Freq/preservation_cat_PRAD$Total)*100

preservation_cat_PRAD <- preservation_cat_PRAD[preservation_cat_PRAD$Category != "Inf_score",]


pdf("Figure_4C.pdf", height=2.75, width=6)
g <- ggplot(preservation_cat_PRAD, aes(x=Category, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Type), position="dodge")+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  ylab("Percentage of\n normal to low and\n low to high modules")+
  facet_grid(. ~ Age)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()


##Panel D
load("age_enrichment_melanoma.Rdata")
load("preservation_melanoma.Rdata")

age_enrichment_df <- vector()
for(type in c("nevus", "primary")){
  local_age_enrichment <- age_enrichment[[type]]
  age_enrichment_df <- rbind(age_enrichment_df, local_age_enrichment)
}
age_enrichment_df$Label <- paste(age_enrichment_df$type, age_enrichment_df$cluster)

all_preservation_melanoma <- preservation_melanoma[,c("Type", "Cluster", "Per_50", "Label", "Size")]
all_preservation_melanoma$Preservation_ratio <- all_preservation_melanoma$Per_50/all_preservation_melanoma$Size

all_preservation_melanoma$Age <- age_enrichment_df$Module_age[match(all_preservation_melanoma$Label, 
                                                                    age_enrichment_df$Label)]
all_preservation_melanoma$Age <- factor(all_preservation_melanoma$Age, levels=c("UC", "Mixed", "MC"))

pdf("Figure_4D.pdf", height=2.75, width=4)
g <- ggplot(all_preservation_melanoma, aes(x=Age, y=Preservation_ratio))+
  geom_boxplot(aes(fill=Age))+
  ylab("Novelty")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age == "Mixed"],
            all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age != "Mixed"])

wilcox.test(all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age == "Mixed"],
            all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age == "UC"], alternative="greater")

wilcox.test(all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age == "Mixed"],
            all_preservation_melanoma$Preservation_ratio[all_preservation_melanoma$Age == "MC"], alternative="greater")


##Panel E - MAML
source("functions.R")
load("all_preservation_pheo_MAML.Rdata")

all_preservation_pheo <- define_categories_preservation_pheo(all_preservation_pheo)

all_preservation_pheo$Age <- factor(all_preservation_pheo$Age,
                                    levels=c("UC", "Mixed", "MC"))
preservation_cat_pheo <- as.data.frame(table(all_preservation_pheo[,c("Type", "Age", "Category")]))

total_mod_benign <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "benign", "Freq"])
total_mod_malignant <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "malignant", "Freq"])

preservation_cat_pheo$Total <- NA
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "benign"] <- total_mod_benign
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "malignant"] <- total_mod_malignant
preservation_cat_pheo$Percentage <- (preservation_cat_pheo$Freq/preservation_cat_pheo$Total)*100

preservation_cat_pheo <- preservation_cat_pheo[preservation_cat_pheo$Category != "Inf_score",]

pdf("Figure_4E_MAML.pdf", height=2.5, width=4)
g <- ggplot(preservation_cat_pheo, aes(x=Category, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Type), position="dodge")+
  facet_grid(. ~ Age)+
  ggtitle("MAML")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "UC"], alternative="greater")

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "MC"], alternative="greater")


##Panel E - SDH 
source("functions.R")
load("all_preservation_pheo_SDH.Rdata")

all_preservation_pheo <- define_categories_preservation_pheo(all_preservation_pheo)

all_preservation_pheo$Age <- factor(all_preservation_pheo$Age,
                                    levels=c("UC", "Mixed", "MC"))
preservation_cat_pheo <- as.data.frame(table(all_preservation_pheo[,c("Type", "Age", "Category")]))

total_mod_benign <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "benign", "Freq"])
total_mod_malignant <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "malignant", "Freq"])

preservation_cat_pheo$Total <- NA
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "benign"] <- total_mod_benign
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "malignant"] <- total_mod_malignant
preservation_cat_pheo$Percentage <- (preservation_cat_pheo$Freq/preservation_cat_pheo$Total)*100

preservation_cat_pheo <- preservation_cat_pheo[preservation_cat_pheo$Category != "Inf_score",]

pdf("Figure_4E_SDH.pdf", height=2.5, width=4)
g <- ggplot(preservation_cat_pheo, aes(x=Category, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Type), position="dodge")+
  facet_grid(. ~ Age)+
  ggtitle("SDH")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "UC"], alternative="greater")

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "MC"], alternative="greater")


#Panel E - VHL
source("functions.R")
load("all_preservation_pheo_VHL.Rdata")

all_preservation_pheo <- define_categories_preservation_pheo(all_preservation_pheo)

all_preservation_pheo$Age <- factor(all_preservation_pheo$Age,
                                    levels=c("UC", "Mixed", "MC"))
preservation_cat_pheo <- as.data.frame(table(all_preservation_pheo[,c("Type", "Age", "Category")]))

total_mod_benign <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "benign", "Freq"])
total_mod_malignant <- sum(preservation_cat_pheo[preservation_cat_pheo$Type == "malignant", "Freq"])

preservation_cat_pheo$Total <- NA
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "benign"] <- total_mod_benign
preservation_cat_pheo$Total[preservation_cat_pheo$Type == "malignant"] <- total_mod_malignant
preservation_cat_pheo$Percentage <- (preservation_cat_pheo$Freq/preservation_cat_pheo$Total)*100

preservation_cat_pheo <- preservation_cat_pheo[preservation_cat_pheo$Category != "Inf_score",]


pdf("Figure_4E_VHL.pdf", height=2.5, width=4)
g <- ggplot(preservation_cat_pheo, aes(x=Category, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Type), position="dodge")+
  facet_grid(. ~ Age)+
  ggtitle("VHL")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "UC"], alternative="greater")

wilcox.test(all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "Mixed"],
            all_preservation_pheo$Preservation_ratio[all_preservation_pheo$Age == "MC"], alternative="greater")
