##Script for Figures S5 with random modules

library(ggplot2)
library(reshape2)
library(knitr)
library(tidyverse)
library(data.table)

load("Objects/age_enrichment.Rdata")

##### Aggregate  module age data from all cohorts into one data.frame for plotting

### Get the counts of each module type in each age class per cohort
mod.age.agg.list <- lapply( age_enrichment, FUN=function(X) { Y=filter(X, cluster!="grey"); aggregate (Y$Module_age, by=list(Y$tissue_type, Y$Module_age), FUN=length ) } )

### Turn into df with cohort names attached
for( i in 1:length( mod.age.agg.list) ) {
  mod.ages <- mod.age.agg.list[[i]]
  num.mods <- nrow( mod.ages )
  cohort = names(mod.age.agg.list)[i]  
  df <- cbind(mod.ages, "Cohort"=rep(cohort,num.mods) )
  mod.age.agg.list[[cohort]] <- df
}

mod.age.agg.df <- do.call(rbind, mod.age.agg.list)
names(mod.age.agg.df) <- c("Source","Age Enrichment Category","Count","Cohort")

##### Add labels and make plot
mod.age.agg.df$`Age Enrichment Category` <- as.factor(mod.age.agg.df$Age)

mod.age.agg.df$Source <- str_replace(mod.age.agg.df$Source, "tum", "Tum")
mod.age.agg.df$Source <- str_replace(mod.age.agg.df$Source, "norm", "Norm")


load("Objects/boot_tumour_result.Rdata")
load("Objects/boot_normal_result.Rdata")

boot_both_result <- rbind(boot_tumour_result, boot_normal_result)
boot_both_result <- as.data.frame(boot_both_result)
boot_both_result$i <- as.numeric(as.character(boot_both_result$i))
boot_both_result$UC_n <- as.numeric(as.character(boot_both_result$UC_n))
boot_both_result$MC_n <- as.numeric(as.character(boot_both_result$MC_n))
boot_both_result$Mixed_n <- as.numeric(as.character(boot_both_result$Mixed_n))
boot_both_result$Total <- as.numeric(as.character(boot_both_result$Total))

#Average across iterations
boot_both_result_agg_UC <- aggregate(data=boot_both_result, UC_n ~ tumour+tissue, mean)
boot_both_result_agg_UC$Age <- "UC"
colnames(boot_both_result_agg_UC) <- c("Cohort", "Source", "Count", "Age")
boot_both_result_agg_MC <- aggregate(data=boot_both_result, MC_n ~ tumour+tissue, mean)
boot_both_result_agg_MC$Age <- "MC"
colnames(boot_both_result_agg_MC) <- c("Cohort", "Source", "Count", "Age")
boot_both_result_agg_Mixed <- aggregate(data=boot_both_result, Mixed_n ~ tumour+tissue, mean)
boot_both_result_agg_Mixed$Age <- "Mixed"
colnames(boot_both_result_agg_Mixed) <- c("Cohort", "Source", "Count", "Age")

boot_both_result_agg <- rbind(boot_both_result_agg_UC, boot_both_result_agg_MC, boot_both_result_agg_Mixed)
boot_both_result_agg$Source[boot_both_result_agg$Source == "normal"] <- "Normal_random"
boot_both_result_agg$Source[boot_both_result_agg$Source == "tumour"] <- "Tumour_random"
colnames(boot_both_result_agg)[4] <- "Age Enrichment Category"

mod.age.agg.df <- rbind(mod.age.agg.df, boot_both_result_agg)

mod.age.agg.df$`Age Enrichment Category` <- factor(mod.age.agg.df$`Age Enrichment Category`,
                                                   levels=c("UC", "Mixed", "MC"))

mod.age.agg.df$Source <- factor(mod.age.agg.df$Source, levels=c("Normal", "Normal_random", "Tumour", "Tumour_random"))

### Supp Fig S5: Fraction of age enrichment categories among tumour and normal modules
pdf("Figure_S5_with_random.pdf", height=5, width=6)
g <- ggplot(mod.age.agg.df, aes(fill=`Age Enrichment Category`, x=Source, y=Count) )+
  geom_bar(position="fill", stat = "summary", fun="mean") + 
  labs(x="Sample Type", y= "Fraction of Co-Expression Modules") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()
