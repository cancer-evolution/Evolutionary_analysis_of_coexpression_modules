##Script for all panels of Figure 1, Figures S1-S5 and Table S1.

library(ggplot2)
library(reshape2)
library(knitr)
library(tidyverse)
library(data.table)

load("Objects/expression_WGCNA_paired.Rdata")

load("Objects/cluster_assignments_high.Rdata")

n_size_sample <- vector()
##Number and size of modules
for(tissue in names(cluster_assignments_high)){
  for(tumour in names(cluster_assignments_high[[tissue]])){
    temp <- cluster_assignments_high[[tissue]][[tumour]]
    temp$grey <- NULL
    n_modules <- length(names(temp)) 
    median_size <- median(sapply(temp, length))
    n_samples <- nrow(expression_WGCNA[[tumour]][[tissue]])
    n_size_sample <- rbind(n_size_sample, c(tissue, tumour, n_modules, median_size, n_samples))
  }
}
colnames(n_size_sample) <- c("Tissue", "Tumour", "N_modules", "Median_module_size", "N_samples")
n_size_sample <- as.data.frame(n_size_sample)
n_size_sample$N_modules <- as.numeric(as.character(n_size_sample$N_modules))
n_size_sample$Median_module_size <- as.numeric(as.character(n_size_sample$Median_module_size))
n_size_sample$N_samples <- as.numeric(as.character(n_size_sample$N_samples))

library(ggplot2)


tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

## Figure 1D: Number of co-expression modules by age enrichment category in each subtype
load("Objects/age_enrichment_high.Rdata")

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

mod.age.agg.df$`Age Enrichment Category` <- str_replace(mod.age.agg.df$`Age Enrichment Category`, "UC", "UC Enriched")
mod.age.agg.df$`Age Enrichment Category` <- str_replace(mod.age.agg.df$`Age Enrichment Category`, "MC", "MC Enriched")
mod.age.agg.df$`Age Enrichment Category` <- str_replace(mod.age.agg.df$`Age Enrichment Category`, "Mixed", "Mixed UC-MC")

mod.age.agg.df$`Age Enrichment Category` <- factor(mod.age.agg.df$`Age Enrichment Category`,
                                                   levels=c("UC Enriched", "Mixed UC-MC", "MC Enriched"))

pdf("Figure_1D_high.pdf", height=3.85, width=5)
g <- ggplot(mod.age.agg.df, aes(fill=`Age Enrichment Category`, x=Source, y=Count) ) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(jitter.width=0.20), size=1) + 
  labs(x="Sample Type", y= "Number of Co-Expression Modules") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

