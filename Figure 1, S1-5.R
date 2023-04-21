##Script for all panels of Figure 1, Figures S1-S5 and Table S1.

library(ggplot2)
library(reshape2)
library(knitr)
library(tidyverse)
library(data.table)

load("expression_WGCNA_paired.Rdata")

load("cluster_assignments.Rdata")

n_size_sample <- vector()
##Number and size of modules
for(tissue in names(cluster_assignments)){
  for(tumour in names(cluster_assignments[[tissue]])){
    temp <- cluster_assignments[[tissue]][[tumour]]
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

##Figure S2
pdf("Figure_S2_top.pdf", height=3.85, width=8)
g <- ggplot(n_size_sample, aes(x=N_samples, y=N_modules))+
  geom_point(aes(colour=Tumour), size=3)+
  facet_grid(.~Tissue)+
  theme_bw()+
  ylab("Number of modules")+
  xlab("Number of samples")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()


##Figure S2
pdf("Figure_S2_bottom.pdf", height=3.85, width=8)
g <- ggplot(n_size_sample, aes(x=N_samples, y=Median_module_size))+
  geom_point(aes(colour=Tumour), size=3)+
  facet_grid(.~Tissue)+
  theme_bw()+
  ylab("Median module size")+
  xlab("Number of samples")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

cor.test(n_size_sample$N_samples[n_size_sample$Tissue == "tumour"], 
         n_size_sample$N_modules[n_size_sample$Tissue == "tumour"], method="sp")

cor.test(n_size_sample$N_samples[n_size_sample$Tissue == "normal"], 
         n_size_sample$N_modules[n_size_sample$Tissue == "normal"], method="sp")


cor.test(n_size_sample$N_samples[n_size_sample$Tissue == "tumour"], 
         n_size_sample$Median_module_size[n_size_sample$Tissue == "tumour"], method="sp")

cor.test(n_size_sample$N_samples[n_size_sample$Tissue == "normal"], 
         n_size_sample$Median_module_size[n_size_sample$Tissue == "normal"], method="sp")


tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

##Figure 1B
pdf("Figure_1B.pdf", height=3.85, width=5)
g <- ggplot(n_size_sample, aes(x=Tissue, y=N_modules)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(aes(colour = Tumour),  width=0.12, size=1.5)  + 
  labs(x="Sample Type", y= "Number of Co-Expression Modules") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

##Figure S1
pdf("Figure_S1.pdf", height=3.85, width=4.25)
g <- ggplot(subset(n_size_sample, Tumour %in% tumours_normal), aes(x=Tissue, y=N_modules))+
  geom_boxplot()+
  geom_jitter(aes(colour=Tumour), width=0.12)+
  ylab("Number of co-expression modules")+
  xlab("Sample type")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

sum(n_size_sample$N_modules[n_size_sample$Tissue == "tumour"])
sum(n_size_sample$N_modules[n_size_sample$Tissue == "normal"])

mean(n_size_sample$N_modules[n_size_sample$Tissue == "tumour"])
mean(n_size_sample$N_modules[n_size_sample$Tissue == "normal"])

wilcox.test(n_size_sample$N_modules[n_size_sample$Tissue == "tumour"], 
            n_size_sample$N_modules[n_size_sample$Tissue == "normal"])


wilcox.test(n_size_sample$N_modules[n_size_sample$Tissue == "tumour" & n_size_sample$Tumour %in% tumours_normal], 
            n_size_sample$N_modules[n_size_sample$Tissue == "normal" & n_size_sample$Tumour %in% tumours_normal])

##Figre 1C
pdf("Figure_1C.pdf", height=3.85, width=5)
g <- ggplot(n_size_sample, aes(x=Tissue, y=Median_module_size))+
  geom_boxplot(outlier.shape=NA)+ 
  geom_jitter(aes(colour = Tumour),  width=0.12, size=1.5)+ 
  labs(x="Sample Type", y= "Median Module Size") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()


##Supplementary Figure S3
pdf("Figure_S3.pdf", height=3.85, width=4.25)
g <- ggplot(subset(n_size_sample, Tumour %in% tumours_normal), aes(x=Tissue, y=Median_module_size))+
  geom_boxplot()+
  geom_jitter(aes(colour=Tumour), width=0.12)+
  ylab("Median module size")+
  xlab("Sample type")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

mean(n_size_sample$Median_module_size[n_size_sample$Tissue == "tumour"])
mean(n_size_sample$Median_module_size[n_size_sample$Tissue == "normal"])

wilcox.test(n_size_sample$Median_module_size[n_size_sample$Tissue == "tumour"], 
            n_size_sample$Median_module_size[n_size_sample$Tissue == "normal"])


wilcox.test(n_size_sample$N_modules[n_size_sample$Tissue == "tumour" & n_size_sample$Tumour %in% tumours_normal], 
            n_size_sample$N_modules[n_size_sample$Tissue == "normal" & n_size_sample$Tumour %in% tumours_normal])

sum(n_size_sample$N_modules[n_size_sample$Tissue == "tumour"])
sum(n_size_sample$N_modules[n_size_sample$Tissue == "normal"])
sum(n_size_sample$N_modules)

### Supp Fig S4: Fraction of genes unassigned to modules

##### Count total number of genes in each module, across all cohorts
##### Calculate proportion that are in the grey (unassigned) module for each cohort
total.genes.in.modules <- function(X, incl.grey=F ) {  
  tl <- sum( unlist( lapply( X, length) ) )
  
  if( !incl.grey ) {
    grey_length <- length( X$grey )
    tl = tl - grey_length
  }
  
  prct.Assigned = ( tl / (tl+grey_length) ) * 100
  
  return( c( "Assigned"=tl, "Unassigned"=grey_length, "Percent.Unassigned"= 100 - prct.Assigned, "Percent.Assigned"=prct.Assigned ) )
}

num.genes.in.modules.tum <- lapply( cluster_assignments$tumour, total.genes.in.modules )
num.genes.in.modules.tum <- t( data.frame( num.genes.in.modules.tum ) )
num.genes.in.modules.norm <- lapply( cluster_assignments$normal, total.genes.in.modules)
num.genes.in.modules.norm <- t( data.frame( num.genes.in.modules.norm ) )

m <- melt(num.genes.in.modules.tum[,-c(1,2)]) 
names(m) = c("Cohort","Status","Num.Genes")

pdf("Figure_S4.pdf", height=3, width=7)
g <- ggplot(m, aes(fill=Status, y=Num.Genes, x=Cohort) ) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("Percentage of transcriptome")
print(g)
dev.off()

## Figure 1D: Number of co-expression modules by age enrichment category in each subtype
load("age_enrichment.Rdata")

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

pdf("Figure_1D.pdf", height=3.85, width=5)
g <- ggplot(mod.age.agg.df, aes(fill=`Age Enrichment Category`, x=Source, y=Count) ) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(jitter.width=0.20), size=1) + 
  labs(x="Sample Type", y= "Number of Co-Expression Modules") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()



summary(mod.age.agg.df)
  
mean(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "Mixed UC-MC" &
                 mod.age.agg.df$Source == "Tumour","Count"])

sum(mod.age.agg.df[mod.age.agg.df$Source == "Tumour","Count"])
sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "Mixed UC-MC" &
                      mod.age.agg.df$Source == "Tumour","Count"])
774/1166

sum(mod.age.agg.df[mod.age.agg.df$Source == "Normal","Count"])
sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "Mixed UC-MC" &
                     mod.age.agg.df$Source == "Normal","Count"])
137/337

sum(mod.age.agg.df[mod.age.agg.df$Source == "Tumour","Count"])
sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "UC Enriched" &
                     mod.age.agg.df$Source == "Tumour","Count"])
185/1166

sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "MC Enriched" &
                     mod.age.agg.df$Source == "Tumour","Count"])
207/1166



sum(mod.age.agg.df[mod.age.agg.df$Source == "Normal","Count"])
sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "UC Enriched" &
                     mod.age.agg.df$Source == "Normal","Count"])
84/337

sum(mod.age.agg.df[mod.age.agg.df$Source == "Normal","Count"])
sum(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "MC Enriched" &
                     mod.age.agg.df$Source == "Normal","Count"])
116/337


mean(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "UC Enriched" &
                      mod.age.agg.df$Source == "Tumour","Count"])
mean(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "MC Enriched" &
                      mod.age.agg.df$Source == "Tumour","Count"])

mean(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "UC Enriched" &
                      mod.age.agg.df$Source == "Normal","Count"])
mean(mod.age.agg.df[mod.age.agg.df$`Age Enrichment Category` == "MC Enriched" &
                      mod.age.agg.df$Source == "Normal","Count"])


##Supp Table 1
module_ages_print <- vector()
for(tumour in names(age_enrichment)){
  module_ages_print <- rbind(module_ages_print, age_enrichment[[tumour]])
}
module_ages_print <- module_ages_print[,c("tumour", "tissue_type", "cluster", "Module_age")]
colnames(module_ages_print) <- c("Tumour_type", "Tissue_type", "Cluster", "Module_age")
module_ages_print <- module_ages_print[module_ages_print$Cluster != "grey",]

temp_UC <- module_ages_print[module_ages_print$Module_age == "UC",]
module_ages_print_UC <- as.data.frame(table(temp_UC[,c("Tumour_type", "Tissue_type")]))
module_ages_print_UC <- module_ages_print_UC[module_ages_print_UC$Freq != 0,]
summary(module_ages_print_UC$Freq[module_ages_print_UC$Tissue_type == "tumour"])
sum(module_ages_print_UC$Freq[module_ages_print_UC$Tissue_type == "tumour"])


summary(module_ages_print_UC$Freq[module_ages_print_UC$Tissue_type == "normal"])
sum(module_ages_print_UC$Freq[module_ages_print_UC$Tissue_type == "normal"])


temp_MC <- module_ages_print[module_ages_print$Module_age == "MC",]
module_ages_print_MC <- as.data.frame(table(temp_MC[,c("Tumour_type", "Tissue_type")]))
module_ages_print_MC <- module_ages_print_MC[module_ages_print_MC$Freq != 0,]
summary(module_ages_print_MC$Freq[module_ages_print_MC$Tissue_type == "tumour"])
sum(module_ages_print_MC$Freq[module_ages_print_MC$Tissue_type == "tumour"])
summary(module_ages_print_MC$Freq[module_ages_print_MC$Tissue_type == "normal"])
sum(module_ages_print_MC$Freq[module_ages_print_MC$Tissue_type == "normal"])


temp_Mixed <- module_ages_print[module_ages_print$Module_age == "Mixed",]
module_ages_print_Mixed <- as.data.frame(table(temp_Mixed[,c("Tumour_type", "Tissue_type")]))
module_ages_print_Mixed <- module_ages_print_Mixed[module_ages_print_Mixed$Freq != 0,]
summary(module_ages_print_Mixed$Freq[module_ages_print_Mixed$Tissue_type == "tumour"])
sum(module_ages_print_Mixed$Freq[module_ages_print_Mixed$Tissue_type == "tumour"])
summary(module_ages_print_Mixed$Freq[module_ages_print_Mixed$Tissue_type == "normal"])
sum(module_ages_print_Mixed$Freq[module_ages_print_Mixed$Tissue_type == "normal"])







## Figure 1E: Most frequently enriched GO terms by modules by age enrichment category
load('module_function_all.Rdata')
load("age_enrichment.Rdata")
load("all_preservation_t_to_n2.Rdata")
load("all_preservation_n_to_t2.Rdata")

##### exclude cancer and other disease related terms
to_exclude <- c("Huntington disease", "Human cytomegalovirus infection","Alzheimer disease","Staphylococcus aureus infection",                                                                                                    
                "Central carbon metabolism in cancer","Hepatitis C","MicroRNAs in cancer","Proteoglycans in cancer","Measles",                                                                                                            
                "Pertussis","Small cell lung cancer","Hepatocellular carcinoma","Viral myocarditis","Parkinson disease",                                                                                                                 
                "Insulin resistance","Human papillomavirus infection","Diseases of signal transduction","Chronic myeloid leukemia",                                                                                                   
                "HIV Infection","Pathways in cancer","Loss of Function of SMAD4 in Cancer","Vibrio cholerae infection",                                                                                               
                "Platinum drug resistance","Hypertrophic cardiomyopathy (HCM)","Gastric cancer","Non-small cell lung cancer",                                                                                                          
                "Toxoplasmosis","Epstein-Barr virus infection","Renal cell carcinoma","Shigellosis","Basal cell carcinoma",                                                                                                                
                "Colorectal cancer","Endometrial cancer","Human immunodeficiency virus 1 infection",                                                                                            
                "Kaposi sarcoma-associated herpesvirus infection","Herpes simplex virus 1 infection",
                "AGE-RAGE signaling pathway in diabetic complications","Defective CFTR causes cystic fibrosis",                                                                              
                "Dilated cardiomyopathy (DCM)","Hepatitis B","Central carbon metabolism in cancer","MicroRNAs in cancer",   
                "Proteoglycans in cancer","Small cell lung cancer","Pathways in cancer","Gastric cancer","Non-small cell lung cancer",              
                "Colorectal cancer","Endometrial cancer","Transcriptional misregulation in cancer","Huntington disease",
                "Parkinson disease","Alzheimer disease","Non-alcoholic fatty liver disease (NAFLD)","Infectious disease",
                "Chagas disease (American trypanosomiasis)","Signaling by FGFR1 in disease",
                "SMAD4 MH2 Domain Mutants in Cancer")

module_function_all2 <- vector()
for(tissue in c("normal", "tumour")){
  for(tumour in names(cluster_assignments[[tissue]])){
    module_function <- vector()
    for(mod in names(cluster_assignments[[tissue]][[tumour]])){
      temp <- module_function_all[module_function_all$Tissue == tissue &
                                    module_function_all$Tumour == tumour &
                                    module_function_all$Module == mod,]
      age <- temp$Module_age[1]
      novelty <- as.character(temp$Novelty[1])
      if(!grepl(paste(to_exclude, collapse="|"), temp$Pathway[1])){
        pathway <- temp$Pathway[1]
      }else if(!grepl(paste(to_exclude, collapse="|"), temp$Pathway[2])){
        pathway <- temp$Pathway[2]
      }else{
        pathway <- "empty"
      }
      module_function_all2 <- rbind(module_function_all2,
                                    c(tissue, tumour, mod, pathway, age, novelty))
    }
  }
}
colnames(module_function_all2) <- c("Tissue", "Tumour", "Module", "Pathway", "Age", "Novelty")

module_function_all2 <- as.data.frame(module_function_all2)
module_function_all2 <- module_function_all2[module_function_all2$Pathway != "empty",]
module_function_all2$Pathway[module_function_all2$Pathway == "Cell Cycle"] <- "Cell cycle"
module_function_all2 <- module_function_all2[!is.na(module_function_all2$Tumour),]


##### Identify most frequently enriched biological terms across age enrichment categories 
sum.topKR.mixed <- summary(as.factor( module_function_all2$Pathway[module_function_all2$Age == "Mixed"] ) ) 

kable( head( sum.topKR.mixed,10) )

topKR.mixed <- names( head( sum.topKR.mixed,10) )

sum.topKR.MC <- summary(as.factor( module_function_all2$Pathway[module_function_all2$Age == "MC"] ) ) 

kable( head( sum.topKR.MC,10) )

topKR.MC <- names( head( sum.topKR.MC,10) )


sum.topKR.UC <- summary(as.factor( module_function_all2$Pathway[module_function_all2$Age == "UC"] ) ) 

kable( head( sum.topKR.UC,10) )

topKR.UC <- names( head( sum.topKR.UC,10) )


##### Aggregate term enrichment counts for plotting

n_times_pathway_age <- vector()
#pathway_order <- unique(module_function_all2$Pathway)

for(tissue in c("normal", "tumour")){
  for(age in c("UC", "Mixed", "MC")){
    temp <- module_function_all2[module_function_all2$Tissue == tissue & module_function_all2$Age == age,]
    temp <- unique(temp[,c("Tissue", "Tumour", "Pathway")])
    temp2 <- table(temp$Pathway)
    n_times_pathway_age <- rbind(n_times_pathway_age,
                                 cbind(names(temp2), temp2, tissue, age ))
  }
}

colnames(n_times_pathway_age) <- c("Pathways", "Count", "Tissue", "Age")
n_times_pathway_age <- as.data.frame(n_times_pathway_age)
n_times_pathway_age$Count <- as.numeric(as.character(n_times_pathway_age$Count))

total_sum <- aggregate(Count ~ Pathways, n_times_pathway_age , sum)
total_sum <- total_sum[total_sum$Count != 1,]
total_sum <- total_sum[total_sum$Count != 0,]
total_sum <- total_sum$Pathways[total_sum$Count != 2]
n_times_pathway_age <- n_times_pathway_age[n_times_pathway_age$Pathways %in% total_sum,]


##### Plot as facets showing # of times each top GO term appeared in enrichment results
top.KR.Pathways <- c( topKR.UC, topKR.mixed, topKR.MC )

n_times_pathway_age$Age <- factor(n_times_pathway_age$Age,
                                  levels=c("UC", "Mixed", "MC"))

pdf("Figure_1E.pdf", height=4, width=9)
g <- ggplot(filter(n_times_pathway_age, Tissue=="tumour", Pathways %in% top.KR.Pathways), aes(x=Pathways, y=Count))+
  geom_bar(stat='identity', aes(fill=Age))+ ## ,col=c("orange","blue")
  facet_grid(.~Age)+
  ylab("Number of tumour types")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()


### Supp Fig S5: Fraction of age enrichment categories among tumour and normal modules
pdf("Figure_S5.pdf", height=5, width=5)
g <- ggplot(mod.age.agg.df, aes(fill=`Age Enrichment Category`, x=Source, y=Count) )+
  geom_bar(position="fill", stat = "summary", fun="mean") + 
  labs(x="Sample Type", y= "Fraction of Co-Expression Modules") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

fisher.test( cbind( c( 774, 1166-774), c(137, 337-137)))$p.value

fisher.test( cbind( c( 185, 774, 207 ), c(84, 137, 116)))
