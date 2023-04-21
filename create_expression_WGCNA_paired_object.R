###WGCNA analysis

library(WGCNA)
library(edgeR)

source("functions.R")

#Will use only paired samples -> maybe difference in sample size leads to differences in
#the ability to detect modules


##read in data
tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

#Keep only paired samples and remove lowly expressed genes (cpm<1 in all tumour or normal samples)
expression <- list()
for(tumour in tumours){
  expression[[tumour]] <- read_expression_paired(tumour)
  print(tumour)
}

#Only 1 SKMC sample, so will remove

expression$SKCM$normal <- NA

#Checking of quality of samples
#Only removes genes with no variance or with too many missing samples.
for(tumour in tumours){
  for(tissue_type in c("normal", "tumour")){
    if(!is.na(expression[[tumour]][[tissue_type]])){
      gsg <- goodSamplesGenes(t(expression[[tumour]][[tissue_type]]), verbose=3)
      if(gsg$allOK != TRUE){
        expression[[tumour]][[tissue_type]] = expression[[tumour]][[tissue_type]][gsg$goodSamples, gsg$goodGenes]
      }
      gsg <- goodSamplesGenes(expression[[tumour]][[tissue_type]], verbose=3)
      print(gsg$allOK)
    }
  }
}
##Ignore warning

#Detecting outliers in the samples
pdf("clustering_of_samples_paired.pdf", width=10)
for(tissue_type in c("normal", "tumour")){
  for(tumour in tumours){
    if(!is.na(expression[[tumour]][[tissue_type]])){
      sampleTree = hclust(dist(t(expression[[tumour]][[tissue_type]])), method = "average");
      plot(sampleTree, main = paste("Sample clustering to detect outliers", tissue_type, tumour, sep="\n"), sub="", xlab="")
    }
  }
}
dev.off()

#Remove outliers of each sample type

expression_WGCNA <- list()
pdf("clustering_of_samples_with_cut_paired.pdf", width=10)
expression_WGCNA[["BLCA"]][["normal"]] <- t(expression[["BLCA"]][["normal"]])
expression_WGCNA[["BRCA"]][["normal"]] <- remove_outlier_samples(t(expression[["BRCA"]][["normal"]]), paste("BRCA", "normal"), 500000, 1)
expression_WGCNA[["COAD"]][["normal"]] <- remove_outlier_samples(t(expression[["COAD"]][["normal"]]), paste("COAD", "normal"), 1000000, 1)
expression_WGCNA[["ESCA"]][["normal"]] <- remove_outlier_samples(t(expression[["ESCA"]][["normal"]]), paste("ESCA", "normal"), 1500000, 1)
expression_WGCNA[["HNSC"]][["normal"]] <- remove_outlier_samples(t(expression[["HNSC"]][["normal"]]), paste("HNSC", "normal"), 1500000, 1:2)
expression_WGCNA[["KICH"]][["normal"]] <- remove_outlier_samples(t(expression[["KICH"]][["normal"]]), paste("KICH", "normal"), 600000, 1)
expression_WGCNA[["KIRC"]][["normal"]] <- remove_outlier_samples(t(expression[["KIRC"]][["normal"]]), paste("KIRC", "normal"), 400000, 1)
expression_WGCNA[["KIRP"]][["normal"]] <- t(expression[["KIRP"]][["normal"]])
expression_WGCNA[["LIHC"]][["normal"]] <- remove_outlier_samples(t(expression[["LIHC"]][["normal"]]), paste("LIHC", "normal"), 4000000, 1)
expression_WGCNA[["LUAD"]][["normal"]] <- remove_outlier_samples(t(expression[["LUAD"]][["normal"]]), paste("LUAD", "normal"), 750000, 1)
expression_WGCNA[["LUSC"]][["normal"]] <- remove_outlier_samples(t(expression[["LUSC"]][["normal"]]), paste("LUSC", "normal"), 800000, 1)
expression_WGCNA[["PRAD"]][["normal"]] <- remove_outlier_samples(t(expression[["PRAD"]][["normal"]]), paste("PRAD", "normal"), 800000, 1)
expression_WGCNA[["READ"]][["normal"]] <- t(expression[["READ"]][["normal"]])
expression_WGCNA[["STAD"]][["normal"]] <- remove_outlier_samples(t(expression[["STAD"]][["normal"]]), paste("STAD", "normal"), 2000000, 1)
expression_WGCNA[["THCA"]][["normal"]] <- remove_outlier_samples(t(expression[["THCA"]][["normal"]]), paste("THCA", "normal"), 1500000, 1)
expression_WGCNA[["UCEC"]][["normal"]] <- t(expression[["UCEC"]][["normal"]])

expression_WGCNA[["ACC"]][["tumour"]] <- remove_outlier_samples(t(expression[["ACC"]][["tumour"]]), paste("ACC", "tumour"), 1500000, 1)
expression_WGCNA[["BLCA"]][["tumour"]] <- t(expression[["BLCA"]][["tumour"]])
expression_WGCNA[["BRCA"]][["tumour"]] <- remove_outlier_samples(t(expression[["BRCA"]][["tumour"]]), paste("BRCA", "tumour"), 850000, 1)
expression_WGCNA[["CHOL"]][["tumour"]] <- t(expression[["CHOL"]][["tumour"]])
expression_WGCNA[["CESC"]][["tumour"]] <- remove_outlier_samples(t(expression[["CESC"]][["tumour"]]), paste("CESC", "tumour"), 700000, 1)
expression_WGCNA[["COAD"]][["tumour"]] <- remove_outlier_samples(t(expression[["COAD"]][["tumour"]]), paste("COAD", "tumour"), 450000, 1)
expression_WGCNA[["ESCA"]][["tumour"]] <- t(expression[["ESCA"]][["tumour"]])
expression_WGCNA[["GBM"]][["tumour"]] <- remove_outlier_samples(t(expression[["GBM"]][["tumour"]]), paste("GBM", "tumour"), 700000, 1)
expression_WGCNA[["HNSC"]][["tumour"]] <- remove_outlier_samples(t(expression[["HNSC"]][["tumour"]]), paste("HNSC", "tumour"), 1250000, 1)
expression_WGCNA[["KICH"]][["tumour"]] <- remove_outlier_samples(t(expression[["KICH"]][["tumour"]]), paste("KICH", "tumour"), 600000, 1)
expression_WGCNA[["KIRC"]][["tumour"]] <- remove_outlier_samples(t(expression[["KIRC"]][["tumour"]]), paste("KIRC", "tumour"), 800000, 1)
expression_WGCNA[["KIRP"]][["tumour"]] <- remove_outlier_samples(t(expression[["KIRP"]][["tumour"]]), paste("KIRP", "tumour"), 600000, 1)
expression_WGCNA[["LGG"]][["tumour"]] <- remove_outlier_samples(t(expression[["LGG"]][["tumour"]]), paste("LGG", "tumour"), 1500000, 1)
expression_WGCNA[["LIHC"]][["tumour"]] <- remove_outlier_samples(t(expression[["LIHC"]][["tumour"]]), paste("LIHC", "tumour"), 2500000, 1)
expression_WGCNA[["LUAD"]][["tumour"]] <- t(expression[["LUAD"]][["tumour"]])
expression_WGCNA[["LUSC"]][["tumour"]] <- remove_outlier_samples(t(expression[["LUSC"]][["tumour"]]), paste("LUSC", "tumour"), 1000000, 1)
expression_WGCNA[["MESO"]][["tumour"]] <- remove_outlier_samples(t(expression[["MESO"]][["tumour"]]), paste("MESO", "tumour"), 1500000, 1)
expression_WGCNA[["OV"]][["tumour"]] <- remove_outlier_samples(t(expression[["OV"]][["tumour"]]), paste("OV", "tumour"), 800000, 1)
expression_WGCNA[["PAAD"]][["tumour"]] <- remove_outlier_samples(t(expression[["PAAD"]][["tumour"]]), paste("PAAD", "tumour"), 900000, 1)
expression_WGCNA[["PCPG"]][["tumour"]] <- remove_outlier_samples(t(expression[["PCPG"]][["tumour"]]), paste("PCPG", "tumour"), 1600000, 1)
expression_WGCNA[["PRAD"]][["tumour"]] <- t(expression[["PRAD"]][["tumour"]])
expression_WGCNA[["READ"]][["tumour"]] <- t(expression[["READ"]][["tumour"]])
expression_WGCNA[["SARC"]][["tumour"]] <- remove_outlier_samples(t(expression[["SARC"]][["tumour"]]), paste("SARC", "tumour"), 1600000, 1)
expression_WGCNA[["SKCM"]][["tumour"]] <- t(expression[["SKCM"]][["tumour"]])
expression_WGCNA[["STAD"]][["tumour"]] <- remove_outlier_samples(t(expression[["STAD"]][["tumour"]]), paste("STAD", "tumour"), 800000, 1)
expression_WGCNA[["TGCT"]][["tumour"]] <- remove_outlier_samples(t(expression[["TGCT"]][["tumour"]]), paste("TGCT", "tumour"), 1000000, 1)
expression_WGCNA[["THCA"]][["tumour"]] <- t(expression[["THCA"]][["tumour"]])
expression_WGCNA[["THYM"]][["tumour"]] <- remove_outlier_samples(t(expression[["THYM"]][["tumour"]]), paste("THYM", "tumour"), 700000, 1)
expression_WGCNA[["UCEC"]][["tumour"]] <- t(expression[["UCEC"]][["tumour"]])
expression_WGCNA[["UCS"]][["tumour"]] <- remove_outlier_samples(t(expression[["UCS"]][["tumour"]]), paste("UCS", "tumour"), 900000, 1)
expression_WGCNA[["UVM"]][["tumour"]] <- remove_outlier_samples(t(expression[["UVM"]][["tumour"]]), paste("UVM", "tumour"), 800000, 1)
dev.off()


#save(expression_WGCNA, file="expression_WGCNA_paired.Rdata")