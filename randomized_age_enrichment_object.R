source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

#Analysing all tumours, not only a subset
tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

load("Objects/cluster_assignments.Rdata")

gene_ages <- read.csv("Text_files/geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_in_clusters <- unique(unname(unlist(cluster_assignments)))
gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 1:3),1]
MC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 4:16),1]

##Randomize genes to clusters

boot_tumour_result <- vector()
boot_normal_result <- vector()
cluster_assignments_random <- list()
age_enrichment_random <- vector()

for(tumour in tumours){
  local_tumour <- cluster_assignments$tumour[[tumour]]
  if("grey" %in% names(local_tumour)){
    local_tumour$grey <- NULL
  }
  
  genes_tumour <- unname(unlist(local_tumour))
  module_size_tumour <- sapply(local_tumour, length)
  
  UC_tumour_total <- genes_tumour[genes_tumour %in% UC_genes]
  MC_tumour_total <- genes_tumour[genes_tumour %in% MC_genes]
  
  
  if(tumour %in% tumours_normal){
    local_normal <- cluster_assignments$normal[[tumour]]
    if("grey" %in% names(local_normal)){
      local_normal$grey <- NULL
    }
    genes_normal <- unname(unlist(local_normal))
    module_size_normal <- sapply(local_normal, length)
    
    UC_normal_total <- genes_normal[genes_normal %in% UC_genes]
    MC_normal_total <- genes_normal[genes_normal %in% MC_genes]
    
  }
  
  cluster_assignments_random[["tumour"]][[tumour]] <- list()
  cluster_assignments_random[["normal"]][[tumour]] <- list()
  for(i in 1:1000){
    age_enrichment_random_tumour <- vector()
    for(mod in names(module_size_tumour)){
      cluster_assignments_random[["tumour"]][[tumour]][[paste(mod, i, sep="_")]] <- sample(genes_tumour, module_size_tumour[mod])
      
      genes_in_cluster <- cluster_assignments_random[["tumour"]][[tumour]][[paste(mod, i, sep="_")]]
      UC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% UC_genes]
      
      MC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% MC_genes]
      
      df <- cbind(c(length(UC_genes_cluster), length(MC_genes_cluster)),
                  c(length(UC_tumour_total), length(MC_tumour_total)))
      UC_p <- fisher.test(df, alternative="greater")$p.value
      MC_p <- fisher.test(df, alternative="less")$p.value
      
      age_enrichment_random_tumour <- rbind(age_enrichment_random_tumour, c(tumour, "tumour", mod, i, UC_p, MC_p))
    }
    
    colnames(age_enrichment_random_tumour) <- c("tumour", "tissue_type", "cluster", "iteration", "p_UC_enrichment", "p_MC_enrichment")
    age_enrichment_random_tumour <- as.data.frame(age_enrichment_random_tumour)
    age_enrichment_random_tumour$p_UC_enrichment <- as.numeric(as.character(age_enrichment_random_tumour$p_UC_enrichment))
    age_enrichment_random_tumour$p_MC_enrichment <- as.numeric(as.character(age_enrichment_random_tumour$p_MC_enrichment))
    age_enrichment_random_tumour$p_UC_adj <- p.adjust(age_enrichment_random_tumour$p_UC_enrichment, method="BH")
    age_enrichment_random_tumour$p_MC_adj <- p.adjust(age_enrichment_random_tumour$p_MC_enrichment, method="BH")
    age_enrichment_random_tumour$Module_age <- ifelse(age_enrichment_random_tumour$p_UC_adj < 0.05, "UC", 
                                                  ifelse(age_enrichment_random_tumour$p_MC_adj < 0.05, "MC",
                                                   "Mixed"))
    age_enrichment_random <- rbind(age_enrichment_random, age_enrichment_random_tumour)
    
    boot_tumour_result <- rbind(boot_tumour_result,
                                c(tumour, "tumour", i, sum(age_enrichment_random_tumour$Module_age == "UC"),
                                  sum(age_enrichment_random_tumour$Module_age == "MC"),
                                  sum(age_enrichment_random_tumour$Module_age == "Mixed"),
                                  nrow(age_enrichment_random_tumour)))
    
    
    if(tumour %in% tumours_normal){
      age_enrichment_random_normal <- vector()
      for(mod in names(module_size_normal)){
        cluster_assignments_random[["normal"]][[tumour]][[paste(mod, i, sep="_")]] <- sample(genes_normal, module_size_normal[mod])
        
        genes_in_cluster <- cluster_assignments_random[["normal"]][[tumour]][[paste(mod, i, sep="_")]]
        UC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% UC_genes]
        
        MC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% MC_genes]
        
        df <- cbind(c(length(UC_genes_cluster), length(MC_genes_cluster)),
                    c(length(UC_normal_total), length(MC_normal_total)))
        UC_p <- fisher.test(df, alternative="greater")$p.value
        MC_p <- fisher.test(df, alternative="less")$p.value
        
        age_enrichment_random_normal <- rbind(age_enrichment_random_normal, c(tumour, "normal", mod, i, UC_p, MC_p))
      }
      
      colnames(age_enrichment_random_normal) <- c("tumour", "tissue_type", "cluster", "iteration", "p_UC_enrichment", "p_MC_enrichment")
      age_enrichment_random_normal <- as.data.frame(age_enrichment_random_normal)
      age_enrichment_random_normal$p_UC_enrichment <- as.numeric(as.character(age_enrichment_random_normal$p_UC_enrichment))
      age_enrichment_random_normal$p_MC_enrichment <- as.numeric(as.character(age_enrichment_random_normal$p_MC_enrichment))
      age_enrichment_random_normal$p_UC_adj <- p.adjust(age_enrichment_random_normal$p_UC_enrichment, method="BH")
      age_enrichment_random_normal$p_MC_adj <- p.adjust(age_enrichment_random_normal$p_MC_enrichment, method="BH")
      age_enrichment_random_normal$Module_age <- ifelse(age_enrichment_random_normal$p_UC_adj < 0.05, "UC", 
                                              ifelse(age_enrichment_random_normal$p_MC_adj < 0.05, "MC",
                                                     "Mixed"))
      age_enrichment_random <- rbind(age_enrichment_random, age_enrichment_random_normal)
      
      boot_normal_result <- rbind(boot_normal_result,
                                  c(tumour, "normal", i, sum(age_enrichment_random_normal$Module_age == "UC"),
                                    sum(age_enrichment_random_normal$Module_age == "MC"),
                                    sum(age_enrichment_random_normal$Module_age == "Mixed"),
                                    nrow(age_enrichment_random_normal))) 
    }
  }
  print(tumour)
}
save(cluster_assignments_random, file="Objects/cluster_assignments_random.Rdata")
save(age_enrichment_random, file="Objects/age_enrichment_random.Rdata")

colnames(boot_tumour_result) <- c("tumour", "tissue", "i", "UC_n", "MC_n", "Mixed_n", "Total") 
colnames(boot_normal_result) <- c("tumour", "tissue", "i", "UC_n", "MC_n", "Mixed_n", "Total")

save(boot_tumour_result,file= "Objects/boot_tumour_result.Rdata")
save(boot_normal_result, file="Objects/boot_normal_result.Rdata")

