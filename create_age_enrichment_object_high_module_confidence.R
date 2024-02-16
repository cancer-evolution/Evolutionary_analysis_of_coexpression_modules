##Create age_enrichment.Rdata object

source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

#Analysing all tumours, not only a subset
tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")

load("Objects/cluster_assignments_high.Rdata")
gene_ages <- read.csv("Text_files/geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")
genes_in_clusters <- unique(unname(unlist(cluster_assignments_high)))

gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 1:3),1]
MC_genes <- gene_ages[which(gene_ages$Phylostrata %in% 4:16),1]


#Enrichment tests (taking into account that each sample type has a diff number of total UC/MC genes)
age_enrichment <- list()

for(tumour in tumours){
  temp <- list()
  for(tissue_type in c("normal", "tumour")){
    if(tissue_type == "tumour" || (tumour %in% tumours_normal & tissue_type == "normal")){
      cluster_names <- names(cluster_assignments_high[[tissue_type]][[tumour]])
      all_genes <- unique(unname(unlist(cluster_assignments_high[[tissue_type]][[tumour]])))
      UC_genes_total <- all_genes[all_genes %in% UC_genes]
      MC_genes_total <- all_genes[all_genes %in% MC_genes]
      
      for(cluster_name in cluster_names){
        genes_in_cluster <- cluster_assignments_high[[tissue_type]][[tumour]][[cluster_name]]
        UC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% UC_genes]
        
        MC_genes_cluster <- genes_in_cluster[genes_in_cluster %in% MC_genes]
        
        df <- cbind(c(length(UC_genes_cluster), length(MC_genes_cluster)),
                    c(length(UC_genes_total), length(MC_genes_total)))
        UC_p <- fisher.test(df, alternative="greater")$p.value
        MC_p <- fisher.test(df, alternative="less")$p.value
        
        temp[[tissue_type]] <- rbind(temp[[tissue_type]],
                                          c(tumour, tissue_type, cluster_name, UC_p, MC_p))
      }
      colnames(temp[[tissue_type]]) <- c("tumour", "tissue_type", "cluster", "p_UC_enrichment", "p_MC_enrichment")
      temp[[tissue_type]] <- as.data.frame(temp[[tissue_type]])
      temp[[tissue_type]]$p_UC_enrichment <- as.numeric(as.character(temp[[tissue_type]]$p_UC_enrichment))
      temp[[tissue_type]]$p_MC_enrichment <- as.numeric(as.character(temp[[tissue_type]]$p_MC_enrichment))
      temp[[tissue_type]]$p_UC_adj <- p.adjust(temp[[tissue_type]]$p_UC_enrichment, method="BH")
      temp[[tissue_type]]$p_MC_adj <- p.adjust(temp[[tissue_type]]$p_MC_enrichment, method="BH")
      temp[[tissue_type]]$Module_age <- ifelse(temp[[tissue_type]]$p_UC_adj < 0.05, "UC", 
                                               ifelse(temp[[tissue_type]]$p_MC_adj < 0.05, "MC",
                                                           "Mixed"))
    }
  }
  for(tissue_type in names(temp)){
    age_enrichment[[tumour]] <- rbind(age_enrichment[[tumour]],
                                      temp[[tissue_type]])
    age_enrichment[[tumour]]$p_UC_enrichment <- as.numeric(as.character(age_enrichment[[tumour]]$p_UC_enrichment))
    age_enrichment[[tumour]]$p_MC_enrichment <- as.numeric(as.character(age_enrichment[[tumour]]$p_MC_enrichment))
    age_enrichment[[tumour]]$p_UC_adj <- as.numeric(as.character(age_enrichment[[tumour]]$p_UC_adj))
    age_enrichment[[tumour]]$p_MC_adj <- as.numeric(as.character(age_enrichment[[tumour]]$p_MC_adj))
  }
}
save(age_enrichment, file="Objects/age_enrichment_high.Rdata")
