library(igraph)
library(reshape2)

source("functions.R")

tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

tumours_normal <- c("BLCA", "BRCA","COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
                    "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA", "UCEC")


load("cluster_assignments.Rdata")


for(tumour1 in tumours){
  number_of_shared_genes <- vector()
  for(tissue_type1 in c("normal", "tumour")){
    if(tissue_type1 == "tumour" || (tumour1 %in% tumours_normal & tissue_type1 == "normal")){
      clusters_of_local_sample_type <- cluster_assignments[[tissue_type1]][[tumour1]]
      for(local_name1 in names(clusters_of_local_sample_type)){
        genes_in_local_cluster <- clusters_of_local_sample_type[[local_name1]]
        for(tumour2 in tumours){
          for(tissue_type2 in c("normal", "tumour")){
            if(tissue_type2 == "tumour" || (tumour2 %in% tumours_normal & tissue_type2 == "normal")){
              clusters_of_other_sample_type <- cluster_assignments[[tissue_type2]][[tumour2]]
              for(local_name2 in names(clusters_of_other_sample_type)){
                genes_in_other_cluster <- clusters_of_other_sample_type[[local_name2]]
                shared_genes <- sum(genes_in_local_cluster %in% genes_in_other_cluster)/length(genes_in_local_cluster)*100
                number_of_shared_genes <- rbind(number_of_shared_genes,
                                                c(local_name1, local_name2, tumour1, tissue_type1,
                                                  tumour2, tissue_type2, shared_genes))
              }
            }
          }
        }
        print(local_name1)
      }
    }
  }
  print(tumour1)
  save(number_of_shared_genes, file=paste("Number_of_shared_genes_", tumour1, ".Rdata", sep=""))
}

number_of_shared_genes_all <- vector()
for(tumour in tumours){
  load(paste("Number_of_shared_genes_", tumour, ".Rdata", sep=""))
  number_of_shared_genes_all <- rbind(number_of_shared_genes_all,
                                      number_of_shared_genes)
  print(tumour)
}


number_of_shared_genes <- as.data.frame(number_of_shared_genes_all)
colnames(number_of_shared_genes) <- c("Cluster1", "Cluster2", "Tumour_type1", "Tissue_type1", "Tumour_type2", "Tissue_type2", "Shared_genes")
number_of_shared_genes$Sample_type1 <- paste(number_of_shared_genes$Tumour_type1, number_of_shared_genes$Tissue_type1)
number_of_shared_genes$Sample_type2 <- paste(number_of_shared_genes$Tumour_type2, number_of_shared_genes$Tissue_type2)

number_of_shared_genes$Shared_genes <- as.numeric(as.character(number_of_shared_genes$Shared_genes))

save(number_of_shared_genes, file="number_of_shared_genes_paired.Rdata")



load("number_of_shared_genes_paired.Rdata")
number_of_shared_genes <- as.data.frame(number_of_shared_genes)
colnames(number_of_shared_genes) <- c("Cluster1", "Cluster2", "Tumour1", "Tissue1", "Tumour2", "Tissue2", "Shared_genes", "Sample_type1", "Sample_type2")

number_of_shared_genes$Shared_genes <- as.numeric(as.character(number_of_shared_genes$Shared_genes))


number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Cluster1 != "grey",]
number_of_shared_genes <- number_of_shared_genes[number_of_shared_genes$Cluster2 != "grey",]

preservation_n_to_t <- vector()
preservation_t_to_n <- vector()

for(tumour in tumours_normal){
  preservation_n_to_t <- rbind(preservation_n_to_t,
                               calculate_preservation2("NormalToTumour", number_of_shared_genes, tumour))
  preservation_t_to_n <- rbind(preservation_t_to_n,
                               calculate_preservation2("TumourToNormal", number_of_shared_genes, tumour))
 
}


cluster_sizes <- vector()
for(tumour in tumours){
  temp <- cluster_assignments$tumour[[tumour]]
  mod <- sapply(temp, length)
  mod <- cbind(Module=names(mod), Size=as.vector(mod), Tumour=tumour)
  cluster_sizes <- rbind(cluster_sizes, mod)
}

cluster_sizes <- as.data.frame(cluster_sizes)
cluster_sizes$Size <- as.numeric(as.character(cluster_sizes$Size))
cluster_sizes$Label <- paste(cluster_sizes$Tumour, cluster_sizes$Module)

preservation_t_to_n$Label <- paste(preservation_t_to_n$Tumour, preservation_t_to_n$Cluster)

preservation_t_to_n$Size <- cluster_sizes[match(preservation_t_to_n$Label, cluster_sizes$Label),"Size"]

preservation_t_to_n$Per_50 <- as.numeric(as.character(preservation_t_to_n$Per_50))
cor(preservation_t_to_n$Size, preservation_t_to_n$Per_50, method="sp")



#Association between module age and preservation in tumours vs. normal
load("age_enrichment.Rdata")

g_n_to_t <- list()
g_t_to_n <- list()
all_preservation_n_to_t <- vector()
all_preservation_t_to_n <- vector()

for(tumour in tumours){
  temp1 <- calculate_association_preservation_age_categorical(preservation_n_to_t,age_enrichment,tumour, "NormalToTumour")
  g_n_to_t[[tumour]] <- temp1$plot
  all_preservation_n_to_t <- rbind(all_preservation_n_to_t, temp1$preservation)
  
  temp2 <- calculate_association_preservation_age_categorical(preservation_t_to_n,age_enrichment,tumour, "TumourToNormal")
  
  g_t_to_n[[tumour]] <- temp2$plot
  all_preservation_t_to_n <- rbind(all_preservation_t_to_n, temp2$preservation)
  
}



##Measuring age of modules by percentage of UC
gene_ages <- read.csv("geneIDs_entrez_final_phylostrata_phy1_phy2_phy3_no_TCGA.txt")

genes_in_clusters <- unique(unname(unlist(cluster_assignments)))

gene_ages <- gene_ages[match(genes_in_clusters, gene_ages[,1]),]

UC_genes <- as.character(gene_ages[which(gene_ages$Phylostrata %in% 1:3),1])
MC_genes <- as.character(gene_ages[which(gene_ages$Phylostrata %in% 4:16),1])

expected_per_UC <- length(UC_genes)/(length(UC_genes)+length(MC_genes))*100

#Modules in t_to_n are tumour modules
#Modules in n_to_t are normal modules

all_preservation_t_to_n$Cluster_name <- paste(all_preservation_t_to_n$Tumour, "tumour", all_preservation_t_to_n$Cluster, sep="_")
all_preservation_n_to_t$Cluster_name <- paste(all_preservation_n_to_t$Tumour, "normal", all_preservation_n_to_t$Cluster, sep="_")

all_preservation_t_to_n2 <- vector()
all_preservation_n_to_t2 <- vector()

for(tumour in tumours_normal){
  all_preservation_t_to_n2 <- rbind(all_preservation_t_to_n2,
                                    add_percentage_UC(all_preservation_t_to_n, tumour, "TumourToNormal"))
  
  all_preservation_n_to_t2 <- rbind(all_preservation_n_to_t2,
                                    add_percentage_UC(all_preservation_n_to_t, tumour, "NormalToTumour"))
}

#Divide preservation score by module size
all_preservation_t_to_n2$Preservation_ratio <- all_preservation_t_to_n2$Per_50/all_preservation_t_to_n2$Cluster_size
all_preservation_n_to_t2$Preservation_ratio <- all_preservation_n_to_t2$Per_50/all_preservation_n_to_t2$Cluster_size

##Converting the preservation scores to categories

all_preservation_t_to_n2 <- define_categories_preservation(all_preservation_t_to_n2, tumours_normal)
all_preservation_n_to_t2 <- define_categories_preservation(all_preservation_n_to_t2, tumours_normal)

save(all_preservation_t_to_n2, file="all_preservation_t_to_n2.Rdata")
save(all_preservation_n_to_t2, file="all_preservation_n_to_t2.Rdata")

