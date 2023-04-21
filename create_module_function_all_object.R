library(gprofiler2)

load("cluster_assignments.Rdata")

module_function <- vector()
for(tissue in c("normal", "tumour")){
  for(tumour in names(cluster_assignments[[tissue]])){
    
    for(mod in names(cluster_assignments[[tissue]][[tumour]])){
      gene_set <- cluster_assignments[[tissue]][[tumour]][[mod]]    
      
      result <- gost(gene_set, organism = "hsapiens", ordered_query = FALSE,
                     multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                     measure_underrepresentation = FALSE, evcodes = FALSE,
                     user_threshold = 0.05, correction_method = c("g_SCS", "bonferroni",
                                                                  "fdr", "false_discovery_rate", "gSCS", "analytical"),
                     domain_scope = c("annotated", "known", "custom"), custom_bg = NULL,
                     numeric_ns = "", sources = c("KEGG", "REAC"))
      
      if(is.null(result)){
        module_function <- rbind(module_function, cbind(tissue, tumour, mod, source="empty", term_name="empty",
                                                        p_value="empty"))
      }else{
        result <- result$result
        result <- result[order(result$p_value, decreasing=FALSE),]
        if(nrow(result) == 1){
          sig <- result[1,c("source", "term_name", "p_value")]
        }else{
          sig <- result[1:2,c("source", "term_name", "p_value")]
        }
        module_function <- rbind(module_function, cbind(tissue, tumour, mod, sig))
      }
      print(mod)
    }
    print(tumour)
    #save(module_function, file=paste("module_function_", type, ".Rdata", sep=""))
  }
}

module_function_all <- vector()
for(tissue in c("normal", "tumour")){
  for(tumour in names(cluster_assignments[[tissue]])){
    load(paste("module_function_", tumour, "_", tissue, ".Rdata", sep=""))
    module_function_all <- rbind(module_function_all, module_function)
  }
}
colnames(module_function_all) <- c("Tissue", "Tumour", "Module", "Pathway_database", "Pathway", "p_value")
module_function_all <- as.data.frame(module_function_all)
module_function_all$Module_name <- paste(module_function_all$Tumour, module_function_all$Tissue, module_function_all$Module, sep="_")

load("age_enrichment.Rdata")
load("all_preservation_t_to_n2.Rdata")
load("all_preservation_n_to_t2.Rdata")

module_function_all2 <- vector()
for(tumour in names(age_enrichment)){
  temp <- age_enrichment[[tumour]]
  temp$Module_name <- paste(temp$tumour, temp$tissue_type, temp$cluster, sep="_")
  
  temp2 <- module_function_all[module_function_all$Tumour == tumour,]
  temp2$Module_age <- temp$Module_age[match(temp2$Module_name, temp$Module_name)]
  
  temp2_normal <- temp2[temp2$Tissue == "normal",]
  temp2_tumour <- temp2[temp2$Tissue == "tumour",]
  
  if(nrow(temp2_normal) > 0){
    temp2_normal$Novelty <- all_preservation_n_to_t2$Category[match(temp2_normal$Module_name, all_preservation_n_to_t2$Cluster_name)]
  }
  temp2_tumour$Novelty <- all_preservation_t_to_n2$Category[match(temp2_tumour$Module_name, all_preservation_t_to_n2$Cluster_name)]
  
  module_function_all2 <- rbind(module_function_all2, temp2_tumour)
  if(nrow(temp2_normal) > 0){
    module_function_all2 <- rbind(module_function_all2, temp2_normal)
  }
}
module_function_all <- module_function_all2
colnames(module_function_all)

#save(module_function_all, file="module_function_all.Rdata")