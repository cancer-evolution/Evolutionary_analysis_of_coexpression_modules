
add_number_cc_genes <- function(preservation_age, tumour, cancer_genes_per_tumour){
  local_preservation_age <- preservation_age[preservation_age$Tumour == tumour,]
  
  local_preservation_age$N_cancer_census <- sapply(as.character(local_preservation_age$Cluster), function(cluster){
    sum(cluster_assignments[["tumour"]][[tumour]][[cluster]] %in% cancer_genes_per_tumour[[tumour]])
  })
  local_preservation_age$Per_cancer_census <- sapply(as.character(local_preservation_age$Cluster), function(cluster){
    (sum(cluster_assignments[["tumour"]][[tumour]][[cluster]] %in% cancer_genes_per_tumour[[tumour]])/length(cluster_assignments[["tumour"]][[tumour]][[cluster]]))*100
  })
  return(local_preservation_age)
}

add_percentage_UC <- function(all_preservation, tumour, direction){
  local_all_preservation <- all_preservation[all_preservation$Tumour == tumour,]
  if(direction == "NormalToTumour"){
    temp_genes <- cluster_assignments$normal[[tumour]]
  }else{
    temp_genes <- cluster_assignments$tumour[[tumour]]
  }
  
  length_temp_genes <- sapply(temp_genes, length)
  
  temp_genes <- sapply(temp_genes, function(set){
    sum(set %in% UC_genes)/sum(set %in% c(UC_genes, MC_genes))*100
  })
  
  per <- temp_genes-expected_per_UC
  
  local_all_preservation$Diff_per_UC <- per[match(local_all_preservation$Cluster, names(per))]
  local_all_preservation$Cluster_size <- unname(length_temp_genes[match(local_all_preservation$Cluster, names(length_temp_genes))])
  
  return(local_all_preservation)
}

calculate_association_preservation_age_categorical <- function(preservation,age_enrichment,tumour, direction){
  all_preservation <- vector()
  UC_inf_all <- 0
  Mixed_inf_all <- 0
  MC_inf_all <- 0
  local_preservation <- preservation[preservation$Tumour == tumour, c("Tumour", "Cluster", "Per_50")]
  if(direction == "TumourToNormal"){
    local_age_enrichment <- age_enrichment[[tumour]][age_enrichment[[tumour]]$tissue_type=="tumour",]  #FOR TUMOUR MODULES
  }else{
    local_age_enrichment <- age_enrichment[[tumour]][age_enrichment[[tumour]]$tissue_type=="normal",]
  }
  local_preservation$Age <- local_age_enrichment[match(local_preservation$Cluster, local_age_enrichment$cluster), "Module_age"]
  local_preservation$Per_50 <- as.numeric(as.character(local_preservation$Per_50))

  local_preservation$Age <- factor(local_preservation$Age, levels=c("UC", "Mixed", "MC"))
  all_preservation <- rbind(all_preservation, local_preservation)
  UC_inf <- sum(is.infinite(local_preservation[local_preservation$Age == "UC", "Per_50"]))
  MC_inf <- sum(is.infinite(local_preservation[local_preservation$Age == "MC", "Per_50"]))
  Mixed_inf <- sum(is.infinite(local_preservation[local_preservation$Age == "Mixed", "Per_50"]))
  UC_inf_all <- UC_inf_all+UC_inf
  MC_inf_all <- MC_inf_all+MC_inf
  Mixed_inf_all <- Mixed_inf_all+Mixed_inf

  g <- ggplot(local_preservation, aes(x=Age, y=Per_50))+
    geom_boxplot(aes(fill=Age))+
    ggtitle(tumour)+
    ylab("Preservation")+
    annotate("text", x = 1:3, y = 30, label = c(UC_inf, Mixed_inf, MC_inf))+
    theme_bw()
  return(list(plot=g, preservation = all_preservation))
}

calculate_degree_in_modules <- function(subnet_local, tumour, tissue_type){
  local_degree_df <- vector()
  for(mod in names(subnet_local)){
    local_net <- subnet_local[[mod]]
    local_net <- replace(local_net, local_net == 1, NA)
    local_degree <- apply(local_net, 1, sum, na.rm=TRUE)
    mean_degree <- median(local_degree, na.rm=TRUE)
    max_degree <- max(local_degree, na.rm=TRUE)
    min_degree <- min(local_degree, na.rm=TRUE)
    local_degree <- data.frame(Genes=names(local_degree),
                               Degree=local_degree,
                               Degree_norm=local_degree/mean_degree,
                               Degree_max=local_degree/max_degree,
                               Degree_min=local_degree/min_degree,
                               Module=mod,
                               Tumour=tumour,
                               Tissue_type = tissue_type)
    local_degree_df <- rbind(local_degree_df, local_degree)
  }
  local_degree_df <- as.data.frame(local_degree_df)
  local_degree_df$Degree <- as.numeric(as.character(local_degree_df$Degree))
  return(local_degree_df)
}

calculate_degree_rank <- function(subnet_local, tissue_type, tumour){
  gene_degree_df <- vector()
  for(mod in names(subnet_local)){
    subnet_local_mod <- subnet_local[[mod]]
    subnet_local_mod[subnet_local_mod ==1] <- NA
    gene_degree <- apply(subnet_local_mod, 1, sum, na.rm=TRUE)
    gene_degree <- data.frame(Gene=names(gene_degree),
                              Degree=gene_degree,
                              Degree_norm=gene_degree/median(gene_degree),
                              Degree_rank=rank(-gene_degree),
                              Degree_rank_norm=rank(-gene_degree)/length(gene_degree),
                              Module=mod,
                              Tissue_type=tissue_type,
                              Tumour=tumour)
    gene_degree_df <- rbind(gene_degree_df, gene_degree)
  }
  return(gene_degree_df)
}

calculate_number_of_connections <- function(subnet_local, tumour, tissue_type){
  number <- vector()
  for(mod in names(subnet_local)){
    local_subnet <- subnet_local[[mod]]
    
    local_subnet[upper.tri(local_subnet, diag=TRUE)] <- NA
    local_subnet <- melt(local_subnet)
    local_subnet <- local_subnet[!is.na(local_subnet$value),]
    
    local_subnet$Age1 <- ifelse(local_subnet$Var1 %in% UC_genes, "UC",
                                ifelse(local_subnet$Var1 %in% MC_genes, "MC", NA))
    local_subnet$Age2 <- ifelse(local_subnet$Var2 %in% UC_genes, "UC",
                                ifelse(local_subnet$Var2 %in% MC_genes, "MC", NA))
    local_subnet$Conn_age <- paste(local_subnet$Age1, local_subnet$Age2, sep="_")
    
    UC_UC_connections <- sum(local_subnet$Conn_age == "UC_UC")
    MC_MC_connections <- sum(local_subnet$Conn_age == "MC_MC")
    UC_MC_connections <- sum(local_subnet$Conn_age %in% c("UC_MC", "MC_UC"))
    
    total <- UC_UC_connections+MC_MC_connections+UC_MC_connections
    
    UC_UC_connections <- (UC_UC_connections/total)*100
    UC_MC_connections <- (UC_MC_connections/total)*100
    MC_MC_connections <- (MC_MC_connections/total)*100
    
    number <- rbind(number, c(tumour, tissue_type, mod, UC_UC_connections, UC_MC_connections, MC_MC_connections))
  }
  colnames(number) <- c("Tumour", "Tissue_type", "Module", "UC_UC_connections", "UC_MC_connections", "MC_MC_connections")
  number <- as.data.frame(number)
  
  number$UC_UC_connections <- as.numeric(as.character(number$UC_UC_connections))
  number$UC_MC_connections <- as.numeric(as.character(number$UC_MC_connections))
  number$MC_MC_connections <- as.numeric(as.character(number$MC_MC_connections))
  
  return(number)
}

calculate_preservation <- function(direction, number_of_shared_genes, tumour){
  preservation <- vector()
  local_sample_type <-  number_of_shared_genes[intersect(which(number_of_shared_genes$Tumour1 == tumour),
                                                         which(number_of_shared_genes$Tumour2 == tumour)),]
  local_sample_type <- local_sample_type[intersect(which(local_sample_type$Tissue1 == "normal"),
                                                   which(local_sample_type$Tissue2 == "tumour")),]
  if(direction == "TumourToNormal"){
    clusters_local_tissue_type <- unique(local_sample_type$Cluster2) ###tumour clusters
  }else{
    clusters_local_tissue_type <- unique(local_sample_type$Cluster1) ###normal clusters
  }
  
  for(local_cluster in clusters_local_tissue_type){
    if(direction == "TumourToNormal"){
      local_sample_type2 <- local_sample_type[which(local_sample_type$Cluster2 == local_cluster),]
    }else{
      local_sample_type2 <- local_sample_type[which(local_sample_type$Cluster1 == local_cluster),]
    }
    shared_genes <- sort(local_sample_type2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c(tumour,local_cluster, sum_to))
  }
  colnames(preservation) <- c("Tumour", "Cluster", paste("Per_", seq(5,100, 5), sep=""))
  preservation <- as.data.frame(preservation)
  return(preservation)
}

calculate_preservation2 <- function(direction, number_of_shared_genes, tumour){
  preservation <- vector()
  local_sample_type <-  number_of_shared_genes[intersect(which(number_of_shared_genes$Tumour1 == tumour),
                                                         which(number_of_shared_genes$Tumour2 == tumour)),]
  
  if(direction == "TumourToNormal"){
    local_sample_type <- local_sample_type[intersect(which(local_sample_type$Tissue1 == "tumour"),
                                                     which(local_sample_type$Tissue2 == "normal")),]
    clusters_local_tissue_type <- unique(local_sample_type$Cluster1) ###tumour clusters
  }else{
    local_sample_type <- local_sample_type[intersect(which(local_sample_type$Tissue1 == "normal"),
                                                     which(local_sample_type$Tissue2 == "tumour")),]
    clusters_local_tissue_type <- unique(local_sample_type$Cluster1) ###normal clusters
  }
  
  for(local_cluster in clusters_local_tissue_type){
    local_sample_type2 <- local_sample_type[local_sample_type$Cluster1 == local_cluster,]
    shared_genes <- sort(local_sample_type2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c(tumour,local_cluster, sum_to))
  }
  colnames(preservation) <- c("Tumour", "Cluster", paste("Per_", seq(5,100, 5), sep=""))
  preservation <- as.data.frame(preservation)
  return(preservation)
}

calculate_preservation_melanoma <- function(number_of_shared_genes){
  preservation <- vector()
  local_sample_type_p_n <-  number_of_shared_genes
  
  local_sample_type_p_n <- local_sample_type_p_n[intersect(which(local_sample_type_p_n$Type1 == "primary"),
                                                           which(local_sample_type_p_n$Type2 == "nevus")),]
  clusters_local_tissue_type_p_n <- unique(local_sample_type_p_n$Module1)
  
  for(local_cluster in clusters_local_tissue_type_p_n){
    local_sample_type_p_n2 <- local_sample_type_p_n[local_sample_type_p_n$Module1 == local_cluster,]
    shared_genes <- sort(local_sample_type_p_n2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c("primary", local_cluster, sum_to))
  }
  
  colnames(preservation) <- c("Type", "Cluster", paste("Per_", seq(5,100, 5), sep=""))
  preservation <- as.data.frame(preservation)
  return(preservation)
}

calculate_preservation_pheo <- function(number_of_shared_genes){
  preservation <- vector()
  local_sample_type_m_b <-  number_of_shared_genes
  
  local_sample_type_m_b <- local_sample_type_m_b[intersect(which(local_sample_type_m_b$Type1 == "malignant"),
                                                           which(local_sample_type_m_b$Type2 == "benign")),]
  clusters_local_tissue_type_m_b <- unique(local_sample_type_m_b$Module1)
  
  for(local_cluster in clusters_local_tissue_type_m_b){
    local_sample_type_m_b2 <- local_sample_type_m_b[local_sample_type_m_b$Module1 == local_cluster,]
    shared_genes <- sort(local_sample_type_m_b2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c("malignant", local_cluster, sum_to))
  }
  
  local_sample_type_b_n <-  number_of_shared_genes
  
  local_sample_type_b_n <- local_sample_type_b_n[intersect(which(local_sample_type_b_n$Type1 == "benign"),
                                                           which(local_sample_type_b_n$Type2 == "normal")),]
  clusters_local_tissue_type_b_n <- unique(local_sample_type_b_n$Module1)
  
  for(local_cluster in clusters_local_tissue_type_b_n){
    local_sample_type_b_n2 <- local_sample_type_b_n[local_sample_type_b_n$Module1 == local_cluster,]
    shared_genes <- sort(local_sample_type_b_n2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c("benign", local_cluster, sum_to))
  }
  
  colnames(preservation) <- c("Type", "Cluster", paste("Per_", seq(5,100, 5), sep=""))
  preservation <- as.data.frame(preservation)
  return(preservation)
}

calculate_preservation_PRAD <- function(number_of_shared_genes){
  preservation <- vector()
  local_sample_type_m_b <-  number_of_shared_genes
  
  local_sample_type_m_b <- local_sample_type_m_b[intersect(which(local_sample_type_m_b$Type1 == "high_grade"),
                                                           which(local_sample_type_m_b$Type2 == "low_grade")),]
  clusters_local_tissue_type_m_b <- unique(local_sample_type_m_b$Module1)
  
  for(local_cluster in clusters_local_tissue_type_m_b){
    local_sample_type_m_b2 <- local_sample_type_m_b[local_sample_type_m_b$Module1 == local_cluster,]
    shared_genes <- sort(local_sample_type_m_b2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c("high_grade", local_cluster, sum_to))
  }
  
  local_sample_type_b_n <-  number_of_shared_genes
  
  local_sample_type_b_n <- local_sample_type_b_n[intersect(which(local_sample_type_b_n$Type1 == "low_grade"),
                                                           which(local_sample_type_b_n$Type2 == "normal")),]
  clusters_local_tissue_type_b_n <- unique(local_sample_type_b_n$Module1)
  
  for(local_cluster in clusters_local_tissue_type_b_n){
    local_sample_type_b_n2 <- local_sample_type_b_n[local_sample_type_b_n$Module1 == local_cluster,]
    shared_genes <- sort(local_sample_type_b_n2$Shared_genes, decreasing=TRUE)
    cum_sum <- cumsum(shared_genes)
    sum_to <- vector()
    for(i in seq(5,100, 5)){
      sum_to <- c(sum_to, min(which(cum_sum > i)))
    }
    preservation <- rbind(preservation, c("low_grade", local_cluster, sum_to))
  }
  
  colnames(preservation) <- c("Type", "Cluster", paste("Per_", seq(5,100, 5), sep=""))
  preservation <- as.data.frame(preservation)
  return(preservation)
}

calculate_strength_of_subnetworks <- function(subnet_local, tumour, tissue_type){
  strength_summary <- vector()
  for(mod in names(subnet_local)){
    local_subnet <- subnet_local[[mod]]
    strength_summary_for_module <- vector()
    for(con in c("All", "UC", "MC", "Mixed")){
      if(con != "Mixed"){
        if(con == "All"){
          genes <- rownames(local_subnet)
        }else if(con == "UC"){
          genes <- rownames(local_subnet)[rownames(local_subnet) %in% UC_genes]
        }else if(con == "MC"){
          genes <- rownames(local_subnet)[rownames(local_subnet) %in% MC_genes]
        }
        local_subnet2 <- local_subnet[rownames(local_subnet) %in% genes,
                                      colnames(local_subnet) %in% genes]
        total_strength <- sum(local_subnet2)
        number_connections <- nrow(local_subnet2)*ncol(local_subnet2)
        ave_strength <- total_strength/number_connections
      }else{
        local_subnet2 <- local_subnet[rownames(local_subnet) %in% UC_genes,
                                      colnames(local_subnet) %in% MC_genes]
        total_strength <- sum(local_subnet2)
        number_connections <- sum(rownames(local_subnet) %in% UC_genes)*sum(rownames(local_subnet) %in% MC_genes)
        ave_strength <- total_strength/number_connections
      }
      strength_summary_for_module <- rbind(strength_summary_for_module, 
                                           c(tumour, tissue_type, mod, total_strength, number_connections, ave_strength, con))
      
    }
    
    colnames(strength_summary_for_module) <- c("Tumour", "Tissue_type", "Module", "Strength", "Number_connections",
                                               "Average_strength", "Connection_type")
    
    strength_summary_for_module <- as.data.frame(strength_summary_for_module)
    
    strength_summary_for_module$Strength <- as.numeric(as.character(strength_summary_for_module$Strength))
    strength_summary_for_module$Number_connections <- as.numeric(as.character(strength_summary_for_module$Number_connections))
    
    total_st <- strength_summary_for_module$Strength[1]
    strength_summary_for_module$Per_exp <- (strength_summary_for_module$Strength/total_st)*100
    
    strength_summary <- rbind(strength_summary,
                              strength_summary_for_module)
    
  }
  return(strength_summary)
}

define_categories_preservation <- function(all_preservation, tumours){
  result <- vector()
  for(tumour in tumours){
    local_subset <- subset(all_preservation, Tumour==tumour)
    temp <- local_subset$Preservation_ratio[is.finite(local_subset$Preservation_ratio)]
    cutoffs <- quantile(temp, c(1/3, 2/3))
    local_subset$Category <- NA
    local_subset[local_subset$Preservation_ratio < cutoffs[1],"Category"] <- "Low_score"
    local_subset[intersect(which(local_subset$Preservation_ratio >= cutoffs[1]),
                           which(local_subset$Preservation_ratio < cutoffs[2])), "Category"] <- "Median_score"
    local_subset[local_subset$Preservation_ratio >= cutoffs[2], "Category"] <- "High_score" 
    local_subset[is.infinite(local_subset$Preservation_ratio), "Category"] <- "Inf_score" 
    result <- rbind(result, local_subset)
  }
  result$Category <- factor(result$Category,
                            levels=c("Low_score", "Median_score", "High_score", "Inf_score"))
  
  return(result)
}

define_categories_preservation_pheo <- function(all_preservation){
  result <- vector()
  local_subset <- all_preservation
  temp <- local_subset$Preservation_ratio[is.finite(local_subset$Preservation_ratio)]
  cutoffs <- quantile(temp, c(1/3, 2/3))
  local_subset$Category <- NA
  local_subset[local_subset$Preservation_ratio < cutoffs[1],"Category"] <- "Low_score"
  local_subset[intersect(which(local_subset$Preservation_ratio >= cutoffs[1]),
                         which(local_subset$Preservation_ratio < cutoffs[2])), "Category"] <- "Medium_score"
  local_subset[local_subset$Preservation_ratio >= cutoffs[2], "Category"] <- "High_score" 
  local_subset[is.infinite(local_subset$Preservation_ratio), "Category"] <- "Inf_score" 
  result <- rbind(result, local_subset)
  
  result$Category <- factor(result$Category,
                            levels=c("Low_score", "Medium_score", "High_score", "Inf_score"))
  return(result)
}

get_metrics_subnetworks <- function(subnet_local, tissue_type, tumour){
  subnet_metrics <- vector()
  for(mod in names(subnet_local)){
    local_subnet <- melt(subnet_local[[mod]])
    colnames(local_subnet) <- c("Gene1", "Gene2", "Strength")
    local_subnet <- data.frame(local_subnet,
                               Module=paste(tissue_type, mod, sep="_"),
                               Tumour = tumour)
    local_subnet <- local_subnet[local_subnet$Strength != 1,]  #with itself
    temp1 <- local_subnet[,c(1,3)]
    temp2 <- local_subnet[,c(2,3)]
    colnames(temp1)[1] <- "Gene"
    colnames(temp2)[1] <- "Gene"
    
    local_subnet_long <- rbind(temp1, temp2)
    
    gene_attributes <- aggregate(Strength~Gene, local_subnet_long, sum)
    colnames(gene_attributes)[2] <- "Degree"
    
    gene_attributes$Age <- genes_phy_categorical[match(gene_attributes$Gene, genes_phy_categorical$GeneID), "Phylostrata"]  
    
    gene_attributes$Mut <- NA
    gene_attributes$Mut[gene_attributes$Gene %in% genes_lof] <- "LoF"
    gene_attributes$Mut[gene_attributes$Gene %in% genes_miss] <- "Miss"
    
    gene_attributes$CNV <- NA
    gene_attributes$CNV[gene_attributes$Gene %in% genes_amp] <- "Amp"
    gene_attributes$CNV[gene_attributes$Gene %in% genes_del] <- "Del"
    
    gene_attributes$Age <- factor(gene_attributes$Age, levels=c("UC", "EM", "MM"))
    
    gene_attributes$Module <- mod
    gene_attributes$Tissue_type <- tissue_type
    
    subnet_metrics <- rbind(subnet_metrics, gene_attributes)
    print(mod)
  }
  
  return(subnet_metrics)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

print_properties_modules <- function(all_preservation, local_modules,
                                     local_clusters_attributes, tumour){
  
  local_module_properties <- vector()
  for(mod in names(local_modules)){
    module_genes <- local_modules[[mod]]
    per_miss <- (sum(module_genes %in% missense)/length(module_genes))*100
    per_lof <- (sum(module_genes %in% lof)/length(module_genes))*100
    per_amp <- (sum(module_genes %in% amp)/length(module_genes))*100
    per_del <- (sum(module_genes %in% del)/length(module_genes))*100
    per <- c(per_miss, per_lof, per_amp, per_del)
    local_clusters <- local_clusters_attributes[local_clusters_attributes$Module == mod,]
    n_cluster_less_10 <- sum(local_clusters$Size < 10)
    n_cluster_more_10 <- sum(local_clusters$Size >= 10)
    age_pre <- all_preservation[all_preservation$Cluster == mod, c("Age", "Category")]
    age_pre <- unlist(age_pre)
    local_module_properties <- rbind(local_module_properties,
                                     c(mod, per, n_cluster_less_10, n_cluster_more_10, as.character(age_pre)))
  }
  local_module_properties <- data.frame(local_module_properties)
  colnames(local_module_properties) <- c("Module", "Per_miss", "Per_lof", "Per_amp", "Per_del", "N_cluster_less_10",
                                         "N_cluster_more_10", "Age", "Preservation")
  
  
  ##Functional enrichment of modules
  
  enrichment_modules <- vector()
  for(mod in names(local_modules)){
    module_genes <- local_modules[[mod]]
    
    n_module_genes <- length(module_genes)
    
    significance_module <- sapply(1:length(all_go), function(i){
      local_go <- all_go[i]
      genes_go_local <- unique(gene_GO_annotations_full[[local_go]])
      n_local_go <- length(genes_go_local)
      
      n_local_go_and_module <- sum(module_genes %in% genes_go_local)
      
      n_module_genes_not_in_go <- n_module_genes-n_local_go_and_module
      
      n_other_go_in_module <- sum(module_genes %in% all_go_genes)
      n_other_go_nomodule <- length(all_go_genes) - n_module_genes_not_in_go
      
      p <- fisher.test(cbind(c(n_local_go_and_module, n_module_genes_not_in_go), c(n_local_go, n_go_genes-n_local_go)), alternative = "greater")$p.value
      return(p)
    })
    enrichment_modules <- rbind(enrichment_modules,
                                cbind(Tumour=tumour, Module=mod, GO=all_go, Enrichment_p=significance_module))
    print(mod)
  }
  
  enrichment_modules_df <- as.data.frame(enrichment_modules)
  enrichment_modules_df$Enrichment_p <- as.numeric(as.character(enrichment_modules_df$Enrichment_p))
  enrichment_modules_df$Enrichment_p_adj <- p.adjust(  enrichment_modules_df$Enrichment_p, method="BH")
  enrichment_modules_df$GO_id <- names(gonames)[match(enrichment_modules_df$GO, gonames)]
  
  
  sig_enrichment_modules_df <- enrichment_modules_df[enrichment_modules_df$Enrichment_p < 0.05,]  #note that these are nominal p-values
  
  table(sig_enrichment_modules_df$Module)
  
  go_summary_module <- vector()
  for(module in unique(sig_enrichment_modules_df$Module)){
    local_enrich <- sig_enrichment_modules_df[sig_enrichment_modules_df$Module == module,]
    summary <- summarize_go_enrichment_results_to_top2(local_enrich$GO_id)
    local_max <- max(summary)
    go_summary_module <- rbind(go_summary_module,
                               c(module, paste(names(summary[summary==local_max]), collapse=",")))
  }
  
  sig_enrichment_modules_df$Module_GO <- go_summary_module[match(sig_enrichment_modules_df$Module, go_summary_module[,1]),2]
  
  summary_significance <- unique(sig_enrichment_modules_df[,c("Module", "Module_GO")])
  
  local_module_properties$GO_enrichment <- summary_significance[match(local_module_properties$Module, summary_significance$Module), "Module_GO"]
  return(local_module_properties)
}

load_CNVs <- function(only_focal){
  CNVs_df <- vector()
  for(tumour in tumours){
    local_amp <- CNVs_curated[[tumour]]$amplifications
    local_del <- CNVs_curated[[tumour]]$deletions
    
    local_amp$Sample <- substr(local_amp$Sample, 9, 12)
    local_del$Sample <- substr(local_del$Sample, 9, 12)
    
    if(only_focal == "FOCAL"){
      focal_amp <- as.character(CNVs_above_fraction_0.25[[tumour]][["Amp"]])
      focal_del <- as.character(CNVs_above_fraction_0.25[[tumour]][["Del"]])
      
      local_amp <- local_amp[local_amp$Gene %in% focal_amp,]
      local_del <- local_del[local_del$Gene %in% focal_del,]
    }else if(only_focal == "GLOBAL"){
      
      focal_amp <- as.character(CNVs_above_fraction_0.25[[tumour]][["Amp"]])
      focal_del <- as.character(CNVs_above_fraction_0.25[[tumour]][["Del"]])
      
      local_amp <- local_amp[!(local_amp$Gene %in% focal_amp),]
      local_del <- local_del[!(local_del$Gene %in% focal_del),]
    }else if(only_focal=="ANY"){
      #No subsetting
    }
    
    n_pat_amp <- aggregate(Segment_Mean ~ Gene, local_amp, length)
    colnames(n_pat_amp) <- c("Gene", "N_patients_amp")
    
    n_pat_del <- aggregate(Segment_Mean ~ Gene, local_del, length)
    colnames(n_pat_del) <- c("Gene", "N_patients_del")
    
    all_genes <- unique(c(n_pat_amp$Gene, n_pat_del$Gene))
    
    n_pat_df <- data.frame(Genes = all_genes,
                           Patients_amp =n_pat_amp[match(all_genes, n_pat_amp$Gene), "N_patients_amp"],
                           Patients_del =n_pat_del[match(all_genes, n_pat_del$Gene), "N_patients_del"])
    
    local_patients <- length(patients_with_CNV_info[[tumour]])
    n_pat_df$Patients_amp <- n_pat_df$Patients_amp/local_patients
    n_pat_df$Patients_del <- n_pat_df$Patients_del/local_patients
    
    n_pat_df$Patients_amp[is.na(n_pat_df$Patients_amp)] <- 0
    n_pat_df$Patients_del[is.na(n_pat_df$Patients_del)] <- 0
    
    n_pat_df$Gene_age <- genes_phy[match(n_pat_df$Genes, genes_phy$GeneID), "Phylostrata"]
    n_pat_df$Gene_age <- factor(n_pat_df$Gene_age, levels=1:16)
    n_pat_df$Tumour <- tumour
    CNVs_df <- rbind(CNVs_df, n_pat_df)
  }
  return(CNVs_df)
}

load_mutations <- function(){
  load("gene_mut_properties.Rdata")
  mut_to_exclude <- read_csv("Genes to exclude.txt", 
                             col_names = FALSE)
  mut_to_exclude <- unname(unlist(mut_to_exclude[,1]))
  
  mutations_df <- vector()
  for(tumour in tumours){
    local_mut <- gene_mut_properties[[tumour]]
    local_mut <- local_mut[!(local_mut$Hugo_Symbol %in% mut_to_exclude),]
    
    local_mut$Replication_time <- NULL
    local_mut$HiC <- NULL
    local_mut$Gene_length <- NULL
    local_mut$GC_percentage <- NULL
    local_mut$LoF_MB <- NULL
    local_mut$Syn_MB <- NULL
    
    local_mut$Gene_age <- genes_phy[match(local_mut$Hugo_Symbol, genes_phy$GeneID), "Phylostrata"]
    local_mut$Tumour <- tumour
    mutations_df <- rbind(mutations_df, local_mut)
    
  }
  return(mutations_df)
}

read_expression_paired <- function(tumour_type){
  path <- paste("TCGA_all_data/", tumour, "/gdac.broadinstitute.org_", tumour, 
                ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/", 
                tumour, 
                ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                sep="")
  
  exp <- read.delim(path) ##Load normalized gene expression data from TCGA
  exp <- exp[-1,]
  
  exp_tumour <- exp[,substr(colnames(exp), 14, 16) %in% c("01A", "01B", "01C", "01D")]
  exp_normal <- exp[,substr(colnames(exp), 14, 16) %in% c("11A", "11B", "11C", "11D")]
  
  all_genes <- as.character(exp[,1])
  
  all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]
  
  genes_remove <- c("?", "SLC35E2")
  genes_remove_index <- which(all_genes %in% genes_remove)
  all_genes <- all_genes[-genes_remove_index]
  
  ##If there are normal samples, keep matched samples
  if(!is.null(ncol(exp_normal))){
    if(ncol(exp_normal) >= 10){
      
      colnames(exp_normal) <- substr(colnames(exp_normal),9,12)
      exp_normal <- exp_normal[-genes_remove_index,]
      colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
      exp_tumour <- exp_tumour[-genes_remove_index,]
      
      rownames(exp_normal) <- all_genes
      rownames(exp_tumour) <- all_genes
      
      patients_normal <- colnames(exp_normal)
      patients_tumour <- colnames(exp_tumour)
      
      patients_paired <- patients_normal[patients_normal %in% patients_tumour]
      exp_normal <- exp_normal[,colnames(exp_normal) %in% patients_paired]
      exp_tumour <- exp_tumour[,colnames(exp_tumour) %in% patients_paired]
      
      exp_normal <- apply(exp_normal, 2, function(x){as.numeric(as.character(x))})
      exp_tumour <- apply(exp_tumour, 2, function(x){as.numeric(as.character(x))})
      
      genes_reasonable_expression_normal <- which(rowSums(cpm(exp_normal)>1) == ncol(exp_normal))
      exp_normal <- exp_normal[genes_reasonable_expression_normal,]
      
      genes_reasonable_expression_tumour <- which(rowSums(cpm(exp_tumour)>1) == ncol(exp_tumour))
      exp_tumour <- exp_tumour[genes_reasonable_expression_tumour,]
      rownames(exp_normal) <- all_genes[genes_reasonable_expression_normal]
      rownames(exp_tumour) <- all_genes[genes_reasonable_expression_tumour]
      
    }else{
      colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
      exp_tumour <- exp_tumour[-genes_remove_index,]
      rownames(exp_tumour) <- all_genes
      
      exp_tumour <- apply(exp_tumour, 2, function(x){as.numeric(as.character(x))})
      
      genes_reasonable_expression_tumour <- which(rowSums(cpm(exp_tumour)>1) == ncol(exp_tumour))
      exp_tumour <- exp_tumour[genes_reasonable_expression_tumour,]
      
      exp_normal <- NA
      rownames(exp_tumour) <- all_genes[genes_reasonable_expression_tumour]
    }
  }else{
    
    colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
    exp_tumour <- exp_tumour[-genes_remove_index,]
    rownames(exp_tumour) <- all_genes
    
    exp_tumour <- apply(exp_tumour, 2, function(x){as.numeric(as.character(x))})
    
    genes_reasonable_expression_tumour <- which(rowSums(cpm(exp_tumour)>1) == ncol(exp_tumour))
    exp_tumour <- exp_tumour[genes_reasonable_expression_tumour,]
    
    exp_normal <- NA
    rownames(exp_tumour) <- all_genes[genes_reasonable_expression_tumour]
    
  }
  return(list(normal = exp_normal, tumour=exp_tumour))
}

remove_outlier_samples <- function(local_exp, local_name, cutHeight, no_clusters_keep){
  sampleTree = hclust(dist(local_exp), method = "average");
  plot(sampleTree, main = paste("Sample clustering to detect outliers", local_name, sep="\n"), sub="", xlab="")
  abline(h = cutHeight, col = "red");
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
  keepSamples = (clust %in% no_clusters_keep)
  final_exp = local_exp[keepSamples, ]
  print(dim(final_exp))
  return(final_exp)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)