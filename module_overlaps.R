load("Objects/age_enrichment.Rdata")
load("Objects/cluster_assignments.Rdata")
load("Objects/cluster_assignments_random.Rdata")


######TUMOUR
# Find highest overlapped module, return the percentage(over cancer1) and module name (in cancer2)
highest_overlap <- function(cancer1, cancer2,TN){
  lists1 <- cluster_assignments[[TN]][[cancer1]][names(cluster_assignments[[TN]][[cancer1]]) != "grey"]
  lists2 <- cluster_assignments[[TN]][[cancer2]][names(cluster_assignments[[TN]][[cancer2]]) != "grey"]
  overlap_percentage <- sapply(lists1, function(x){
    sapply(lists2, function(y) {
      return(length(intersect(x,y))/length(x))
    })
  })
  top_percentages <- apply(overlap_percentage, 2, max)
  return(top_percentages)
}


# Between each pair of cancer types, get the top percentages
tumour_type <- names(cluster_assignments$tumour)
all_top_overlaps <- lapply(tumour_type, function(cancer1){
  col_cancer2 <- tumour_type[tumour_type != cancer1]
  row_cluster1 <- names(cluster_assignments$tumour[[cancer1]])
  row_cluster1 <- row_cluster1[row_cluster1 != "grey"]
  cancer1_df <- data.frame(matrix(ncol = length(col_cancer2), nrow = length(row_cluster1)))
  colnames(cancer1_df) <- col_cancer2
  rownames(cancer1_df) <- row_cluster1
  
  return(cancer1_df)
})
names(all_top_overlaps) <- tumour_type

for (i in tumour_type){
  for (j in colnames(all_top_overlaps[[i]])){
    print(j)
    cancer2_overlaps <- highest_overlap(i, j, 'tumour')
    all_top_overlaps[[i]][j] <- unlist(cancer2_overlaps)
  }
}


# plot the heatmaps
library(ComplexHeatmap)
library(circlize)

#Getting ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(3)
cols

Plot_Overlap_Heatmapx <- function(cancer1) {
  temp <- age_enrichment[[cancer1]]
  temp <- temp[temp$cluster != "grey",]
  age_df <- data.frame(age = temp$Module_age[temp$tissue_type == 'tumour'])
  HA <- HeatmapAnnotation(df=age_df, which = "row",
                          col = list(age=c("UC" = "#F8766D", "Mixed" = "#00BA38", "MC" = "#619CFF")))
  col_fun = colorRamp2(c(0, 0.25, 0.75, 1), c("blue", "white", "yellow", "red"))
  print(Heatmap(as.matrix(all_top_overlaps[[cancer1]]), 
                rect_gp = gpar(col = "white", lwd = 1),
                right_annotation = HA,
                row_split = age_df,
                column_title = cancer1,
                col = col_fun))
}

pdf("Overlaps_Age_Modules_tumour.pdf", width = 10, height = 10)
for (cancer1 in names(cluster_assignments$tumour)) {
  Plot_Overlap_Heatmapx(cancer1)
}
dev.off()


############NORMAL
# Between each pair of cancer types, get the top percentages
normal_type <- names(cluster_assignments$normal)
all_top_overlaps_normal <- lapply(normal_type, function(cancer1){
  col_cancer2 <- normal_type[normal_type != cancer1]
  row_cluster1 <- names(cluster_assignments$normal[[cancer1]])
  row_cluster1 <- row_cluster1[row_cluster1 != "grey"]
  cancer1_df <- data.frame(matrix(ncol = length(col_cancer2), nrow = length(row_cluster1)))
  colnames(cancer1_df) <- col_cancer2
  rownames(cancer1_df) <- row_cluster1
  return(cancer1_df)
})
names(all_top_overlaps_normal) <- normal_type

for (i in normal_type){
  for (j in colnames(all_top_overlaps_normal[[i]])){
    print(j)
    cancer2_overlaps <- highest_overlap(i, j, 'normal')
    all_top_overlaps_normal[[i]][j] <- unlist(cancer2_overlaps)
  }
}


Plot_Overlap_Heatmapx_Normal <- function(cancer1) {
  temp <- age_enrichment[[cancer1]]
  temp <- temp[temp$cluster != "grey",]
  age_df <- data.frame(age = temp$Module_age[temp$tissue_type == 'normal'])
  HA <- HeatmapAnnotation(df=age_df, which = "row",
                          col = list(age=c("UC" = "#F8766D", "Mixed" = "#00BA38", "MC" = "#619CFF")))
  col_fun = colorRamp2(c(0, 0.25, 0.75, 1), c("blue", "white", "yellow", "red"))
  print(Heatmap(as.matrix(all_top_overlaps_normal[[cancer1]]), 
                rect_gp = gpar(col = "white", lwd = 1),
                right_annotation = HA,
                row_split = age_df,
                column_title = cancer1,
                col = col_fun))
}

pdf("Overlaps_Age_Modules_normal.pdf", width = 10, height = 10)
for (cancer1 in names(cluster_assignments$normal)) {
  Plot_Overlap_Heatmapx_Normal(cancer1)
}
dev.off()


################RANDOM

highest_overlap_random <- function(cancer1, cancer2, TN){
  lists1 <- cluster_assignments[[TN]][[cancer1]][names(cluster_assignments[[TN]][[cancer1]]) != "grey"]
  top_percentages_all <- vector()
  for(it in 1:50){
    lists2 <- cluster_assignments_random[[TN]][[cancer2]]
    reg <- paste("_", it, "$", sep="")
    lists2 <- lists2[grep(reg, names(cluster_assignments_random[[TN]][[cancer2]]))]
    lists2 <- lists2[names(lists2) != paste("grey", i, sep="_")]
    overlap_percentage <- sapply(lists1, function(x){
      sapply(lists2, function(y) {
        return(length(intersect(x,y))/length(x))
      })
    })
    top_percentages <- apply(overlap_percentage, 2, max)
    top_percentages_all <- cbind(top_percentages_all, top_percentages)
    #print(it)
  }
  colnames(top_percentages_all) <- paste("Iteration", 1:50, sep="_")
  top_percentages_mean <- apply(top_percentages_all, 1, mean)
  return(top_percentages_mean)
}


#############Tumour with random modules
tumour_type <- names(cluster_assignments$tumour)
all_top_overlaps_random <- lapply(tumour_type, function(cancer1){
  col_cancer2 <- tumour_type[tumour_type != cancer1]
  row_cluster1 <- names(cluster_assignments$tumour[[cancer1]])
  row_cluster1 <- row_cluster1[row_cluster1 != "grey"]
  cancer1_df <- data.frame(matrix(ncol = length(col_cancer2), nrow = length(row_cluster1)))
  colnames(cancer1_df) <- col_cancer2
  rownames(cancer1_df) <- row_cluster1
  
  return(cancer1_df)
})
names(all_top_overlaps_random) <- tumour_type

for (i in tumour_type){
  for (j in colnames(all_top_overlaps_random[[i]])){
    print(j)
    cancer2_overlaps <- highest_overlap_random(i, j, 'tumour')
    all_top_overlaps_random[[i]][j] <- unlist(cancer2_overlaps)
  }
  print(i)
}
save(all_top_overlaps_random, file="Objects/all_top_overlaps_random.Rdata")
load("Objects/all_top_overlaps_random.Rdata")

Plot_Overlap_Heatmapx_random <- function(cancer1) {
  temp <- age_enrichment[[cancer1]]
  temp <- temp[temp$cluster != "grey",]
  age_df <- data.frame(age = temp$Module_age[temp$tissue_type == 'tumour'])
  HA <- HeatmapAnnotation(df=age_df, which = "row",
                          col = list(age=c("UC" = "#F8766D", "Mixed" = "#00BA38", "MC" = "#619CFF")))
  col_fun = colorRamp2(c(0, 0.25, 0.75, 1), c("blue", "white", "yellow", "red"))
  print(Heatmap(as.matrix(all_top_overlaps_random[[cancer1]]), 
                rect_gp = gpar(col = "white", lwd = 1),
                right_annotation = HA,
                row_split = age_df,
                column_title = cancer1,
                col = col_fun))
}

pdf("Overlaps_Age_Modules_tumour_random.pdf", width = 10, height = 10)
for (cancer1 in names(cluster_assignments$tumour)) {
  Plot_Overlap_Heatmapx_random(cancer1)
}
dev.off()


##normal with random modules

# Between each pair of cancer types, get the top percentages
normal_type <- names(cluster_assignments$normal)
all_top_overlaps_normal_random <- lapply(normal_type, function(cancer1){
  col_cancer2 <- normal_type[normal_type != cancer1]
  row_cluster1 <- names(cluster_assignments$normal[[cancer1]])
  row_cluster1 <- row_cluster1[row_cluster1 != "grey"]
  cancer1_df <- data.frame(matrix(ncol = length(col_cancer2), nrow = length(row_cluster1)))
  colnames(cancer1_df) <- col_cancer2
  rownames(cancer1_df) <- row_cluster1
  
  return(cancer1_df)
})
names(all_top_overlaps_normal_random) <- normal_type

for (i in normal_type){
  for (j in colnames(all_top_overlaps_normal_random[[i]])){
    print(j)
    cancer2_overlaps <- highest_overlap_random(i, j, 'normal')
    all_top_overlaps_normal_random[[i]][j] <- cancer2_overlaps
  }
  print(i)
}

save(all_top_overlaps_normal_random, file="Objects/all_top_overlaps_normal_random.Rdata")
load("Objects/all_top_overlaps_normal_random.Rdata")

Plot_Overlap_Heatmapx_Normal_random <- function(cancer1) {
  temp <- age_enrichment[[cancer1]]
  temp <- temp[temp$cluster != "grey",]
  age_df <- data.frame(age = temp$Module_age[temp$tissue_type == 'normal'])
  HA <- HeatmapAnnotation(df=age_df, which = "row",
                          col = list(age=c("UC" = "#F8766D", "Mixed" = "#00BA38", "MC" = "#619CFF")))
  col_fun = colorRamp2(c(0, 0.25, 0.75, 1), c("blue", "white", "yellow", "red"))
  print(Heatmap(as.matrix(all_top_overlaps_normal_random[[cancer1]]), 
                rect_gp = gpar(col = "white", lwd = 1),
                right_annotation = HA,
                row_split = age_df,
                column_title = cancer1,
                col = col_fun))
}

pdf("Overlaps_Age_Modules_normal_random.pdf", width = 10, height = 10)
for (cancer1 in names(cluster_assignments$normal)) { ##check normal object
  Plot_Overlap_Heatmapx_Normal_random(cancer1)
}
dev.off()



###Extract the level of overlap - tumour
mean_overlap_both <- vector()
for(tumour in names(cluster_assignments$tumour)){
  mean_overlap <- apply(all_top_overlaps[[tumour]], 1, mean)
  temp <- age_enrichment[[tumour]]
  temp <- temp[temp$tissue_type == "tumour",]
  temp <- temp[match(names(mean_overlap), temp$cluster),]
  mean_overlap_both <- rbind(mean_overlap_both,
                             cbind(tumour, "tumour", names(mean_overlap), temp$Module_age, unname(mean_overlap)))
}

for(tumour in names(cluster_assignments$tumour)){
  mean_overlap <- apply(all_top_overlaps_random[[tumour]], 1, mean)
  temp <- age_enrichment[[tumour]]
  temp <- temp[temp$tissue_type == "tumour",]
  temp <- temp[match(names(mean_overlap), temp$cluster),]
  mean_overlap_both <- rbind(mean_overlap_both,
                             cbind(tumour, "tumour_random", names(mean_overlap), temp$Module_age, unname(mean_overlap)))
}

#In normal
for(tumour in names(cluster_assignments$normal)){
  mean_overlap <- apply(all_top_overlaps_normal[[tumour]], 1, mean)
  temp <- age_enrichment[[tumour]]
  temp <- temp[temp$tissue_type == "normal",]
  temp <- temp[match(names(mean_overlap), temp$cluster),]
  mean_overlap_both <- rbind(mean_overlap_both,
                             cbind(tumour, "normal", names(mean_overlap), temp$Module_age, unname(mean_overlap)))
}

for(tumour in names(cluster_assignments$normal)){
  mean_overlap <- apply(all_top_overlaps_normal_random[[tumour]], 1, mean)
  temp <- age_enrichment[[tumour]]
  temp <- temp[temp$tissue_type == "normal",]
  temp <- temp[match(names(mean_overlap), temp$cluster),]
  mean_overlap_both <- rbind(mean_overlap_both,
                             cbind(tumour, "normal_random", names(mean_overlap), temp$Module_age, unname(mean_overlap)))
}
mean_overlap_both <- data.frame(mean_overlap_both)
colnames(mean_overlap_both) <- c("Tumour", "Tissue_type", "Module", "Age", "Mean_overlap")
mean_overlap_both$Mean_overlap <- as.numeric(as.character(mean_overlap_both$Mean_overlap))

mean_overlap_both$Age <- factor(mean_overlap_both$Age, levels=c("UC", "Mixed", "MC"))

library(ggplot2)
ggplot(mean_overlap_both, aes(x=Mean_overlap))+
  geom_density(aes(colour=Tissue_type))+
  theme_bw()

ks.test(x = mean_overlap_both$Mean_overlap[mean_overlap_both$Tissue_type == "normal"],
        y = mean_overlap_both$Mean_overlap[mean_overlap_both$Tissue_type == "tumour"],
        alternative = c("two.sided"))

ggplot(mean_overlap_both[mean_overlap_both$Tissue_type == "tumour",], aes(x=Mean_overlap))+
  geom_density(aes(colour=Age))+
  theme_bw()

ggplot(mean_overlap_both[mean_overlap_both$Tissue_type == "normal",], aes(x=Mean_overlap))+
  geom_density(aes(colour=Age))+
  theme_bw()

pdf("Density_Overlaps_Age_Modules_Tumour.pdf", height=4, width=6)
g1 <- ggplot(mean_overlap_both[mean_overlap_both$Tissue_type %in% c("tumour", "tumour_random"),], aes(x=Mean_overlap))+
  geom_density(aes(colour=Age, linetype=Tissue_type))+
  theme_bw()
print(g1)
dev.off()


pdf("Density_Overlaps_Age_Modules_Normal.pdf", height=4, width=6)
g1 <- ggplot(mean_overlap_both[mean_overlap_both$Tissue_type %in% c("normal", "normal_random"),], aes(x=Mean_overlap))+
  geom_density(aes(colour=Age, linetype=Tissue_type))+
  theme_bw()
print(g1)
dev.off()
