setwd("single_cell_analysis")
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(scRNAseq)
library(celldex)

# list of files 
untar("GSE176078_RAW.tar", exdir = "single_cell_analysis/Breast/data") ##26 samples
file_list1 <- untar("GSE176078_RAW.tar", list = TRUE, exdir = "single_cell_analysis/Breast/data")
file_list1
obj_df1 <- data.frame(file = file_list1)
obj_df1$obj <- substr(obj_df1$file,1,10)
obj_df1

lst <- scRNAseq::listDatasets()
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
hESCS <- LaMannoBrainData("human-es")
Add_cell.type <- function(seurat_obj) {
  pred.df <- SingleR(test = as.SingleCellExperiment(seurat_obj), 
                     ref = hpca.se, assay.type.test = 1, 
                     labels = hpca.se$label.main)
  seurat_obj$cell.type <- pred.df$labels
  return(seurat_obj)
}

# object: QC, Normalisation, Dim Reduction, Annotation
setup_obj <- function(x, file_name, mt_cutoff) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < mt_cutoff)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- row.names(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- JackStraw(x, num.replicate = 100)
  x <- ScoreJackStraw(x, dims = 1:20)
  x <- FindNeighbors(x, dims = 1:20)
  x <- FindClusters(x, resolution = 1)
  x <- RunUMAP(x, dims = 1:20)
  x <- Add_cell.type(x)
  saveRDS(x, file = file_name)
  return(x)
}

# unzip all files
for (i in c(1:nrow(obj_df1))) {
  path_dir <- "single_cell_analysis/Breast/data"
  file_dir = paste0(path_dir,"/", obj_df1$file[i])
  print(file_dir)
  untar(file_dir, exdir = "single_cell_analysis/Breast/data/data")
}
obj_df1$file_alt <- list.files("single_cell_analysis/Breast/data/data")
obj_df1

library(ggplot2)
pdf("single_cell_analysis/Breast/Breast_Data_Quality.pdf")
for (i in c(1:nrow(obj_df1))) {
  path_dir <- "single_cell_analysis/Breast/data/data"
  file_dir = paste0(path_dir,"/", obj_df1$file_alt[i])
  print(file_dir)
  files <- list.files(file_dir)
  data <- ReadMtx(mtx = paste0(file_dir, "/", files[3]),
                  cells = paste0(file_dir, "/", files[1]), 
                  features = paste0(file_dir, "/", files[2]), feature.column = 1)
  obj <- CreateSeuratObject(counts = data, project = "SeuratProject", 
                            assay = "RNA", min.cells = 3, min.features = 200)
  obj_name <- paste0(obj_df1$obj[i], ".Rds")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  g <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+
    ggtitle(obj_name)
  print(g)
}
dev.off()

datasets_breast <- list()
for (i in c(1:nrow(obj_df1))) {
  path_dir <- "single_cell_analysis/Breast/data/data"
  file_dir = paste0(path_dir,"/", obj_df1$file_alt[i])
  print(file_dir)
  files <- list.files(file_dir)
  data <- ReadMtx(mtx = paste0(file_dir, "/", files[3]),
                  cells = paste0(file_dir, "/", files[1]), 
                  features = paste0(file_dir, "/", files[2]), feature.column = 1)
  obj <- CreateSeuratObject(counts = data, project = "SeuratProject", 
                            assay = "RNA", min.cells = 3, min.features = 200)
  obj_name <- paste0(obj_df1$obj[i], ".Rds")
  datasets_breast[[i]] <- setup_obj(obj, obj_name, 15)
}
names(datasets_breast) <- obj_df1$file_alt
save(datasets_breast, file="single_cell_analysis/Breast/datasets_breast.Rdata")

load("single_cell_analysis/Breast/datasets_breast.Rdata")
load("Objects/age_enrichment.Rdata")
load("Objects/cluster_assignments.Rdata")

pdf("single_cell_analysis/Breast/Breast_Per_Mito.pdf")
for (i in c(1:length(datasets_breast))) {
  print(obj_df1$obj[i])
  data <- datasets_breast[[i]]
  obj_name <- paste0(obj_df1$obj[i], ".Rds")
  g <- FeaturePlot(datasets_breast[[i]], "percent.mt") + ggtitle(obj_name)
  print(g)
}
dev.off()


# add module scores
module_list <- unlist(cluster_assignments,recursive = FALSE, use.names = TRUE)
class(module_list)
module_list_breast <- c(module_list$tumour.BRCA)
module_type <- c(rep("tumour_", length(module_list$tumour.BRCA)))
names(module_list_breast) <- paste0(module_type, "_",names(module_list_breast))

## for reference
module_list_df <- data.frame(Geneset = 1:length(module_list_breast), 
                             Module = paste0(names(module_list_breast)))
module_list_df



for (i in c(1:length(datasets_breast))) {
  print("--------------Adding Module Score-------------")
  print(obj_df1$obj[i])
  data <- datasets_breast[[i]]
  Ctrl <- length(data)
  data <- AddModuleScore(data, features = module_list$tumour.BRCA, ctrl = Ctrl, 
                         name = "Geneset_", search = TRUE)
  file_name <- paste0(obj_df1$obj[i], "_M.Rds")
  saveRDS(data@meta.data, file_name)
  print(obj_df1$obj[i])
}


library(pheatmap)
library(patchwork)
plot_ModuleScore <- function(mtx, obj_name, srt_obj) {
  info = obj_name
  unsorted_mtx <- mtx[,c(7:ncol(mtx))]
  sorted_mtx <- unsorted_mtx[order(unsorted_mtx[,"cell.type"]),]
  cell_label <- data.frame(cell.type = sorted_mtx[,1])
  rownames(cell_label) <- rownames(sorted_mtx)
  score_mtx <- sorted_mtx[,c(2:ncol(sorted_mtx))]
  for (cluster in names(cluster_assignments$tumour$BRCA)) {
    print(VlnPlot(srt_obj, features = cluster, group.by = "cell.type")+plot_annotation(info))
    print(FeaturePlot(srt_obj, features = cluster))
  }
}


for (i in c(1:nrow(obj_df1))) {
  obj_name <- obj_df1$obj[i]
  print(obj_name)
  obj <- datasets_breast[[i]]
  print(DimPlot(obj, group.by = "cell.type"))
  print("Dimplot done!")
  m_obj <- readRDS(paste0("single_cell_analysis/", obj_name,"_M.Rds"))
  colnames(m_obj)[8:55] <- names(cluster_assignments$tumour$BRCA)
  obj@meta.data <- m_obj
  print("Obj ready!")
  pdf(paste0("single_cell_analysis/Breast/Breast_All_Plots_", obj_name, ".pdf"))
  plot_ModuleScore(m_obj, obj_name, obj)
  dev.off()
}


FeaturePlot(GSM5354515, features = "EPCAM")

# wilcoxon between tumour and normal
library(ggplot2)
pdf("single_cell_analysis/Breast/Breast_Tumour_Cells.pdf")
for (i in c(1:nrow(obj_df1))) {
  if (i != 5){ ##5 has no EPCAM or MKI67
    obj_name <- obj_df1$obj[i]
    print(obj_name)
    obj <- datasets_breast[[i]]
    m_obj <- readRDS(paste0("single_cell_analysis/", obj_name,"_M.Rds"))
    colnames(m_obj)[8:55] <- names(cluster_assignments$tumour$BRCA)
    obj@meta.data <- m_obj
    print("Obj ready!")
    print(table(obj$cell.type))
    print(DimPlot(obj, label = TRUE, label.box = TRUE)+ggtitle(obj_name))
    print(FeaturePlot(obj, features = c('EPCAM'))) 
  }
}
dev.off()

##Clusters of tumour cells based on the visualization
obj_df1$tumour_cluster <- list(
  c(14,15,16,9),
  NA,
  c(10),
  c(2,7),
  NA,
  c(4),
  c(NA),
  c(NA),
  c(10,12,14,15),
  c(9,1,3,2,5,6,14),
  c(0,2,1,6,3,13),
  NA,
  c(11,12),
  c(0),
  c(4,2,0),
  c(5),
  c(12,9,5,10, 16),
  c(5,9,15),
  c(8,12,21),
  c(10,13,5,0,11),
  NA,
  c(1,2,5),
  c(1,8),
  c(8,2,0,3,5,10),
  c(1,6,13,12),
  c(1,3,13,6)
)

obj_df2 <- obj_df1[is.na(obj_df1$tumour_cluster) == FALSE,]

##Determine number of "tumour" cells
n_tumour_cells <- vector()
T_cells <- list()
N_cells <- list()
pdf("single_cell_analysis/Breast/Breast_Tumour_Cells_Assignment.pdf")
for (i in c(1:nrow(obj_df1))) {
  obj_name <- obj_df1$obj[i]
  if(obj_name %in% obj_df2$obj){
    print(obj_name)
    obj <- datasets_breast[[i]]
    obj@meta.data$is.tumour <- "normal"
    obj@meta.data$is.tumour[obj@meta.data$seurat_clusters %in% obj_df1$tumour_cluster[[i]]] <- "tumour"
    
    T_cells[[obj_name]] <- rownames(obj@meta.data[obj@meta.data$is.tumour == "tumour",])
    N_cells[[obj_name]] <- rownames(obj@meta.data[obj@meta.data$is.tumour == "normal",])
    
    print(DimPlot(obj, group.by="is.tumour")+ggtitle(obj_name))
    print(VlnPlot(obj, features = c('EPCAM'), group.by="is.tumour")) 
    n_tumour_cells <- rbind(n_tumour_cells, c(obj_name, table(obj@meta.data$is.tumour)))
  }
}
dev.off()

n_tumour_cells <- as.data.frame(n_tumour_cells)
n_tumour_cells$normal <- as.numeric(as.character(n_tumour_cells$normal))
n_tumour_cells$tumour <- as.numeric(as.character(n_tumour_cells$tumour))

hist(n_tumour_cells$tumour, breaks=100)

sample_to_keep <- n_tumour_cells[n_tumour_cells$tumour >= 500,]



wilcox_df <- data.frame(row.names = names(cluster_assignments$tumour$BRCA))

for (i in c(1:nrow(obj_df1))) {
  obj_name <- obj_df1$obj[i]
  if(obj_name %in% sample_to_keep$V1){
    print(obj_name)
    m_obj <- readRDS(paste0("single_cell_analysis/", obj_name,"_M.Rds"))
    colnames(m_obj)[8:55] <- names(cluster_assignments$tumour$BRCA)
    m_obj_T <- m_obj[rownames(m_obj) %in% T_cells[[obj_name]],8:ncol(m_obj)]
    m_obj_N <- m_obj[rownames(m_obj) %in% N_cells[[obj_name]],8:ncol(m_obj)]
    mean_score_T <- apply(m_obj_T, 2, mean)
    mean_score_N <- apply(m_obj_N, 2, mean)
    diff_scores <- mean_score_T-mean_score_N
    wilcox_results_greater <- sapply(names(cluster_assignments$tumour$BRCA), function(cluster){
      return(wilcox.test(m_obj_T[,cluster], m_obj_N[,cluster], alternative = "greater")$p.value)
    })
    wilcox_results_less <- sapply(names(cluster_assignments$tumour$BRCA), function(cluster){
      return(wilcox.test(m_obj_T[,cluster], m_obj_N[,cluster], alternative = "less")$p.value)
    })
    
    wilcox_results_two <- sapply(names(cluster_assignments$tumour$BRCA), function(cluster){
      return(wilcox.test(m_obj_T[,cluster], m_obj_N[,cluster])$p.value)
    })
    
    wilcox_results <- cbind(obj_name, names(wilcox_results_greater), unname(mean_score_T),
                            unname(mean_score_N), unname(diff_scores), 
                            unname(wilcox_results_greater), unname(wilcox_results_less), 
                            unname(wilcox_results_two))
    
    wilcox_df <- rbind(wilcox_df, wilcox_results)
  }
}

colnames(wilcox_df) <- c("Sample", "Module", "Mean_tumour_score", "Mean_normal_score", "Diff", "P.Value_greater", "P.Value_less", "P.Value_two")
wilcox_df$Diff <- as.numeric(as.character(wilcox_df$Diff))
wilcox_df$P.Value_greater <- as.numeric(as.character(wilcox_df$P.Value_greater))
wilcox_df$P.Value_less <- as.numeric(as.character(wilcox_df$P.Value_less))
wilcox_df$P.Value_two <- as.numeric(as.character(wilcox_df$P.Value_two))
wilcox_df <- wilcox_df[wilcox_df$Module != "grey",]
wilcox_df$Adj.P.Value_greater <- p.adjust(wilcox_df$P.Value_greater, method="BH")
wilcox_df$Adj.P.Value_less <- p.adjust(wilcox_df$P.Value_less, method="BH")
wilcox_df$Adj.P.Value_two <- p.adjust(wilcox_df$P.Value_two, method="BH")
wilcox_df$Sig_greater <- ifelse(wilcox_df$Adj.P.Value_greater < 0.01, "Y", "N")
wilcox_df$Sig_less <- ifelse(wilcox_df$Adj.P.Value_less < 0.01, "Y", "N")
wilcox_df$Sig_two <- ifelse(wilcox_df$Adj.P.Value_two < 0.01, "Y", "N")


load("Objects/cluster_assignments.Rdata")
load("Objects/age_enrichment.Rdata")

ages_tumour <- age_enrichment$BRCA[age_enrichment$BRCA$tissue_type == "tumour",]


wilcox_df$Age <- ages_tumour$Module_age[match(wilcox_df$Module, ages_tumour$cluster)]

temp1 <- cbind(wilcox_df[,c("Sample", "Module", "Mean_tumour_score", "Age")], Tissue="Tumour")
temp2 <- cbind(wilcox_df[,c("Sample", "Module", "Mean_normal_score", "Age")], Tissue="Normal")
colnames(temp1) <- c("Sample","Module","Mean_score","Age","Tissue")
colnames(temp2) <- c("Sample","Module","Mean_score","Age","Tissue")
wilcox_df2 <- rbind(temp1, temp2)
wilcox_df2$Mean_score <- as.numeric(as.character(wilcox_df2$Mean_score))

ggplot(wilcox_df2, aes(x=Age, y=Mean_score))+
  geom_boxplot(aes(fill=Tissue))

temp3 <- wilcox_df[,c("Sample", "Module", "Diff", "Age")]
temp3$Diff
temp3$Age <- factor(temp3$Age, levels=c("UC", "Mixed", "MC"))

pdf("single_cell_analysis/Breast/Breast_Diff_Scores.pdf", height=4, width=5)
g <- ggplot(temp3, aes(x=Age, y = Diff))+
  geom_boxplot(aes(fill=Age))+
  geom_hline(yintercept = 0, linetype = 2, colour="grey20")+
  geom_boxplot(aes(fill=Age))+
  ylab("Difference in expression score")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

##How many times each is significant

recurrence_greater <- table(wilcox_df$Module, wilcox_df$Sig_greater)
recurrence_greater_df <- data.frame(Module = rownames(recurrence_greater), 
                                    N = recurrence_greater[,1], Y = recurrence_greater[,2])
recurrence_less <- table(wilcox_df$Module, wilcox_df$Sig_less)
recurrence_less_df <- data.frame(Module = rownames(recurrence_less), 
                                    N = recurrence_less[,1], Y = recurrence_less[,2])
recurrence_two <- table(wilcox_df$Module, wilcox_df$Sig_two)
recurrence_two_df <- data.frame(Module = rownames(recurrence_two), 
                                 N = recurrence_two[,1], Y = recurrence_two[,2])

recurrence_df <- recurrence_greater_df
colnames(recurrence_df) <- c("Module", "N_greater", "Y_greater")
recurrence_df$N_less <- recurrence_less_df$N
recurrence_df$Y_less <- recurrence_less_df$Y
recurrence_df$N_two <- recurrence_two_df$N
recurrence_df$Y_two <- recurrence_two_df$Y

tumour_specific_modules <- recurrence_df[recurrence_df$Y_greater >= 6,] ##more than 50%

##33/47 modules are tumour-specific

save(tumour_specific_modules, file = "single_cell_analysis/Breast/tumour_specific_modules.Rdata")

tumour_specific_modules$Age <- ages_tumour$Module_age[match(tumour_specific_modules$Module, ages_tumour$cluster)]


