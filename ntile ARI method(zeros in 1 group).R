library(readxl)
melanoma_ori <- read.delim("~/Downloads/Research/data/melanoma.txt", header = TRUE, stringsAsFactors = FALSE)
#MelanomaDatasetWithClusters <- read.csv("~/Downloads/data/MelanomaDatasetWithClusters", stringsAsFactors = FALSE)
cell_type_ori <- read_excel("~/Downloads/Research/data/aad0501_Table_S3.xlsx")
excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")
GSE72056_melanoma_single_cell_revised_v2 <- read.delim("~/Downloads/Research/data/GSE72056_melanoma_single_cell_revised_v2.txt")
excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")

melanoma <- melanoma_ori
genes <- melanoma[,1]
melanoma <- melanoma[,-1]
melanoma <- as.matrix(melanoma)
rownames(melanoma) <- genes

col_index <- as.numeric(GSE72056_melanoma_single_cell_revised_v2[3,-1])
col_index[as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 2] <- 7

melanoma <- melanoma[, !(col_index == 0)]
col_index <- col_index[col_index != 0]

# #filter out 95% 0's
# ninety_five <- ncol(melanoma) * 0.95
# index <- numeric(0)
# for(i in 1:nrow(melanoma)){
#   if(sum(melanoma[i,] == 0) > ninety_five){
#     index <- c(index, i)
#   }
# }
# melanoma <- melanoma[-index,]

#直接把quantile改参数
#把0分出来
#每个gene做6次，选出最大的ARI
#利用cell type information
library(dplyr)
library(aricode)
set.seed(10000)
ARIs <- numeric(0)
max_indexes <- numeric(0)
#matrix(ncol = 2:7, nrow= # of genes)
m <- matrix(ncol = 6, nrow = nrow(melanoma))
colnames(m) <- as.character(1:6)
for(i in 1:nrow(melanoma)){
  expressions <- melanoma[i,]
  for(j in 1:6){
    zero_ind <- (expressions == 0)
    non_zero_ind <- !zero_ind
    clusters <- numeric(length(expressions))
    clusters[non_zero_ind] <- ntile(expressions[non_zero_ind], j);
    clusters[zero_ind] <- 0
    m[i,j] <- ARI(clusters, col_index)
  }
  max_indexes <- c(max_indexes, which(m[i,] == max(m[i,])))
  ARIs <- c(ARIs, max(m[i,]))
  print(i)
}
#orders <- order(m[,1])
orders <- order(ARIs)
#saveRDS(orders, file = "orders.rds")
#orders <- readRDS("orders.rds")
genes <- rownames(melanoma)
#marker <- genes[orders[(nrow(melanoma)-901):nrow(melanoma)]]
marker <- genes[orders[(nrow(melanoma)-101):nrow(melanoma)]]

library(Seurat)
data=melanoma
gene=row.names(data)
gene=unique(gene)
index=match(gene,row.names(data))
data=data[index,]
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,
                           project = "project")
pbmc <- NormalizeData(object = pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- RunPCA(object = pbmc, features = marker)
pbmc <- FindNeighbors(object = pbmc,dims = 1:8)
pbmc <- FindClusters(object = pbmc, reduction= "pca", resolution = 0.6)

# pbmc <- RunTSNE(object = pbmc, dims.use = 1:15, do.fast = TRUE)
# TSNE = Embeddings(pbmc[['tsne']])
# 
# plot(TSNE, col = c('red','pink','orange','blue',"yellow","white", "brown")[col_index], cex = 0.5)
# legend("topleft", c("T", "B", "Macro", "Endo", "CAF", "NK", "Tumor"), 
#        col = c('red','pink','orange','blue',"yellow","white", "brown"),
#        lty = 1, cex = 0.35 )
library(aricode)
ARI(pbmc$seurat_clusters, col_index)


