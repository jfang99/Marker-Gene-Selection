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

library(dplyr)
library(aricode)
set.seed(10000)
m <- matrix(ncol = 6, nrow = nrow(melanoma))
colnames(m) <- as.character(2:7)
for(i in 1:nrow(melanoma)){
  expressions <- melanoma[i,]
  for(j in 1:6){
    clusters <- ntile(expressions, j+1)
    m[i,j] <- ARI(clusters, col_index)
  }
  print(i)
}

#scale m
scaled_m <- m
for(i in 1:6){
  scaled_m[,i] <- scaled_m[,i] / max(scaled_m[,i])
}
ARIs <- numeric(0)
max_indexes <- numeric(0)
for(i in 1:nrow(melanoma)){
  max_indexes <- c(max_indexes, which(scaled_m[i,] == max(scaled_m[i,])) + 1)
  ARIs <- c(ARIs, max(scaled_m[i,]))
}

orders <- order(ARIs)
genes <- rownames(melanoma)
marker <- genes[orders[(nrow(melanoma)-901):nrow(melanoma)]]

library(Seurat)
data=melanoma
gene=row.names(data)
gene=unique(gene)
index=match(gene,row.names(data))
data=data[index,]
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,
                           project = "project")
pbmc <- NormalizeData(object = pbmc, verbose = FALSE)
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

ARI(pbmc$seurat_clusters, col_index)


