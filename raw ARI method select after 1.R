# library(readxl)
# melanoma_ori <- read.delim("~/Downloads/Research/data/melanoma.txt", header = TRUE, stringsAsFactors = FALSE)
# #MelanomaDatasetWithClusters <- read.csv("~/Downloads/data/MelanomaDatasetWithClusters", stringsAsFactors = FALSE)
# cell_type_ori <- read_excel("~/Downloads/Research/data/aad0501_Table_S3.xlsx")
# excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")
# GSE72056_melanoma_single_cell_revised_v2 <- read.delim("~/Downloads/Research/data/GSE72056_melanoma_single_cell_revised_v2.txt")
# excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")

melanoma <- melanoma_ori
genes <- melanoma[,1]
melanoma <- melanoma[,-1]
melanoma <- as.matrix(melanoma)
rownames(melanoma) <- genes

col_index <- as.numeric(GSE72056_melanoma_single_cell_revised_v2[3,-1])
col_index[as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 2] <- 7
# 
# melanoma <- melanoma[, !(col_index == 0)]
# col_index <- col_index[col_index != 0]
#train_col_index <- col_index[1:2100]
train_col_index <- col_index
test_col_index <- col_index[2101:length(col_index)]

#train <- melanoma[,1:2100]
train <- melanoma
test <- melanoma[,1:ncol(melanoma)]

library(dplyr)
library(aricode)
ARIs <- numeric(0)
for(i in 1:nrow(train)){
  expressions <- train[i,]
  clusters <- ntile(expressions, 7)
  ARIs <- c(ARIs, ARI(clusters, train_col_index))
}

orders <- order(ARIs)
genes <- rownames(melanoma)
marker <- genes[orders <= 902]

library(Seurat)
data=melanoma
gene=row.names(data)
gene=unique(gene)
index=match(gene,row.names(data))
data=data[index,]
#marker=obj$marker
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,
                           project = "project")
pbmc <- NormalizeData(object = pbmc, verbose = FALSE)
#pbmc <- FindVariableFeatures(object = pbmc, selection.method = "disp", nfeatures = 902)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- RunPCA(object = pbmc, features = marker)
pbmc <- FindNeighbors(object = pbmc,dims = 1:8)
pbmc <- FindClusters(object = pbmc, reduction= "pca", resolution = 0.6)
#pbmc <- FindClusters(object = pbmc, resolution = 0.6)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:15, do.fast = TRUE)
#TSNE=pbmc@dr$tsne@cell.embeddings
TSNE = Embeddings(pbmc[['tsne']])

train_col_index <- train_col_index[col_index != 0]
seurat_clusters <- pbmc$seurat_clusters[col_index != 0]

TSNE <- TSNE[col_index != 0,]
plot(TSNE, col = c('red','pink','orange','blue',"yellow","white", "brown")[train_col_index], cex = 0.5)
legend("topleft", c("T", "B", "Macro", "Endo", "CAF", "NK", "Tumor"), 
       col = c('red','pink','orange','blue',"yellow","white", "brown"),
       lty = 1, cex = 0.35 )

ARI(seurat_clusters, train_col_index) #0.3759673

# overlap between markers & known markers

## dbscan
dbscanCluster=dbscan(TSNE,eps=1.2,minPts=15)$cluster
ARI(dbscanCluster, train_col_index) # 0.116988


