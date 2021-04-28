library(readxl)
melanoma_ori <- read.delim("~/Downloads/Research/data/melanoma.txt", header = TRUE, stringsAsFactors = FALSE)
#MelanomaDatasetWithClusters <- read.csv("~/Downloads/data/MelanomaDatasetWithClusters", stringsAsFactors = FALSE)
cell_type_ori <- read_excel("~/Downloads/Research/data/aad0501_Table_S3.xlsx")
excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")
GSE72056_melanoma_single_cell_revised_v2 <- read.delim("~/Downloads/Research/data/GSE72056_melanoma_single_cell_revised_v2.txt")
excludeGene_ori <- read.csv("~/Downloads/Research/data/excludeGene.txt", sep="")

# malignant <-melanoma[,as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 2]
# non_malignant <-melanoma[,as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 1]
# 
# col_index <- as.numeric(GSE72056_melanoma_single_cell_revised_v2[3,-1])[as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 1]

#library(Rtsne)

excludeGene <- excludeGene_ori

melanoma <- melanoma_ori
genes <- melanoma[,1]
melanoma <- melanoma[,-1]
melanoma <- as.matrix(melanoma)
rownames(melanoma) <- genes

res=ModalFilter(data=melanoma,geneK=20,cellK=20,width=2)
res=GeneFilter(obj=res)
ress=getMarker(obj=res,k=50,n=80)
length(ress$marker)

marker=ress$marker

####
#sum(marker_genes %in% marker) #171
####

# t_melanoma <- t(melanoma)
# reduced_melanoma <- t_melanoma[,colnames(t_melanoma) %in% marker]
# 
# ts <- Rtsne(reduced_melanoma)
# plot(ts$Y[,1], ts$Y[,2], col = as.factor(col_index), cex = 0.5)

col_index <- as.numeric(GSE72056_melanoma_single_cell_revised_v2[3,-1])
col_index[as.numeric(GSE72056_melanoma_single_cell_revised_v2[2,-1]) == 2] <- 7

library(Seurat)
library(dbscan)

obj <- res

data=obj$rawdata
gene=row.names(data)
gene=unique(gene)
index=match(gene,row.names(data))
data=data[index,]
marker=obj$marker
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,
                           project = "project")
pbmc <- NormalizeData(object = pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- RunPCA(object = pbmc, features = marker)
pbmc <- FindNeighbors(object = pbmc,dims = 1:8)
pbmc <- FindClusters(object = pbmc, reduction= "pca", resolution = 0.6)

train_col_index <- col_index[col_index != 0]
seurat_clusters <- pbmc$seurat_clusters[col_index!=0]
library(aricode)
ARI(seurat_clusters, train_col_index)  #0.37249


