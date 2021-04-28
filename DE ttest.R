#helper function
vec_if_contains <- function(vec, pattern){
  ret <- logical(0)
  for(i in 1:length(vec)){
    if(str_contains(vec[i], pattern)){
      ret <- c(ret, TRUE)
    }else{
      ret <- c(ret, FALSE)
    }
  }
  ret
}


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

m <- matrix(ncol=21, nrow=nrow(melanoma))
for(i in 1:nrow(melanoma)){
  expressions <- melanoma[i,]
  ind <- 0
  for(j in 1:6){
    for(k in (j+1):7){
      ind <- ind + 1
      m[i,ind] <- t.test(expressions[col_index == j], expressions[col_index == k])$p.value
    }
  }
  print(i)
}

colnames(m) <- c("1:2","1:3","1:4","1:5","1:6","1:7","2:3","2:4","2:5","2:6","2:7","3:4","3:5","3:6","3:7",
                 "4:5","4:6","4:7","5:6","5:7","6:7")
rownames(m) <- rownames(melanoma)
scale_m <- (m < 0.05)

#remove NA
index <- numeric(0)
for(i in 1:nrow(scale_m)){
  if(any(is.na(scale_m[i,]))){
    index <- c(index, i)
  }
}
#from here
scale_m_no_na <- scale_m[-index,]

col_names <- colnames(scale_m_no_na)

genes <- character(0)
#all 7 types: found 3 genes
marker_7 <- numeric(0)
for(i in 1:nrow(scale_m_no_na)){
  if(all(scale_m_no_na[i,], na.rm=T)){
    marker_7 <- c(marker_7, i)
  }
}
genes <- c(genes, rownames(scale_m_no_na)[marker_7])

#only 6 types
#remove already-found genes
scale_m_no_na <- scale_m_no_na[-marker_7,]
marker_6 <- numeric(0)
for(i in 1:7){
  current_ind <- !vec_if_contains(col_names, as.character(i))
  scale_m_to_examine <- scale_m_no_na[, current_ind]
  for(q in 1:nrow(scale_m_to_examine)){
    if(all(scale_m_to_examine[q,])){
      marker_6 <- c(marker_6, q)
    }
  }
}
genes <- c(genes, rownames(scale_m_no_na)[marker_6])

#only 5 types
#remove already-found genes
scale_m_no_na <- scale_m_no_na[-marker_6,]
marker_5 <- numeric(0)
for(i in 1:6){
  for(j in (i+1):7){
    current_ind <- !(vec_if_contains(col_names, as.character(i)) | vec_if_contains(col_names, as.character(j)))
    scale_m_to_examine <- scale_m_no_na[, current_ind]
    for(q in 1:nrow(scale_m_to_examine)){
      if(all(scale_m_to_examine[q,])){
        marker_5 <- c(marker_5, q)
      }
    }
  }
}
genes <- c(genes, rownames(scale_m_no_na)[marker_5])

marker <- genes[1:902]

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

ARI(pbmc$seurat_clusters, col_index)













#only 4 types
#remove already-found genes
scale_m_no_na <- scale_m_no_na[-marker_5,]
marker_4 <- numeric(0)
for(i in 1:5){
  for(j in (i+1):6){
    for(k in (j+1):7){
      current_ind <- !(vec_if_contains(col_names, as.character(i)) | vec_if_contains(col_names, as.character(j)) | vec_if_contains(col_names, as.character(k)))
      scale_m_to_examine <- scale_m_no_na[, current_ind]
      for(q in 1:nrow(scale_m_to_examine)){
        if(all(scale_m_to_examine[q,])){
          marker_4 <- c(marker_4, q)
        }
      }
    }
  }
}

  
  
  
  

#(which(unlist(strsplit(colnames(m),split=":")) %in% 1) + 1) / 2
#floor((which(unlist(strsplit(colnames(m),split=":")) %in% 2) + 1) / 2)

result <- 0
for(i in c(3,2,1,3.5,6.1,1.4,5.3,6,2.4)){
  result <- result + log(i)
}
mu = result/9

result <- 0
for(i in c(3,2,1,3.5,6.1,1.4,5.3,6,2.4)){
  result <- result + (log(i)-1.058246)^2
}
sig = sqrt(result/9)
 
result <- 0
for(i in c(3,2,1,3.5,6.1,1.4,5.3,6,2.4)){
  result <- result + (log(i)-mu)*(-2)*sig^(-3)
} 
result

result <- 9*sig^(-2)
for(i in c(3,2,1,3.5,6.1,1.4,5.3,6,2.4)){
  result - 3*(log(i)-mu)^2*sig^(-4)
} 
result

A=matrix(c(-24.4,-9.103829e-15,-9.103829e-15,24.4),ncol=2)
eigen(A)

