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
# same as table(col_index)
each_type_num <- numeric(0)
for(i in 1:7){
  each_type_num <- c(each_type_num, sum(col_index == i))
}

n <- 60
sim_times <- 10000
vec <- numeric(0)
for(i in 1:6){
  vec <- c(vec, rep(i, n))
}
ARIs <- numeric(sim_times)
for(i in 1:sim_times){
  rand_vec <- vec[sample(length(vec))]
  ARIs[i] <- ARI(vec, rand_vec)
}
mean(ARIs)





