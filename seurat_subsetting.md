

Subsetting the cells based on expression of specific genes or for clusters of interest



```R
###################################################
###################################################
#### Written by Karthigayini Sivaprakasam
#### 2019-12-30
###################################################
###################################################

module load seurat/3.0.0

library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

setwd("SettheWorkingDirectory")
fileName="subsetting"

##########################
##########################
#### subsetting based on genes of interest 
load("noMito.RData")
cells_dspn=as.data.frame(dspn_obj[["SCT"]]@data[c(5904,5903,4706,4707),]) #normalized counts of cells from genes of interest
list_dspn_nonzerocells=apply(cells_dspn,1,function(x) names(cells_dspn[x!=0])) #get the cell names for gene
table(sapply(strsplit(list_dspn_nonzerocells$Foxp1,"_"),"[",1))
                             
drd1_cells=colnames(drds[,which(colSums(drds)!=0)])
#get the names of the cells that have the genes of interest expressed. 
drd_counts=merge_mat[,colnames(merge_mat) %in% drd1_cells]
dim(drd_counts)

# continue with seurat workflow
                           
#############################
#############################
#### subset for clusters
load("filtr_obj.RData")
dspn_cells=colnames(subset(filtr_obj,seurat_clusters==0 |seurat_clusters==1))
dspn_counts=merge_mat[,colnames(merge_mat) %in% dspn_cells]

#follow the full linear/sct method from here 
dspn_obj=CreateSeuratObject(dspn_counts)

```

