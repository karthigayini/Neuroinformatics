**Saunder's data**

Dropseq data from paper - Saunders A*, Macosko E.Z*, Wysoker A, Goldman M, Krienen, F, de Rivera H, Bien E, Baum M, Wang S, Bortolin L, Goeva A, Nemesh J, Kamitaki N, Brumbaugh S, Kulp D and McCarroll, S.A. 2018. Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain. 2018. Cell. 174(4) P1015-1030.E16 ([DOI](https://doi.org/10.1016/j.cell.2018.07.028))

In this code Saunder's data in the count table format is downloaded and re-analysed in the way the house-data is analysed. By re-analyzing, the reference data - Saunder's can be compared with the query- in-houe data for annotating.  



- download the dropsequtil.tar into the location of the dge, data from http://dropviz.org/
- wget https://storage.googleapis.com/dropviz-downloads/static/DropSeq.util_2.0.tar.gz

```R
###################################################
###################################################
#### Written by Karthigayini Sivaprakasam
#### 2019-12-30
###################################################
###################################################

module load seurat/3.0.0

install.packages("DropSeq.util_2.0.tar.gz",repos=NULL)
library(Seurat)
library(DropSeq.util)
library(ggplot2)
library(dplyr)

################################################
################################################
#### Input

#read the different Saunder's data tables. This one is specific for STRIATUM data. But can be applied to other regions. 
saunder_dge=loadSparseDge("F_GRCm38.81.P60Striatum.raw.dge.txt.gz")
saunder_cluster=readRDS("F_GRCm38.81.P60Striatum.cluster.assign.RDS")
saunder_subclus=readRDS("F_GRCm38.81.P60Striatum.subcluster.assign.RDS")
saunder_annot=data.frame(readRDS("annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))
##subset the required cell type - here it is striatum 
saunder_annot=subset(saunder_annot,tissue=="STR")
saunder_outcome=data.frame(readRDS("F_GRCm38.81.P60Striatum.cell_cluster_outcomes.RDS"))
saunder_outcome$cluster=as.numeric(levels(saunder_outcome$cluster))[saunder_outcome$cluster]

#add the cluster  number and sub-cluster nos. 
saunder_annot$cluster=as.numeric(sapply(strsplit(saunder_annot$subcluster,'-'),`[`,1))

#Filter the cells that does not have any annotation as per the reference data table. 
saunder_outcome=saunder_outcome[saunder_outcome$cluster %in% saunder_annot$cluster,]
saunder_dge=saunder_dge[,colnames(saunder_dge) %in% rownames(saunder_outcome)]

#from 25645 90628 to  25645 77454

################################################
################################################
############SEURAT WORKFLOW##################

fileName="saunder"
saunder_obj=CreateSeuratObject(counts=saunder_dge)
saunder_obj[["percent.mt"]] <- PercentageFeatureSet(saunder_obj, pattern = "^mt")

apply(saunder_obj@meta.data[,2:4],2,range)
#     nCount_RNA nFeature_RNA percent.mt
#[1,]        459          401  0.3623188
#[2,]      40864         7051 31.8953324

###option 1 -linear transformation
### method is chosen based on the transformation of the input(house) data. 

saunder_obj=NormalizeData(saunder_obj,normalization.method = "LogNormalize") 
#identify highest variable feature
saunder_obj=FindVariableFeatures(saunder_obj, selection.method = "vst", nfeatures = 2000)

all.genes =rownames(saunder_obj)
all.genes=head(VariableFeatures(saunder_obj), 2000)
#load the data and run the RScript until the clusters
saunder_obj=ScaleData(saunder_obj, features = all.genes,model.use = "linear",vars.to.regress = c("percent.mito","batch","nCount_RNA"))

###option 2:sct method 
saunder_obj=SCTransform(saunder_obj,verbose=T,conserve.memory=T)
save.image("saunders_midanalysis.RData")

#######################################
#### for both options : pca visualization
saunder_obj=RunPCA(saunder_obj, features = VariableFeatures(object = saunder_obj))

#determine the dimensionality of the dataset use jackstraw 
saunder_obj= JackStraw(saunder_obj, num.replicate = 100,dims = 50)
saunder_obj= ScoreJackStraw(saunder_obj, dims = 1:50)
setEPS()
postscript(paste0(fileName,"_jackstraw.eps"))
JackStrawPlot(saunder_obj, dims = 1:50)
dev.off()

#set dims based on pvals
pvals=sort(which(saunder_obj@reductions$pca@jackstraw$overall.p.values[,2]>0.05))
dimension=ifelse(length(pvals)==0,50,pvals[1])

####cluster cells
saunder_obj<- FindNeighbors(saunder_obj, dims = 1:dimension)
saunder_obj <- FindClusters(saunder_obj, resolution = 0.5)

####UMAP and replace the cluster the no with the cluster no from the annotation 
saunder_obj <- RunUMAP(saunder_obj, dims = 1:dimension)

###replace the cluster nos in the object with the ones in the annotation
clusters=saunder_outcome$cluster
names(clusters)=rownames(saunder_outcome)
Idents(saunder_obj,var=active.ident)=clusters
saunder_obj$seurat_clusters=clusters

setEPS()
postscript(paste(fileName,"_umap.eps",sep=""))
#p1=DimPlot(saunder_obj, reduction = "umap", group.by = "orig.ident")
p2=DimPlot(saunder_obj, reduction = "umap", label = TRUE)
p2
dev.off()
save.image("allsamples.RData")

##### differential expression
cluster.markers =FindAllMarkers(saunder_obj,only.pos = T,logfc.threshold = 0.2) 
head(subset(cluster.markers,avg_logFC>0))
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(cluster.markers,paste0(fileName,"_clusterMarkers.csv"),sep="\t",row.names = F,quote = F)
save.image("saunder_str_cluster.RData")

################################################
################################################
##########overlap with the user data############ 

#load the saunder_str_cluster.RData here, having the saunder clusters, annotation table
###IMP USE R3.5 AND MODULE SEURAT 3


#read the user data
library(dplyr)
house_data=data.frame(read.table("/from/directory/clusterMarkers.csv",sep="\t", header=T))

#df 1
house_geneClus=subset(house_data,pct.1>=0.5& p_val_adj<=0.005,select=c("cluster","gene"))

#df 2
house_geneNo=as.data.frame(table(house_geneClus$cluster))
colnames(house_geneNo)=c("cluster","noOfGenes")

#df 3
saun_data=subset(cluster.markers,pct.1>=0.5 & p_val_adj <=0.005,select=c("cluster","gene"))
dist_saun_ann=data.frame(read.csv2("AK_cluster.csv",header=T,sep="\t"))

#df 4
ann_saunData=merge(dist_saun_ann,saun_data,by="cluster",all=T)
saun_data$cluster=as.numeric(levels(saun_data$cluster))[saun_data$cluster]

#check the no. of genes in each cluster before and afer the annotation
table(ann_saunData$cluster)
table(saun_data$cluster)

#the 3 frames we need house DE data, house_geneNo and the ann_saunData data frame. 

#convert the two gene pval DE tables into list 
annot_list=split(ann_saunData$gene,ann_saunData$class)
house_geneList=split(house_geneClus$gene,house_geneClus$cluster)

###overlap
#intersect the two gene list to get a list of list where every house cluster is intersected with the annotated cluster 

#to get the no. of overlaping genes 
intersection=lapply(house_geneList,function(x) lapply(annot_list,function(y) intersect(x,y)))
housedata_ncol=length(intersection)

#lengths of intersection , nrow is the no. of cluster in the annotation and ncol is the no. of clusters in the house data
overlaps=matrix(as.data.frame(rapply(intersection,length,how="list")),nrow=12,ncol=housedata_ncol)

rownames(overlaps)=names(intersection[[2]]) #annotation names
colnames(overlaps)=paste0("C_",(seq(0:(housedata_ncol-1))-1)) #house cluster no. 

###########################
###fishers exact test
##########################

#the goal is the compute the similarity of every house cluster with the annotated cluster
#dfs needed are 1. overlaps (list), annot_list

len_anno_list=lengths(annot_list)
names(len_anno_list)=NULL

#create a matrix for every cluster and generate a data frame with rows as the elements for fishers test
contingency_tbl=list()
for(i in 1:ncol(overlaps)){
	temp=matrix(overlaps[,i])
	temp=cbind(temp,rep(house_geneNo[i,2],12),len_anno_list,rep(15585,12))
	colnames(temp)=c("overlap","freq_clus","freq_ref","bs")
	temp=apply(temp,2,function(x) unlist(x))
	temp[,2:4]=sweep(temp[,2:4],1,temp[,1])
	contingency_tbl=append(contingency_tbl,list(temp),0)
}

#once you have this - do fishers for every row 
fishers_tbl=list()
for(i in 1:length(contingency_tbl)){
	df=as.data.frame(contingency_tbl[i])
	fisher=data.frame()
	for(j in 1:nrow(df)){
		df1=matrix(unlist(df[j,]),nrow=2)
		fisher[j,1]=fisher.test(df1,conf.level=0.95)$p.value
		fisher[j,2]=fisher.test(df1,conf.level=0.95)$estimate
	}
	colnames(fisher)=c("pval","or")
	rownames(fisher)=names(intersection[[2]])
	fishers_tbl=append(fishers_tbl,list(fisher),0)
}



#create a data frame of p values to use pheatmap
library(gplots)

pval_tbl=as.data.frame(matrix(,nrow=12,ncol=length(fishers_tbl)))
colnames(pval_tbl)=paste0("C_",(seq(0:(housedata_ncol-1))-1)) #house cluster no.
for(i in 1:length(fishers_tbl)){
	pval_tbl[,i]=as.vector(fishers_tbl[[i]][[1]])
}
pval_tbl[pval_tbl>0.05]=NA
pval_tbl$annot=rownames(fishers_tbl[[1]])

pval_ncol=ncol(pval_tbl)
pdf("heatmap_pval.pdf")
heatmap.2(as.matrix(pval_tbl[,-pval_ncol]),colsep=1:ncol(pval_tbl)-1,rowsep=1:nrow(pval_tbl),sepcolor="white",
cellnote=round(as.matrix(pval_tbl[,-pval_ncol]),2),notecex=0.3, Rowv=NA,Colv=NA,labRow=pval_tbl[,pval_ncol],cexRow=0.5,
margins=c(6,7),labCol=colnames(pval_tbl[1:pval_ncol]),key=TRUE,trace='none',na.color="lightgrey",lhei=c(2,8),lwid=c(2,8),
key.title=NA,sepwidth=c(0.03, 0.03),notecol="black")
dev.off()


#OR table
pval_tbl=as.data.frame(matrix(,nrow=12,ncol=length(fishers_tbl)))
colnames(pval_tbl)=paste0("C_",(seq(0:(housedata_ncol-1))-1)) #house cluster no.
for(i in 1:length(fishers_tbl)){
	pval_tbl[,i]=as.vector(fishers_tbl[[i]][[2]])
}
pval_tbl[pval_tbl==0]=NA
pval_tbl$annot=rownames(fishers_tbl[[1]])

pval_ncol=ncol(pval_tbl)
pdf("heatmap_or.pdf")
heatmap.2(as.matrix(pval_tbl[,-pval_ncol]),colsep=1:ncol(pval_tbl)-1,rowsep=1:nrow(pval_tbl),sepcolor="white",
cellnote=round(as.matrix(pval_tbl[,-pval_ncol]),2),notecex=0.3, Rowv=NA,Colv=NA,labRow=pval_tbl[,pval_ncol],cexRow=0.5,
margins=c(6,7),labCol=colnames(pval_tbl[1:pval_ncol]),key=TRUE,trace='none',na.color="lightgrey",lhei=c(2,8),lwid=c(2,8),
key.title=NA,sepwidth=c(0.03, 0.03),notecol="black")
dev.off()
################################################
################################################
```



