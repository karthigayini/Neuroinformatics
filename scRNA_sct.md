Single cell RNA processing using SCTransform/linear regression with either manually aligned samples or count table from cell ranger.

vignette - https://satijalab.org/seurat/v3.0/sctransform_vignette.html

This script is modified to fit the NA-2019 project

```R
module load R/3.5.1
library(Seurat)
library(sctransform)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)

rm(list=ls())
setwd(#setTheDirectory#)
#############################################################
#set input variables
#input sex genes as a csv or txt file, tab separated with column-1 as chromosomeName and column-2 as genename
sex_genes=data.frame(read.table("/work/Neuroinformatics_Core/s186964/NA_scRna/texts/sex_genes_mm10Vm17.txt",sep="\t",header=F,stringsAsFactors = F))
colnames(sex_genes)=c("chrom","geneName")

#list the sample names
sampleList=c("NA1","NA2","NA3","NA4","NA5","NA6","NA7","NA8","NA9","NA10","NA11","NA12")
#############################################################
#############################################################
####merging the count matrices into one 

#create an empty matrix
merge_mat=data.frame(stringsAsFactors = F)

#option 1: for manually aligned samples
setwd("/work/Neuroinformatics_Core/s186964/NA_scRna/texts")
for(iSample in 1:length(sampleList)){
  filename=paste(sampleList[iSample], "_Counts_NovaSeq.txt", sep="")
  if (file.exists(filename) & file.info(filename)$size != 0) #if the file exists and not empty
    print(filename)
  count_mat=data.frame(read.table(filename,header = T, stringsAsFactors = F, sep = "\t"))
  count_mat=count_mat[setdiff(rownames(count_mat),sex_genes$geneName),] 
  colnames(count_mat)=paste(sampleList[iSample],colnames(count_mat),sep="_")
  count_mat$gene=rownames(count_mat)
  rownames(count_mat)=NULL
  if(iSample==1){
    merge_mat=count_mat
  }else{
    merge_mat=merge(merge_mat,count_mat,by="gene",all.x=T,all.y=T)
  }
}


#option 2:for cell ranger data
for(iSample in 1:length(sampleList)){
  print(iSample)
  count_mat=as.data.frame(Read10X(data.dir=file.path(sampleList[iSample],"filtered_feature_bc_matrix")))
  count_mat=count_mat[setdiff(rownames(count_mat),sex_genes$geneName),] 
  colnames(count_mat)=paste(sampleList[iSample],colnames(count_mat),sep="_")
  count_mat$gene=rownames(count_mat)
  rownames(count_mat)=NULL
  if(iSample==1){
    merge_mat=count_mat
  }else{
    merge_mat=merge(merge_mat,count_mat,by="gene",all.x=T,all.y=T)
  }
}


#for both the options
merge_mat[is.na(merge_mat)]=0
rownames(merge_mat)=merge_mat$gene
merge_mat$gene=NULL

#############################################################
#############################################################

####create seurat object and calculate % mito. 
first_complete_obj=CreateSeuratObject(counts = merge_mat)
first_complete_obj[["percent.mt"]] <- PercentageFeatureSet(first_complete_obj, pattern = "^mt") # "MT for humans"

first_complete_obj@meta.data %>% group_by(orig.ident) %>% summarize(mean(nFeature_RNA))
##Usually, the no. of Mito genes in MM10 is 13 but, cell ranger gives 37 genes

genenames=rownames(first_complete_obj)
no_mito_genes=genenames[!grepl("^mt",genenames)]
length(genenames)
length(no_mito_genes)

####create a named vector of %mito values with names as cell barcodes
keepMito_all=first_complete_obj@meta.data$percent.mt
names(keepMito_all)=rownames(first_complete_obj@meta.data)

######input meta data
scx_df=data.frame(matrix(NA,nrow=nrow(first_complete_obj@meta.data)))
scx_df$sex="N"
scx_df$genotype="N"
scx_df$batch="N"
scx_df$orig.ident=as.character(first_complete_obj$orig.ident)
str(scx_df)
scx_df[,1]=NULL
rownames(scx_df)=rownames(first_complete_obj@meta.data)

#based on sample replace the value
scx_df=mutate(scx_df,sex = case_when(
	orig.ident=="NA3" ~ "F",
	orig.ident=="NA4" ~ "F",
	orig.ident=="NA5" ~ "F",
	orig.ident=="NA7" ~ "F",
	TRUE ~ "M"
	)) 
scx_df=mutate(scx_df,batch=case_when(
	orig.ident=="NA1" ~ "1",
	orig.ident=="NA2" ~ "1",
	orig.ident=="NA3" ~ "1",
	orig.ident=="NA4" ~ "2",
	orig.ident=="NA5" ~ "2",
	orig.ident=="NA6" ~ "3",
	orig.ident=="NA7" ~ "3",
	orig.ident=="NA8" ~ "3",
	TRUE ~'4'
))

scx_df=mutate(scx_df,genotype=case_when(
	orig.ident=="NA3" ~ "control",
	orig.ident=="NA5" ~ "control",
	orig.ident=="NA7" ~ "control",
	orig.ident=="NA1" ~ "Foxp1",
	orig.ident=="NA6" ~ "Foxp1",
	orig.ident=="NA10" ~ "Foxp1",
	orig.ident=="NA8" ~ "Foxp2",
	orig.ident=="NA11" ~ "Foxp2",
	orig.ident=="NA12" ~ "Foxp2",
	TRUE ~'double_cko'
))

#create named vector of the meta data
sex=scx_df$sex
batch=scx_df$batch
genotype=scx_df$genotype
names(sex)=rownames(scx_df)
names(batch)=rownames(scx_df)
names(genotype)=rownames(scx_df)

#############################################################
#############################################################

####remove Mito genes and genes with no counts in all the cells. 
merge_mat=merge_mat[rowSums(merge_mat)!=0,]
merge_mat=merge_mat[rownames(merge_mat) %in% no_mito_genes,]

### create seurat object from the above data frame and add the meta data to it.
noMito_obj=CreateSeuratObject(merge_mat)
noMito_obj=AddMetaData(noMito_obj,keepMito_all,"percent.mito")
noMito_obj=AddMetaData(noMito_obj,sex,"sex")
noMito_obj=AddMetaData(noMito_obj,batch,"batch")
noMito_obj=AddMetaData(noMito_obj,genotype,"genotype")

#############################################################
#############################################################
####thresholding

apply(noMito_obj@meta.data[,c(2:4)],2,range,na.rm=T) #will give total range
table(noMito_obj$orig.ident) 												# total no. of cells (rows) for every sample
noMito_obj@meta.data %>% group_by(orig.ident) %>% summarise(min_umi=min(nCount_RNA),max_umi=max(nCount_RNA)) #this is will range of Umi count from every sample


#gene frequency can be plotted to fix the cut-off
#change the breaks and vline intercept based on the histogram
#change lower and upper count based on sample

percent.mit=1.5
upper_feature=6000


ggplot(noMito_obj@meta.data,aes(x=nFeature_RNA))+geom_histogram(binwidth=100)+ggtitle("Gene frequency before filterin")+
	theme(axis.text.x=element_text(size=10,angle=45,hjust=1))+
	geom_vline(xintercept=upper_feature,linetype="dotted",color="red")+xlab("Gene")
ggsave("noMito_gene_hist.pdf")

ggplot(noMito_obj@meta.data,aes(y=nCount_RNA,x=nFeature_RNA))+geom_point()+geom_vline(xintercept = upper_feature,colour="red")+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),axis.line=element_line((colour="black")))+
  xlab("genes")+ylab("UMI")
ggsave("noMito_umiCutoff.pdf")

pdf("vln_beforeFilt.pdf")
VlnPlot(noMito_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size=0,group.by="orig.ident")
dev.off()

ggplot(noMito_obj@meta.data,aes(y=nCount_RNA,x=nFeature_RNA,color=genotype))+geom_point()+geom_vline(xintercept = upper_feature,colour="red")+
	theme_bw()+theme(panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),axis.line=element_line((colour="black")))
ggsave("genotype_umiCutoff.pdf") 

ggplot(noMito_obj@meta.data,aes(y=nCount_RNA,x=nFeature_RNA,color=sex))+geom_point()+geom_vline(xintercept = upper_feature,colour="red")+
	theme_bw()+theme(panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),axis.line=element_line((colour="black")))
ggsave("sex_umiCutoff.pdf") 

ggplot(noMito_obj@meta.data,aes(y=nCount_RNA,x=nFeature_RNA,color=batch))+geom_point()+geom_vline(xintercept = upper_feature,colour="red")+
	theme_bw()+theme(panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),axis.line=element_line((colour="black")))
ggsave("batch_umiCutoff.pdf") 

ggplot(noMito_obj@meta.data,aes(y=nCount_RNA,x=nFeature_RNA,color=orig.ident))+geom_point()+geom_vline(xintercept = upper_feature,colour="red")+
	theme_bw()+theme(panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),axis.line=element_line((colour="black")))
ggsave("sample_umiCutoff.pdf") 


#######seurat object TO BE USED FOR FURTHER ANLAYSIS

filtr_obj=subset(noMito_obj, subset = nFeature_RNA < upper_feature & percent.mito<percent.mit& orig.ident !="NA12" )
filtr_obj$orig.ident=as.character(filtr_obj$orig.ident)

##mean no. of UMI and genes after filtering

filtr_obj@meta.data %>% group_by(orig.ident) %>% summarize(mean(nCount_RNA))
filtr_obj@meta.data %>% group_by(orig.ident) %>% summarize(mean(nFeature_RNA))

pdf("vln_afterFilt.pdf")
VlnPlot(filtr_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size=0,group.by="orig.ident")
dev.off()

#############################################################
#############################################################
#HOUSEKEEPING
#check size of the objects and remove if needed
#for(obj in ls()) {message(obj);print(object.size(get(obj)),units='auto')}
save.image("mid_analysis.RData")
#############################################################
#############################################################
#####1- sc transform - regularized negative binomial regression
#This code will calculate top 3000 variable features (default), scale.data matrices will have only variable features, 

#add the no. of cells from every sample to include in the regression analysis if the cell no. varies a lot. 
scx_df=mutate(scx_df, cell_no=case_when(
	orig.ident=="NA1" ~'21499',
	orig.ident=="NA2" ~ '7396',
	orig.ident=="NA3" ~ '24711',
	orig.ident=="NA4" ~ '27002',
	orig.ident=="NA5" ~ '25906',
	orig.ident=="NA6" ~ '23633',
	orig.ident=="NA7" ~ '19905',
	orig.ident=="NA8" ~ '16516',
	orig.ident=="NA9" ~ '14223',
	orig.ident=="NA10" ~ '11174',
	orig.ident=="NA11" ~ '2924',
	TRUE ~'NA'
	))

rownames(scx_df)=rownames(first_complete_obj@meta.data)	
cell_no=scx_df$cell_no
names(cell_no)=rownames(scx_df)

filtr_obj=AddMetaData(filtr_obj,cell_no,"cell_no")

# CAN BE AUTOMATED FROM HERE
fileName="integr" 
####1. regression 
#########option1 - SCtransform
filtr_obj=SCTransform(filtr_obj,vars.to.regress=c("percent.mito","batch","nCount_RNA","cell_no"),verbose=T,conserve.memory=T)
#log normalized counts are in filtr_obj[["SCT"]]@data

######option 2: linear regression
#normalization and find the most variable feature
filtr_obj=NormalizeData(filtr_obj,normalization.method = "LogNormalize")

#identify highest variable feature
filtr_obj=FindVariableFeatures(filtr_obj, selection.method = "vst", nfeatures = 2000)

#top 10 features and plotting
top20 =head(VariableFeatures(filtr_obj), 20)
plot1= VariableFeaturePlot(filtr_obj)
plot2=LabelPoints(plot = plot1, points = top20, repel = TRUE)
setEPS()
postscript(paste(fileName,"_variable.eps",sep=""))
plot2
dev.off()

### scale data regression
#all.genes = rownames(filtr_obj)
all.genes=head(VariableFeatures(filtr_obj), 2000) #top 2000
filtr_obj=ScaleData(filtr_obj, features = all.genes,model.use = "linear",vars.to.regress = c("percent.mito","batch","nCount_RNA","cell_no"))

#for both the option
####2.  pca visualization

filtr_obj=RunPCA(filtr_obj, features = VariableFeatures(object = filtr_obj))

#determine the dimensionality of the dataset using jackstraw 
filtr_obj= JackStraw(filtr_obj, num.replicate = 100,dims = 50)
filtr_obj= ScoreJackStraw(filtr_obj, dims = 1:50)
pdf(paste0(fileName,"_jackstraw.pdf"))
JackStrawPlot(filtr_obj, dims = 1:50)
dev.off()

#set dims based on pvals
pvals=sort(which(filtr_obj@reductions$pca@jackstraw$overall.p.values[,2]>0.05))
dimension=ifelse(length(pvals)==0,50,pvals[1]-1)

####3. cluster cells

filtr_obj<- FindNeighbors(filtr_obj, dims = 1:dimension)
save.image(paste0(fileName,"_cluster.RData"),compress=T)
  
for(iRes in seq(0.5,1,0.1)){ 
	print (iRes) 
	filtr_obj <- FindClusters(filtr_obj, resolution = iRes)
	#modulation - nstar-20,  niter500, force recalc T,epochs =30

	####4.UMAP
	filtr_obj <- RunUMAP(filtr_obj, dims = 1:dimension)
	#filtr_obj <- RunUMAP(filtr_obj, dims = 1:dimension,n.neighbors=50, min.dist=0.001)

	pdf(paste(fileName,iRes,"umap.pdf",sep="_"))
	#p1=DimPlot(filtr_obj, reduction = "umap", group.by = "orig.ident")
	p2=DimPlot(filtr_obj, reduction = "umap", label = TRUE)
	print(p2)
	dev.off()
	
	pdf(paste(fileName,iRes,"umap_genotype.pdf",sep="_"))
	print(DimPlot(filtr_obj,reduction="umap",split.by="genotype"))
	dev.off()
	
	hist_data=filtr_obj@meta.data[,c(7,12)] #7,10 for linear and 7,12 for sct
	head(hist_data)
	rownames(hist_data)=NULL
	hist_data_plot=melt(table(hist_data))
	hist_final=as.data.frame(group_by(hist_data_plot, seurat_clusters) %>% mutate(percent = value/sum(value)))
  hist_final=hist_final[order(hist_final$seurat_clusters),]
	head(hist_final)
	
	lim_plot=max(as.numeric(dspn_obj@meta.data$seurat_clusters))-1
	ggplot(hist_final,aes(x=seurat_clusters,y=percent,fill=genotype))+geom_bar(stat="identity")+ylab("% cells")+
	xlab("clusters")+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7))+
	scale_y_continuous(labels=scales::percent)+theme_bw()+scale_x_continuous(breaks=seq(0,lim_plot,1))
	ggsave(paste(fileName,iRes,"cellsperc_gen.pdf",sep="_"))


	#####5. differential expression
	cluster.markers =FindAllMarkers(filtr_obj,only.pos = T,logfc.threshold = 0.2,test.use="MAST") 
	write.table(cluster.markers,paste(fileName,iRes,"mast_deg.csv",sep="_"),sep="\t",row.names=F,quote = F)
	cluster.markers =FindAllMarkers(filtr_obj,only.pos = T,logfc.threshold = 0.2,test.use="wilcox") 
	write.table(cluster.markers,paste(fileName,iRes,"wilcox_deg.csv",sep="_"),sep="\t",row.names=F,quote = F)
	cluster.markers =FindAllMarkers(filtr_obj,only.pos = T,logfc.threshold = 0.2,test.use="roc") 
	write.table(cluster.markers,paste(fileName,iRes,"roc_deg.csv",sep="_"),sep="\t",row.names=F,quote = F)
}

save.image(paste0(fileName,"_allsample.RData"),compress=T)

#####differential exp by genotypes
Idents(filtr_obj)=filtr_obj$genotype
nam="bulk"

  foxp1VsCont=FindMarkers(filtr_obj,ident.1="Foxp1",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(foxp1VsCont,paste0(nam,"_p1VsCont_mast.txt"),sep="\t",quote = F)

  foxp2VsCont=FindMarkers(filtr_obj,ident.1="Foxp2",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(foxp2VsCont,paste0(nam,"_p2VsCont_mast.txt"),sep="\t",quote = F)

  foxp1Vsdoubl=FindMarkers(filtr_obj,ident.1="Foxp1",ident.2="double_cko",verbose = FALSE,test.use="MAST")
  write.table(foxp1Vsdoubl,paste0(nam,"_p1Vsdoubl_mast.txt"),sep="\t",quote = F)

  foxp2Vsdoubl=FindMarkers(filtr_obj,ident.1="Foxp2",ident.2="double_cko",verbose = FALSE,test.use="MAST")
  write.table(foxp2Vsdoubl,paste0(nam,"_p2Vsdoubl_mast.txt"),sep="\t",quote = F)

  doublVsCont=FindMarkers(filtr_obj,ident.1="double_cko",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(doublVsCont,paste0(nam,"_doublVsCont_mast.txt"),sep="\t",quote = F)

  foxp1Vsfoxp2=FindMarkers(filtr_obj,ident.1="Foxp1",ident.2="Foxp2",verbose = FALSE,test.use="MAST")
  write.table(foxp1Vsfoxp2,paste0(nam,"_p1Vsp2_mast.txt"),sep="\t",quote = F)

save.image("drd1_0.2_linear",compress=T)

######differential exp by genotypes within clusters
for(i in 0:21){ 
  nam=paste("cls",i,sep="_"); 
  print(nam)
  #assign(nam, subset(filtr_obj,seurat_clusters==i))
  cls_sub=subset(filtr_obj,seurat_clusters==i)
  Idents(cls_sub)=cls_sub$genotype
  
  foxp1VsCont=FindMarkers(cls_sub,ident.1="Foxp1",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(foxp1VsCont,paste0(nam,"_foxp1VsCont_deg.txt"),sep="\t",quote = F)

  foxp2VsCont=FindMarkers(cls_sub,ident.1="Foxp2",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(foxp2VsCont,paste0(nam,"_foxp2VsCont_deg.txt"),sep="\t",quote = F)

  foxp1Vsdoubl=FindMarkers(cls_sub,ident.1="Foxp1",ident.2="double_cko",verbose = FALSE,test.use="MAST")
  write.table(foxp1Vsdoubl,paste0(nam,"_foxp1Vsdoubl_deg.txt"),sep="\t",quote = F)

  foxp2Vsdoubl=FindMarkers(cls_sub,ident.1="Foxp2",ident.2="double_cko",verbose = FALSE,test.use="MAST")
  write.table(foxp2Vsdoubl,paste0(nam,"_foxp2Vsdoubl_deg.txt"),sep="\t",quote = F)

  doublVsCont=FindMarkers(cls_sub,ident.1="double_cko",ident.2="control",verbose = FALSE,test.use="MAST")
  write.table(doublVsCont,paste0(nam,"_doublVsCont_deg.txt"),sep="\t",quote = F)
}

#############################################################
#############################################################
```

