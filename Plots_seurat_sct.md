Generate Plots from the previously analyzed Single cell RNA seq data

Requirement:Already analysed scRNA data using Seurat/3.0.0. The RData file should have all the meta data and the final complete seurat object. The analysis could be with either ScT method or other regressions. 

If other regression forms are used, double check the column nos. 

This code is modifed to fit the NA-2019 project

```R
###################################################
###################################################
#### Written by Karthigayini Sivaprakasam
#### 2019-12-30
###################################################
###################################################

module load Seurat/3.0.0 #works using R/3.5.2

library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(Seurat)

fileName="integr"
load(paste0(fileName,"allsamples.RData"))

#############################################################
#############################################################
#### heatmap with top 5 DEG from every cluster
my_levels=as.factor(seq(0,max(as.numeric(filtr_obj$seurat_clusters))-1))  #reorder cluster nos
filtr_obj@active.ident=factor(filtr_obj@active.ident,levels=my_levels)

cluster.markers$cluster=factor(cluster.markers,levels=my_levels)
top_5=cluster.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
DoHeatmap(filtr_obj,features=top_5$gene,size=3,angle=90,hjust=0.5)+NoLegend()
ggsave("top5_clusterwise_res1.1.pdf")

#############################################################
#############################################################
##########VIOLIN PLOTS##########

#general plot from Seurat 
pdf("vlnplot.pdf")
VlnPlot(filtr_obj,features=c("nFeature_RNA","nCount_RNA"),ncol=1,group.by="genotype",slot="counts",pt.size=0,log=F)
dev.off()

#### Violin plot here is generated between two genotypes - control and experiment
foxp1_sub=subset(filtr_obj,genotype=="Foxp1" | genotype=="control")
unique(foxp1_sub$genotype)

pdf(paste0(fileName,"_upset_p1_fox.pdf"))
plots=VlnPlot(foxp1_sub, features = c("Foxp1","Foxp2"),slot="counts",ncol=1,pt.size=0,log=F, split.by="genotype",combine=F)
CombinePlots(plots=plots,legend="right",ncol=1)
dev.off()

#for the same set of genotypes, violin plots are generated for two genes
pdf(paste(fileName,iRes,"vln_p1_drd.pdf",sep="_"))
plots=VlnPlot(foxp1_sub, features = c("Drd","Drd2"),slot="counts",ncol=1,pt.size=0,log=F, split.by="genotype",combine = F)
CombinePlots(plots=plots,legend="right",ncol=1)
dev.off()

##### custom violin plots for genes, UMIs and percent mito content
g=ggplot(filtr_obj,aes(x=genotype,y=nFeature_RNA,fill=genotype))+geom_violin()+scale_fill_manual(values=c("black","red","#008000"))+
ylab("count")+ggtitle("Genes")+theme(legend.position="none",axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=12))

u=ggplot(filtr_obj,aes(x=genotype,y=nCount_RNA,fill=genotype))+geom_violin()+scale_fill_manual(values=c("black","red","#008000"))+
ylab("count")+ggtitle("UMI")+theme(legend.position="none",axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=12))

p=ggplot(filtr_obj,aes(x=genotype,y=percent.mito,fill=genotype))+geom_violin()+scale_fill_manual(values=c("black","red","#008000"))+
ylab("count")+ggtitle("Percent Mitochondrial Content")+
theme(axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.x=element_blank(),plot.title=element_text(size=12),legend.title=element_text(size=10),legend.text=element_text(size=10))

plot_grid(g,u,p,ncol=2,nrow=2)
ggsave("vlnplot_style.pdf")

## Sample-wise
p1=ggplot(filtr_obj@meta.data,aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+geom_violin()+theme_classic()+xlab("Sample")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=orig.ident,y=nCount_RNA,fill=orig.ident))+geom_violin()+theme_classic()+xlab("Sample")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p3=ggplot(filtr_obj@meta.data,aes(x=orig.ident,y=percent.mito,fill=orig.ident))+geom_violin()+theme_classic()+xlab("Sample")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("vln_",fileName,"_sample.eps"))  plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),p3+theme(legend.position="none"),labels="AUTO",nrow=1)
dev.off()

## Batch-wise
p1=ggplot(filtr_obj@meta.data,aes(x=batch,y=nFeature_RNA,fill=batch))+geom_violin()+theme_classic()+xlab("Batch")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=batch,y=nCount_RNA,fill=batch))+geom_violin()+theme_classic()+xlab("Batch")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p3=ggplot(filtr_obj@meta.data,aes(x=batch,y=percent.mito,fill=batch))+geom_violin()+theme_classic()+xlab("Batch")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("vln_",fileName,"_batch.eps"))
plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),p3+theme(legend.position="none"),labels="AUTO",nrow=1)
dev.off()

## sex-wise
p1=ggplot(filtr_obj@meta.data,aes(x=sex,y=nFeature_RNA,fill=sex))+geom_violin()+theme_classic()+xlab("sex")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=sex,y=nCount_RNA,fill=sex))+geom_violin()+theme_classic()+xlab("sex")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p3=ggplot(filtr_obj@meta.data,aes(x=sex,y=percent.mito,fill=sex))+geom_violin()+theme_classic()+xlab("sex")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("vln_",fileName,"_sex.eps"))
plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),p3+theme(legend.position="none"),labels="AUTO",nrow=1)
dev.off()

## genotype-wise
p1=ggplot(filtr_obj@meta.data,aes(x=genotype,y=nFeature_RNA,fill=genotype))+geom_violin()+theme_classic()+xlab("genotype")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90))
p2=ggplot(filtr_obj@meta.data,aes(x=genotype,y=nCount_RNA,fill=genotype))+geom_violin()+theme_classic()+xlab("genotype")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90))
p3=ggplot(filtr_obj@meta.data,aes(x=genotype,y=percent.mito,fill=genotype))+geom_violin()+theme_classic()+xlab("genotype")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90))

setEPS()
postscript(paste0("vln_",fileName,"_genotype.eps"))
plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),p3+theme(legend.position="none"),labels="AUTO",nrow=1)
dev.off()

#############################################################
#############################################################
##########SCATTER PLOTS##########

## sample-wise
p1=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=orig.ident))+geom_point(shape=19)+theme_classic()+xlab("Sample")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=percent.mito,color=orig.ident))+geom_point(shape=19)+theme_classic()+xlab("Sample")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("scatter_",fileName,"_sample.eps"))
plot_grid(p1+theme(legend.position="none"),p2,labels="AUTO",nrow=1)
dev.off()

## sex-wise
p1=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=sex))+geom_point(shape=19)+theme_classic()+xlab("Sex")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=percent.mito,color=sex))+geom_point(shape=19)+theme_classic()+xlab("Sex")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("scatter_",fileName,"_sex.eps"))
plot_grid(p1+theme(legend.position="none"),p2,labels="AUTO",nrow=1)
dev.off()

## batch wise
p1=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=batch))+geom_point(shape=19)+theme_classic()+xlab("batch")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=percent.mito,color=batch))+geom_point(shape=19)+theme_classic()+xlab("batch")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("scatter_",fileName,"_batch.eps"))
plot_grid(p1+theme(legend.position="none"),p2,labels="AUTO",nrow=1)
dev.off()

## genotype-wise
p1=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=genotype))+geom_point(shape=19)+theme_classic()+xlab("genotype")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggplot(filtr_obj@meta.data,aes(x=nCount_RNA,y=percent.mito,color=genotype))+geom_point(shape=19)+theme_classic()+xlab("genotype")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

setEPS()
postscript(paste0("scatter_",fileName,"_genotype.eps"))
plot_grid(p1+theme(legend.position="none"),p2,labels="AUTO",nrow=1)
dev.off()
#############################################################
#############################################################

#####plot top 5 and the least 5 variable genes
top10=head(VariableFeatures(filtr_obj),10)
new_obj=filtr_obj
new_obj$seurat_clusters=as.numeric(new_obj$seurat_clusters)

pdf("vln_top5.pdf",width=8,height=16)
 VlnPlot(new_obj,features=top10[1:5],ncol=1,pt.size=0.1,log=T, slot="counts")
dev.off()

pdf("vln_tail5.pdf",width=8,height=16)
VlnPlot(new_obj,features=top10[6:10],ncol=1,pt.size=0.4,log=T, slot="counts")
dev.off()
#############################################################
#############################################################
################CELL PERCENTAGE PLOTS################

aggregate(nCount_RNA ~genotype, range_mat,sum)

hist_data=filtr_obj@meta.data[,c(1,9)] 
# to get the 5,6,7 instead of 1 to creat plots batch-wise, sex-wise and genotype-wise
head(hist_data)
rownames(hist_data)=NULL
hist_data_plot=melt(table(hist_data))
hist_final=as.data.frame(group_by(hist_data_plot, seurat_clusters) %>% mutate(percent = value/sum(value)))
hist_final=hist_final[order(hist_final$seurat_clusters),]
head(hist_final)

ggplot(hist_final,aes(x=seurat_clusters,y=percent,fill=orig.ident))+geom_bar(stat="identity")+ylab("% cells")+xlab("clusters")+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=7))+scale_y_continuous(labels=scales::percent)+theme_bw()
ggsave(paste0(fileName,"_cells_perc_sample.pdf"))

#############################################################
#############################################################
################UPSET PLOTS################

ibrary(UpSetR)

#read the DEG tables generated from Seurat workflow

doubl=data.frame(read.table("bulk_no0-6-21_doublVsCont_deg.txt",sep="\t",header=T))
foxp1=data.frame(read.table("bulk_no0-6-21_foxp1VsCont_deg.txt",sep="\t",header=T))
foxp2=data.frame(read.table("bulk_no0-6-21_foxp2VsCont_deg.txt",sep="\t",header=T))

doubl$gene=rownames(doubl)
foxp1$gene=rownames(foxp1)
foxp2$gene=rownames(foxp2)

apply(foxp2[,c(1:5)],2,range)
#did this for all the 3 data frames and saw that the Pval and FDR are within range But not Pcts.  

#so based on Idents (1 and 2) choose which are up and downregulated. 
down_doubl=doubl[which(doubl$avg_logFC<0),6]
up_doubl=doubl[which(doubl$avg_logFC>0),6]

down_foxp1=foxp1[which(foxp1$avg_logFC<0),6]
up_foxp1=foxp1[which(foxp1$avg_logFC>0),6]

down_foxp2=foxp2[which(foxp2$avg_logFC<0),6]
up_foxp2=foxp2[which(foxp2$avg_logFC>0),6]

#meta_data=data.frame(sets=c("up_doubl","up_foxp1","up_foxp2","down_doubl","down_foxp1","down_foxp2"),
#reg=c("upregulated","upregulated","upregulated","downregulated","downregulated","downregulated "))
listInput=list(up_doubl,up_foxp1,up_foxp2,down_doubl,down_foxp1,down_foxp2)
names(listInput)=c("Foxp1/2 up","Foxp1 up","Foxp2 up","Foxp1/2 dn","Foxp1 dn","Foxp2 dn")

pdf("upset_all_mast_noFilt.pdf")
upset(fromList(listInput), order.by="freq",nsets=6)
#sets.bar.color=c("pink","green","pink","green","pink","green"),sets.x.label="no. of DEGs",matrix.color="blue",sets=c("up_doubl","down_doubl","up_foxp1","down_foxp1","up_foxp2","down_foxp2"),keep.order=T)
dev.off()
#note if u use keep.order, the intersection happens by exclusivity. 

#############################################################
#############################################################
################CLUSTER MARKERS################

plots=FeaturePlot(filtr_obj,features=c("Aqp4","Gja1","Mobp","Mbp","Olig1","Olig2","Cx3cr1","Ppp1r1b"),combine=F)
#plots=FeaturePlot(filtr_obj,features=c("Drd","Drd2","Adarb2","Tac1","Foxp1","Foxp2","Isl1","Dscam","Ebf1"),combine=F)
CombinePlots(plots=plots,legend="right",ncol=3)
ggsave("cliuster_markers.pdf",dpi="retina")

#############################################################
#############################################################
################UMAP MODIFICATION################

#Mention specific color for every cell type. 

#first rename the culster nos with names
new_cluster.ids=c("iSPNs","dSPNs","dSPNs","eSPNs","oligodendrocytes","dSPNs","eSPNs","oligodendrocytes","endothelial","dSPNs","oligodendrocytes",
"oligodendrocytes","iSPNs","eSPNs","eSPNs","iSPNs","oligodendrocytes","astrocytes","eSPNs","eSPNs","dSPNs",
"dSPNs","oligodendrocytes","iSPNs","oligodendrocytes","oligodendrocytes","eSPNs","eSPNs","endothelial","oligodendrocytes","interneurons",
"eSPNs","oligodendrocytes","dSPNs")

f=levels(filtr_obj)
names(new_cluster.ids)=sort(as.numeric(as.character(f)))
filtr_obj=RenameIdents(filtr_obj,new_cluster.ids)

pdf(paste0(fileName,"_umap_cellColors.pdf")) 										
p2=DimPlot(filtr_obj,reduction="umap",cols=c('dSPNs'="dark green",'iSPNs'="Red",'eSPNs'="Purple",'oligodendrocytes'="grey",'endothelial'="dark blue",'interneurons'="pink",'astrocytes'="orange"))
p2
dev.off()


```

