
###################################################
###################################################
#### Written by Karthigayini Sivaprakasam
#### 2020-1-6
###################################################
###################################################

######################################################
####Deseq2 analysis and generation of plots
#reference :https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/#
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#deseq2-import-functions
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq
######################################################
#####DESEQ2
library(DESeq2)
directory="/set/the/directory" # set the working directory
setwd(directory)
outputPrefix="DE_nmda" # the analysis title

#sampleTable: tab delimited file with 3 columns -sampleName, Filename with extension and condition
sampleTable=data.frame(read.table("sampleTable.csv",header=F,stringsAsFactors=F,sep=","))
colnames(sampleTable)=c("sampleName","fileName","condition")

treatments=c("NMDA","Control")
treatments=c("HET","WT")
treatments=c("Emx_CRE","control")
ddsHTseq=DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design = ~ condition)
colData(ddsHTseq)$condition=factor(colData(ddsHTseq)$condition,levels=treatments)

dds <- DESeq(ddsHTseq)
res <- results(dds)

## transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#########plot PCA of the data to see which samples to use. 
#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))

# set condition
scores <- data.frame(pc$x, sampleTable$condition)
colnames(scores)[9]="condition"

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5) 
  +geom_text(label=rownames(scores),hjust=0.4, vjust=-1.2,size=3)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "_ggplot2.pdf"))

####################
#removing genes with normalized counts less than 5 in less than 3 samples
est_dds=estimateSizeFactors(dds)
idx <- rowSums( counts(est_dds, normalized=TRUE) >= 5 ) >= 3
dds_filter <- dds[idx,]
#or genes that have a SUM of 10 reads in total 
#keep=rowSums(counts(dds))>=5, ds_1=dds[keep,]

#remove the samples based on pCA - for ctx_het
#based on the plot ctx21 nad 8 are removed, ctx_ko= remove 18_emx and 11_ctr
dds_filter=dds[,-c(1,8)] #dds_filter=dds[,-c(1,7)]
dds_filtered=DESeq(dds_filter)
##############
res <- results(dds_filtered)
summary(res)

#write the resutls
res_all=merge(as.data.frame(res), as.data.frame(counts(dds_filtered,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(res_all)[1] <- 'gene'
write.csv(res_all, file = paste0(outputPrefix, "_all_results-with-normalized.csv"),row.names = F)


#order the results by padj value
res_05= subset(res, padj<0.05)
res_05 <- res_05[order(res_05$padj),]
resdata <- merge(as.data.frame(res_05), as.data.frame(counts(dds_filtered,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "_results-with-normalized.csv"),row.names = F)

#for GSEA
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_gsea_normalized_counts.txt"), sep = '\t',quote = F)

mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "_test-conditions.txt"))

# heatmap of data
library("RColorBrewer")
library("pheatmap")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds_filtered,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
pdf(file = paste0(outputPrefix,"_data_heatmap.pdf"))
#my heatmap
pheatmap(as.matrix(resdata[,c(8:13)]),cluster_cols = F,color = my_palette,scale = "row",
         main = "padj genes_heatmap")
dev.off()
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")

#####################################################################
#optional analysis

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))


####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
pdf(file = paste0(outputPrefix,"_maplot.pdf"))
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.off()

#####plot the counts of gene with lowest pvalue 
plotCounts(dds_filtered,gene=which.min(res$padj), intgroup="condition")

###########
###heatmap of the count matrix
rld_filtered=rlogTransformation(dds_filtered, blind=T)
select <- order(rowMeans(counts(dds_filtered,normalized=TRUE)),decreasing=TRUE)[1:40]
pheatmap(assay(rld_filtered)[select,],cluster_cols = F,cluster_rows = F,show_rownames = F)

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds_filtered), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf(file = paste0(outputPrefix,"_heatmap.pdf"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
