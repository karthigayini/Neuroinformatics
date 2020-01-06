For tmp-Seq data, follow the following steps

1. Trim the illumina adaptors and cut 3Bp from the 5' end of the read1 - as per email conversations with standford group and the reference paper 

   ```Bash
   module load trimgalore/0.4.1
   ls Sample_*/*R1_001.fastq.gz|sed "s/_R1_001.fastq.gz//g"|xargs -I % -P 12 sh -c 'trim_galore --paired --fastqc --illumina --clip_R1 3 --stringency 5  %_R1_001.fastq.gz %_R2_001.fastq.gz'
   ```

> The above steps does the following
>
> \#trim galore command will do the following
>
> SUMMARISING RUN PARAMETERS
>
> ==========================
>
> Input filename: Sample_Control_1/Control_1_S1_R1_001.fastq.gz
>
> Trimming mode: paired-end
>
> Trim Galore version: 0.4.1
>
> Cutadapt version: 1.9.1
>
> Quality Phred score cutoff: 20
>
> Quality encoding type selected: ASCII+33
>
> Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; user defined)
>
> Maximum trimming error rate: 0.1 (default)
>
> Minimum required adapter overlap (stringency): 5 bp
>
> Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
>
> All Read 1 sequences will be trimmed by 3 bp from their 5' end to avoid poor qualities or biases
>
> Running FastQC on the data once trimming has completed
>
> Output file(s) will be GZIP compressed



**I made a new folder called new_method and the results from this are stored there** 

2. Move the reports and the fq to the respective folders

3. convert fq to bam, sorted, indexed, duplicates are only marked using the script 

   ```bash
   #!/bin/bash
   
   #SBATCH -D /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq
   #SBATCH -N 1
   #SBATCH -t 80:00:00
   #SBATCH --job-name="fqToBam"
   #SBATCH -o /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq/%x%j.out
   #SBATCH -e  /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq/%x%j.err
   #SBATCH --mail-user=karthigayini.sivaprakasam@utsouthwestern.edu
   #SBATCH --mail-type=ALL
   
   cd /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq
   
   module load samtools/gcc/1.8 bwa/intel/0.7.17
   for i in `ls *_R1_001_val_1.fq.gz`; do
           fn=${i%_R1_001_val_1.fq.gz}
           bwa mem -t 12 -M /project/shared/bicf_workflow_ref/mouse/GRCm38/genome.fa $fn"_R1_001_val_1.fq.gz" $fn"_R2_001_val_2.fq.gz" |samtools view -bS - > ../bam/$fn".bam"
   done
   
   mkdir ../sortedBam
   module load samtools/1.6
   cd /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/bam
   for i in `ls *.bam`; do
          fn=${i%.bam}
           samtools view -bq 10 $i |samtools sort - > ../sortedBam/$fn"_sorted.bam"
           java -jar ~/tools/picard_2.18.29/picard.jar MarkDuplicates I=../sortedBam/$fn"_sorted.bam" O=../sortedBam/$fn"_noDup.bam" M=../sortedbam/$fn"_metrics.txt" REMOVE_DUPLICATES=false;
           samtools index ../sortedBam/$fn"_noDup.bam"
   done
   ```

   or use this code for parallelization inside the slurm script  ###this was used 

   ```bash
   #!/bin/bash
   
   #SBATCH -D /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq
   #SBATCH -N 1
   #SBATCH -t 80:00:00
   #SBATCH --job-name="fqToBam"
   #SBATCH -o /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq/%x%j.out
   #SBATCH -e  /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq/%x%j.err
   #SBATCH --mail-user=karthigayini.sivaprakasam@utsouthwestern.edu
   #SBATCH --mail-type=ALL
   #SBATCH --export=ALL
   #SBATCH -p super 
   
   #######################################
   #USUAGE: 
   #sbatch --export=ALL,SAMPLE=Control.bam ~/scripts/fqToBam_2.slurm 
   #######################################3
   
   echo "START_fq"
   date
   echo ${SAMPLE}
   sh /home2/s186964/scripts/fqToBam.sh ${SAMPLE}
   date
   echo "DONE"
   ```

   And the fqToBam.sh is 

   ```bash
   #!/bin/bash
   
   cd /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/trimmed_fq
   
   module load samtools/gcc/1.8 bwa/intel/0.7.1
   file=$1
   echo $file
   
   bwa mem -t 12 -M /project/shared/bicf_workflow_ref/mouse/GRCm38/genome.fa $file"_R1_001_val_1.fq.gz" $file"_R2_001_val_2.fq.gz" |samtools view -bS - > ../bam/$file".bam"
   
   cd /work/Neuroinformatics_Core/s186964/RMadbushi/round_2/new_method/bam
   samtools view -bq 10 $file".bam" |samtools sort - > ../sortedBam/$file"_sorted.bam"
   java -jar ~/tools/picard_2.18.29/picard.jar MarkDuplicates I=../sortedBam/$file"_sorted.bam" O=../sortedBam/$file"_noDup.bam" M=../sortedbam/$file"_metrics.txt" REMOVE_DUPLICATES=false
   ```

   4. Optional - merge the bams from the replicates to get the required read depth. 

      ```bash
      samtools merge Con_merged.bam 1Con_S1_noDup.bam 2Con_S6_noDup.bam 3Con_S11_noDup.bam
      ```
      
   5. Find the 5; end of each read - use Merged bam if the replicates were merged

      ```bash
      for i in `ls *_noDup.bam`; do fn=${i%_noDup.bam}; echo $fn; samtools view $i|awk '{print $3"\t"$4}' > $fn"_leftend.csv";done
      ```

   6. extend 20bp on both the sides of the this 5' or the leftmost position. Before that made a folder called beds and moved the .csv files into it. 

   ```bash
   for i in `ls *.csv`; do fn=${i%_leftend.csv}; echo $fn; awk 'BEGIN{OFS=FS="\t"}{print $1,$2-20,$2+20}' $i > $fn"_20bp.bed";done
   ```


7. To convert bam files to bigwile.

   ```bash
   samtools index xx.noDup.bam
   module load deeptools/2.5.0.1
   bamCoverage -b xx.noDup.bam -o xx.bw
   ```

   

###### ##########################

Quality control 

1. get the no. of reads mapped, and % duplicates - *.txt are the metrics file from mark duplicates. 

   ```bash
   for i in `ls *.txt`; do echo $i; grep -v '#' $i|head -3|tail -1|awk  '{print $4,$10}';done
   ```

2. Get the length vs MAPq quality - to check if there are very many fragments and the distribution of the scores - for 3 samples Control_3, gDNA1 and NMDA50-3

```bash
samtools view gDNA_1_S2.bam|awk '{print $5 "\t"length($10)}' > gDNA_1_length.csv
```

###### ####################

Getting the gtf file  in bed format with the required upstream and downstream window. 

- Grep for transcripts from GRCm38 gtf gencode version M20. And then get its start and end.

```bash
awk '{if($3=="transcript"){print $0}}' gencode.vM20.chr_patch_hapl_scaff.annotation.gtf |awk -F '[\t;]' 'BEGEIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$10,$12}' > grcm38_transcripts_length.bed
```

- Remove all those that do not start with chr and those start with chrM. Remove "transcript_id" and "gene_name" as well. In vi 

  ```vi
  :g/^chrM/d 
  :g!/^chr/d
  ```

- there are genes with multiple start codons. The multiple ones have different reading frames. (0,1,2,3) in such cases, we need to check if there is a UTR before the particular start region. In. Most of the cases, a UTR is present in the case with ‘0’ . Hence:

```bash
grep 'start_codon' gencode.vM20.chr_patch_hapl_scaff.annotation.gtf |awk -F '[\t;]' 'begin{OFS="\t"}{if($8==0) {print $1,$4,$5,$10,$12}}'  >grcm38_startcodon.bed
grep 'stop_codon' gencode.vM20.chr_patch_hapl_scaff.annotation.gtf  |awk -F '[\t;]' '{if($8==0) print $1,$4,$5,$10,$12}' > grcm38_stopcodon.bed
```

Then using Vi remove "transcript_id" and "gene_name" and codons in chrM and scaffolds with no names. 

- Using R get  

  ```R
  #######step1 :#finding the  longest transcript
  setwd("/Users/ksivaprakasam/Documents/RMadabushi/analysis/resources/")
  trans_table=data.frame(read.csv("grcm38_transcripts_length.bed",sep="\t",header=F),stringsAsFactors=F)
  head(trans_table)
  trans_table=trans_table[,-5]
  colnames(trans_table)=c("chr","start","end","length","transID","geneName")
  
  #get the transcript with max length and then its start and end location
  trans.agg=aggregate(length~geneName, data=trans_table, max) 
  trans.max=merge(trans.agg,trans_table,by.x=c("geneName","length"))
  
  #if you dim the the two DF we can see that after merging there are 2000 more rows. 54138 and 56138 transcripts
  # some genes have 2 or more transcripts with same length. Therefore remove the duplicates
  
  trans.max.uniq=trans.max[with(trans.max,order(geneName,-length)),] #order 
  trans.max.uniq=trans.max[!duplicated(trans.max[,c(1,2)]),]					#remove duplciates
  #trans.max.uniq final transcripts table
  
  #######################################################
  ##step 2: locate 3kbTSS1kb 1KBTES3kb for the longest transcripts
  tss=data.frame(read.csv("grcm38_startcodon.bed",sep="\t",header=F),stringsAsFactors=F)
  tss=tss[,-c(4,6)]
  colnames(tss)=c("chr","start","end","transID","geneName")
  
  tes=data.frame(read.csv("grcm38_stopcodon.bed",sep="\t",header=F),stringsAsFactors=F)
  tes=tes[,-c(4,6)]
  colnames(tes)=c("chr","start","end","transID","geneName")
  
  #merging the tes and tss with trans and then the both
  tss.max=merge(trans.max.uniq,tss,by="transID",all.x = T,suffixes = c("_trans","_ss"))
  tes.max=merge(trans.max.uniq,tes,by="transID",all.x = T,suffixes = c("_trans","_es"))
  codons=merge(tes.max[,c(1,2,7:9)],tss.max[,c(1,7:9)],by="transID",all = T)
  
  #identify the columns with NA in the chr and remove whole row. This way the trans with both start and stop are only retained reducing from 54621 transcripts to 20560
  codons_use=codons[complete.cases(codons),] 
  
  #Add 3k,1k,1k 3k bases
  #start codon
  codons_use$start_ss=codons_use$start_ss+1-3000
  codons_use$end_ss=codons_use$end_ss-1+1000
  
  #end codon
  codons_use$start_es=codons_use$start_es+1-1000
  codons_use$end_es=codons_use$end_es-1+3000
  
  write.table(codons_use,"codons_tobeused.csv",sep="\t",quote = F, row.names = F)
  #final table to be used. 
  ```

- ```
  awk 'BEGIN{OFS=FS="\t"}{print $3,$4,$5,$1":"$2}' codons_tobeused.csv > tes_3k.bed
  awk 'BEGIN{OFS=FS="\t"}{print $6,$7,$8,$1":"$2}' codons_tobeused.csv > tss_3k.bed
  vi to remove header
  ```

- Make bed files with a window id 100bp and sliding at 50b (window > sliding). Mathematically, 4000BP are separated into100bp windows = i.e 400. but they have 50bp sliding window = leading to 40x2 = 80 lines for every gene. 80x20560=1644800

  ```bash
  module load bedtools/2.29.0
  sort -k 1,1 -k2,2n tes_3k.bed |bedtools makewindows -b stdin -w 100 -s 50 -i src > downstrm_100bin_50bSliding.bed
  
  sort -k 1,1 -k2,2n tss_3k.bed |bedtools makewindows -b stdin -w 100 -s 50 -i src > upstream_100bin_50bSliding.bed
  ```

- change bam to bedPE. This script is named as BamtoBed.slurm

  ```bash
  #SBATCH -D /work/Neuroinformatics_Core/s186964/RMadbushi
  #SBATCH -N 1
  #SBATCH -t 100:00:00
  #SBATCH -o /work/Neuroinformatics_Core/s186964/RMadbushi/round_3/2019_10_03_CPL362_8541_0/raw/bamToBedPe/%j.out
  #SBATCH -e /work/Neuroinformatics_Core/s186964/RMadbushi/round_3/2019_10_03_CPL362_8541_0/raw//bamToBedPe/%j.err
  #SBATCH --export=ALL
  #SBATCH --mail-type=BEGIN,END,FAILED
  #SBATCH --mail-user=karthigayini.sivaprakasam@utsouthwestern.edu
  
  cd /work/Neuroinformatics_Core/s186964/RMadbushi/round_3/2019_10_03_CPL362_8541_0/raw/sortedBam
  
  module load samtools/1.6
  module load bedtools/2.29.0
  
  #sorted the bam file
  for i in `ls *noDup.bam`; do
          name=${i%.bam}
          samtools sort $i > ../bamToBedPe/$name"_sorted.bam"
          samtools index ../bamToBedPe/$name"_sorted.bam"
  done
  
  #sorting and intersecting with the 10bin reference file
  cd ../bamToBedPe/
  ls *sorted.bam|xargs -I {} -n 1 -P 3 sh -c 'file={};name=${file%_noDup_sorted.bam}; bedtools intersect -a ~/resources/downstrm_5kb_10bin.bed -b $file -f 0.7 -c >$name"_down_intersect.bed"'
  
  ls *sorted.bam|xargs -I {} -n 1 -P 3 sh -c 'file={};name=${file%_noDup_sorted.bam};bedtools intersect -a ~/resources/upstrm_5kb_10bin.bed -b $file -f 0.7 -c >$name"_up_intersect.bed"'
  
  #concatenate the intersected files by column - start from here. make awk print nth column
   paste *_up_intersect.bed| awk 'BEGIN{OFS=FS="\t"}{DL="";for(i=1;i<=NF;i+=(i<5?1:5)) {printf "%s%s",DL,$i;DL="\t"}; printf "\n"}' > up_all_intersects.bed
   
   paste *_down_intersect.bed|  awk 'BEGIN{OFS=FS="\t"}{DL="";for(i=1;i<=NF;i+=(i<5?1:5)) {printf "%s%s",DL,$i;DL="\t"}; printf "\n"}'  > down_all_intersects.bed
  ```

###### ##############

BACK TO R

Intersected files are normalized

```R
##############################################
###reading the aggregated matrix - intersect 10bp window with bamTobedPE 
###normalizing control and treatment signal counts with gDNA

module load R/3.6.1-gccmkl
setwd("/work/Neuroinformatics_Core/s186964/RMadbushi/round_3/2019_10_03_CPL362_8541_0/raw/bamToBedPe/R_ops")
up=data.frame(read.csv2("up_all_intersects.bed",header = F,sep="\t",stringsAsFactors = F))
down=data.frame(read.csv2("down_all_intersects.bed",header = F,sep="\t",stringsAsFactors = F))

colnames(up)=c("chr","start","end","id","Con","CPT","DRB","ETO","gen","Unx")
colnames(down)=c("chr","start","end","id","Con","CPT","DRB","ETO","gen","Unx")

#add one to all the values to avoid INF
down[,c(5:10)]=down[,c(5:10)]+1
up[,c(5:10)]=up[,c(5:10)]+1

#normalize the values with Unx (the control sample)
nc_up=sweep(up[,c(5:9)],1,FUN = "/",STATS = up$Unx)
nc_up=cbind(up[,c(1:4)],nc_up)
nc_down=sweep(down[,c(5:9)],1,FUN = "/",STATS = down$Unx)
nc_down=cbind(down[,c(1:4)],nc_down)

#basic reference matrix
write.table(nc_down,"NormalizedCount_down.csv",sep="\t",quote = F, row.names = F)
write.table(nc_up,"NormalizedCount_up.csv",sep="\t",quote = F, row.names = F)
##############################################
##############################################
###average counts at every 100bp window

#This is done by adding a column called group - which has the group ID, that is the no.odf lines for every gene. There are 80 lines for every gene(ID) and each line is 100BP. So in total there are 8000BP. 
#For every gene, counts in the first 100BP are averaged. This averaging goes on for the entire 8000Bp length of a gene. This concept is called binning here. 

nc_up$group=unlist(lapply(table(nc_up$id),seq.int))
nc_down$group=unlist(lapply(table(down$id),seq.int))

avg_100_up=aggregate(cbind(Con,CPT,DRB,ETO,gen)~group,nc_up,mean)
avg_100_down=aggregate(cbind(Con,CPT,DRB,ETO,gen)~group,nc_down,mean)

write.table(avg_100_down,"avg_100bp_down.csv",sep = "\t",row.names = F,quote = F)
write.table(avg_100_up,"avg_100bp_up.csv",sep = "\t",row.names = F,quote = F)

##############################################
####plotting
##############################################
library(reshape2)
library(ggplot2)
library(gridExtra)

lab=seq(-3,0.99,0.05) 
# the seq here is determined based on no. of rows in Avg_100. 
#the length region is 4kb, divided into 100bp region. But they have a sliding window of 50bp, making it (4000/100)x2
avg_100_up=cbind(avg_100_up,lab)
avg_100_up=avg_100_up[,-1]
long_avg_up=melt(avg_100_up,id="lab")

lab=seq(-0.99,3,0.05) 
avg_100_down=cbind(avg_100_down,lab)
avg_100_down=avg_100_down[,-1]
long_avg_down=melt(avg_100_down,id="lab")

#to exclude "gen" sample
long_avg_up=subset(long_avg_up,variable !="gen")
long_avg_down=subset(long_avg_down,variable !="gen")

up_plot=ggplot(long_avg_up,aes(x=lab,y=value,colour=variable))+
  theme_classic()+geom_smooth(method = "loess",span=0.1,se=F)+
  theme(legend.position = "none")+
  scale_x_continuous(name="Distance from TSS (KB)", limits=c(-3,0.95),breaks = c(-3,-2,-1,0,0.5))+
  scale_y_continuous(name="Normalized counts",limits=c(1,1.15),breaks=c(1,1.05,1.1,1.15))

down_plot=ggplot(long_avg_down,aes(x=lab,y=value,colour=variable))+
  theme_classic()+geom_smooth(method = "loess",span=0.1,se=F)+
  theme(plot.title = element_text(hjust=0.5),legend.position = c(0.7,0.9),legend.title = element_blank(),
        legend.background = element_rect(color="black",linetype = "solid",size=0.3))+
  scale_x_continuous(name="Distance from TES (KB)", limits=c(-0.99,2.95),breaks = c(-0.5,0,1,2,2.5))+
  scale_y_continuous(name="Normalized counts",limits=c(1,1.15),breaks=c(1,1.05,1.1,1.15),position="right")

pdf("avg_100.pdf")
grid.arrange(up_plot,down_plot,nrow=1,top="Normalized counts of treatment and control samples")
dev.off()

##############################################
##############################################
###smoothening the NC to calculate TI and TD

#since every line is 100BP, a modulo of 5 will take into account every 5th line, meaning every 5*100=500bp 
nc_down$gr2=(seq(nrow(nc_down))-1) %/%5
nc_up$gr2=(seq(nrow(nc_up))-1) %/%5 # -1 is needed, otherwise the first one will be only 4 rows and not 5

#calculate the average of running window of 500bp for the NCs
smooth_up=aggregate(cbind(Con,CPT,DRB,ETO,gen)~gr2,nc_up,mean)
smooth_down=aggregate(cbind(Con,CPT,DRB,ETO,gen)~gr2,nc_down,mean)

smooth_up=merge(smooth_up,nc_up,by="gr2")
smooth_up=smooth_up[!duplicated(smooth_up$gr2),c(2:10)]
smooth_up$end=smooth_up$start+500	#because we are averaging for every 500BP
colnames(smooth_up)[1:5]=c("Con","CPT","DRB","ETO","gen")
smooth_up=smooth_up[,c(6:9,1:5)]

smooth_down=merge(smooth_down,nc_down,by="gr2")
smooth_down=smooth_down[!duplicated(smooth_down$gr2),c(2:10)]
smooth_down$end=smooth_down$start+50

write.table(smooth_down,"smooth_down.csv",sep = "\t",row.names = F,quote = F)
write.table(smooth_up,"smooth_up.csv",sep = "\t",row.names = F,quote = F)

###################################################

#TI and TD calculation
#creat group 3 as row no. based on gene
ti_up=smooth_up[,c(4:9)]
avg_nc_up=aggregate(cbind(control_norm,nmda_norm)~id,ti_up,mean)
ti_down=smooth_down[,c(4:9)]
avg_nc_down=aggregate(cbind(control_norm,nmda_norm)~id,ti_down,mean)

merged_nc=merge(avg_nc_up,avg_nc_down,"id")
colnames(merged_nc)=c("id","cntrl_up","nmda_up","cntrl_down","nmda_down")
merged_nc$ti_cntrl=merged_nc$cntrl_up-merged_nc$cntrl_down
merged_nc$ti_nmda=merged_nc$nmda_up-merged_nc$nmda_down
merged_nc$td=merged_nc$ti_nmda-merged_nc$ti_cntrl

write.table(merged_nc,"ti.csv",sep="\t",row.names = F,quote = F)
#TD calculation remaining


```

######

#NGSPLOT

Installation

- Download NGSplot from Github. CD into the folder and unzip it. Install the required packages R. 
- Export the path as given by opening and editing bash_profile file
- Download the mm10 ngsplotdb.py form their google drive. Then use that and ngsplotdb.py install ref.tar.gz to install genome. 



IMP::: Rdata is saved as filename.RData. Change the filename after every iteration. Also , just to be safe, go to R and clear the env (R -> rm(list=ls()) )

To denormalize and get counts :is to multiply the RPM values by the library sizes and then divide by one million and you'll have the raw counts.

```R
#general workflow
module load R/3.2.1-intel
cat config.txt
Unx_merged_sorted.bam:Con_merged_sorted.bam     -1      Control
Unx_merged_sorted.bam:CPT_merged_sorted.bam     -1      CPT
Unx_merged_sorted.bam:DRB_merged_sorted.bam     -1      DRB
Unx_merged_sorted.bam:ETO_merged_sorted.bam     -1      ETO

~/tools/ngsplot-develop/bin/ngs.plot.r -G mm10 -R genebody -C config.txt -SS same -SE 0 -MW 2  -O pair_genebody
~/tools/ngsplot-develop/bin/ngs.plot.r -G mm10 -R tss -D refseq -Al bin -C config.txt -SE 0 -MW 5  -O Tss_normalized_mw5
~/tools/ngsplot-develop/bin/ngs.plot.r -G mm10 -R tss -L 5000 -D refseq -Al bin -C config.txt -SE 0 -MW 5  -O unnormalized_tss_mw5

```

Code for custom normalization - From discussions forum
 Make a copy of zip file (containing "avgprof.RData" and "avgprof.txt") in ngsplot output folder

```R
## Load Rdata file created by ngsplot preferably in an empty environment.
load("avgprof.RData")                           
ls()                                                                             
v.lib.size  # Vector containing total mapped reads per samples used for the default plot normalisation.
regcovMat   # "avgprof.txt" contains the output of this matrix.
ctg.tbl     # Table containing the bam file to sample name mappings.

## Multiply values in regcovMat by (MappedReads/1e6) to get raw counts.
raw.counts <- t(t(regcovMat)*(v.lib.size/1e6))  

## Create a vector containing the custom normalisation factors you wish to use per sample.
## The length and ordering of this has to match that of the sample names in the columns of regcovMat.
norm.factors <- c(0.06058716,0.06058716,0.06058716,0.06058716)  
                                                                                        
regcovMat <- t(t(raw.counts)/(norm.factors)) 
## Output matrix created using custom normalisation factors. 
## Have to name it "regcovMat" because replot.R will use it to plot profiles again with custom normalisation factors.

write.table(regcovMat, file="avgprof.txt", row.names=F, sep="\t", quote=F)     
save(list=ls(all=TRUE),file="avgprof.RData")                              
## Overwrite avgprof.RData to contain custom normalised regcovMat.

## Zip the folder containing these files
## Run replot.r on the new zip file
~/tools/ngsplot-develop/bin/replot.r prof OR heatmap -I normalized_mw5.zip -O renormed
```

The norm.factors were decided by running ngsplot for only the control sample. Followed by loading the avgprof.RData. Raw counts for that sample was calculated as given in the above code and the median was obtained. Since the mean and median were almost same, mean was used. Range was 0.04503140 0.07258116

Even after re-normalising the ratios remained in negative. I suspect that vlibsize is where the error could be. 
