#### RNASEQ workflow

```
###################################################
###################################################
#### Written by Karthigayini Sivaprakasam
#### 2020-1-6
###################################################
###################################################
```

STEPS:

1. ALIGNMENT

2. QUALITY METRICS

3. FILTERING

4. GET COUNTS

5. DIFFERENTIAL EXPRESSION

\############################################

##### **REFERENCE FILES**

1. **GENCODE ANNOTATION** - download the comprehensive annotation and grep for protein coding genes and 

those not in chrM,X or Y

2. **rRNA BED FILE**: get it from UCSC table browser 

- Select "All Tables" from the group drop-down list
- Select the "rmsk" table from the table drop-down list
- Choose "bed" as the output format
- Type a filename in "output file" so your browser downloads the result
- Click "create" next to filter
- Next to "repClass," type rRNA
- Next to free-form query, select "OR" and type repClass = "5S"
- Click submit on that page, then get output on the main page

the above bed will give a 4 column Bed file- but Rseqc need a 12 columns. Hence:

```bash
awk 'BEGIN{OFS=FS="\t"}{print $0,$2,$3,"0","1",$3-$2,"0"}' mm10_rRNA.bed >mm10_rRNA_rseqc.bed
```

3. **KNOWN GENE LIST**: go to downloads page in ucsc -> select the organism -> annotations database -> download knowngene.txt.gz

   

##### **WORKFLOW**

1. Run the alignment through ASTROCYTE pipeline and get the BAMs

2. Get the stats for the bams 

```bash
for i in `ls */bams/*.bam`; do fn=`basename -s .bam $i`; echo $fn; java -jar ~/tools/picard_2.18.29/picard.jar CollectRnaSeqMetrics I=$i O=$fn"_picard.txt" REF_FLAT= ~/reference_genomes/mm10_Grch38/gtf_gff/refFlat.txt.gz STRAND=FIRST_READ_TRANSCRIPTION_STRAND;done

module load samtools/1.6
ls *.bam|xargs -P 14 -I % sh -c 'samtools stats %|grep ^SN > %_samStats.txt'
```

3. Filtering BAMs

-   Filter for primary alignments only - remove unmapped and chimeric

```bash
module load samtools/1.6
ls *.bam| xargs -n 1 -P 12 -I % sh -c 'samtools view -F 256 -b % > %_primary.bam;' 
```

- Split Primary alignment in rRNA (in.bam) and non rRNA (ex.bam).

```bash
module load RSeQC/2.6.4

ls *primary.bam |xargs -n 1 -P 12 -I % sh -c 'split_bam.py -i % -r /home2/s186964/reference_genomes/mm10/gtf_gff/mm10_rRNA_rseqc.bed -o %".rrna"'
rm *.junk.bam
rm *.in.bam
```

-  Get uniquely mapped

```bash
ls *ex.bam|xargs -I % -n 1 -P 12 sh -c '(samtools view -H %; samtools view -F 2308 % | grep -w 'NH:i:1') | samtools view -bS - > %".uniq.bam"'
```

-  index the uniq bams

```bash
for i in `ls *uniq.bam`; do echo $i; samtools index $i $i.bai; done
rename _primary.bam.rrna.ex.bam.uniq _primary_rrnaEx_uniq *_primary.bam.rrna.ex.bam.uniq*
```

- count the no. of mapped reads

```bash
 samtools view -c -F 4 $i in a for loop
```

​	4. Generate counts

```bash
#!/bin/bash

#SBATCH -N 1
#SBATCH -t 20:00:00
#SBATCH --mail-type=BEGIN,END,FAILED
#SBATCH --mail-user=karthigayini.sivaprakasam@utsouthwestern.edu
#SBATCH -o /work/Neuroinformatics_Core/s186964/klab_rna/%j_klab_htseq.log
#SBATCH -e /work/Neuroinformatics_Core/s186964/klab_rna/%j_klab_htseq.err

cd /work/Neuroinformatics_Core/s186964/klab_rna/ctx_ko/bams

module load HTSeq/0.6.1
ls *uniq.bam|xargs -P 8 -I % sh -c 'htseq-count -f bam -t exon -m intersection-strict -i gene_name -s reverse % /home2/s186964/reference_genomes/mm10/gtf_gff/gencode.vM21.annotation_protCod_noMXY_2.gtf > %"_htseq.txt"'
```

5. Generate DE using deseq2 - script is in the deseq.slurm

**NOTE:** If by any chance HTSeq does not return any results, there might be a strand issue. So, just check it using the following command.

```bash
infer_experiment.py -r /home2/s186964/reference_genomes/mm10/gtf_gff/download_GRCm38_mm10_RefSeq_2.bed -i CTX_KO_S11.bam_primary_rrnaEx_uniq.bam

Reading reference gene model /home2/s186964/reference_genomes/mm10/gtf_gff/download_GRCm38_mm10_RefSeq_2.bed ... Done

Loading SAM/BAM file ... Total 200000 usable reads were sampled


This is SingleEnd Data
Fraction of reads failed to determine: 0.0087
Fraction of reads explained by "++,--": 0.0074
Fraction of reads explained by "+-,-+": 0.9839
```

For this the bed file was downloaded from the same site (HTSeq) and the chr was replaced with nothing since the bam file does not have chr (noted from head of the bam)

The above result tells that the strand is reversed. 

