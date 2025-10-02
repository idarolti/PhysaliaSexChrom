# Day 4 Practical

This practical will cover:

1. Differential gene expression analysis
2. Detecting Y gene activity decay and dosage compensation

## 00. Prepare work folder for day 4

```
mkdir day4
cd day4
conda activate sexchr
```

## 01. Differential gene expression analysis

For this part of the practical, we will use an example for quantifying gene expression in the absence of a GTF file or transcriptome reference. Several tools support genome-guided quantification, where RNA-seq reads are aligned to the genome and expression is estimated from assembled transcript structures. Here, we will focus on a pipeline based on HISAT2 and StringTie, and the model system is the willow (Salix viminalis), a dioecious species for which we have male and female RNA-seq data from both reproductive tissue (catkin) and somatic tissue (leaf).

### Map RNA-seq reads

Mapping reads can be done with **[HISAT2](https://daehwankimlab.github.io/hisat2/)**, a fast and sensitive alignment program for mapping next-generation sequencing reads. 

Generate genome index.

```
mkdir differential_gene_expression
cd differential_gene_expression
mkdir genome
cd genome
cp /home/ubuntu/Share/day4/willow/genome/genome_assembly_1k.fa ./
hisat2-build -f genome_assembly_1k.fa genome_assembly_1k
```

Align paired-end reads to the genome, then sort.

```
cd ..
cp -r /home/ubuntu/Share/day4/willow/reads/ ./
mkdir hisat
cd hisat

hisat2 ../genome/genome_assembly_1k -1 ../reads/female1_catkin_R1.fastq -2 ../reads/female1_catkin_R2.fastq -q --no-discordant --no-mixed --no-unal --dta -S female1_catkin.sam

samtools view -bS female1_catkin.sam | samtools sort -o female1_catkin_coordsorted.bam
rm female1_catkin.sam
```

HISAT2 options:

--no-discordant: suppress discordant alignments for paired reads

--no-mixed: suppress unpaired alignments for paired reads

--no-unal: suppresses the output of SAM records for reads that fail to align

--dta/--downstream-transcriptome-assembly: report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computationa and memory usage.

Run hisat2 and sorting for all the samples.

```
reads_dir="../reads"
# Genome index prefix
genome_index="../genome/genome_assembly_1k"

# Loop over all R1 fastq files
for r1 in ${reads_dir}/*_R1.fastq; do
    # Extract base name (e.g., female1_catkin)
    base=$(basename "$r1" "_R1.fastq")
    r2="${reads_dir}/${base}_R2.fastq"
    sam="${base}.sam"
    bam_coordsorted="${base}_coordsorted.bam"

    # Run HISAT2 alignment
    hisat2 "$genome_index" -1 "$r1" -2 "$r2" -q --no-discordant --no-mixed --no-unal --dta -S "$sam"

    # Convert SAM to sorted BAM
    samtools view -bS "$sam" | samtools sort -o "$bam_coordsorted"

    # Optional: remove intermediate SAM file to save space
    rm "$sam"
done
```

### Extract gene coordinates

**[StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)** - A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

```
cd ../
mkdir stringtie
cd stringtie
mkdir subset
cd subset

stringtie ../../hisat/female1_catkin_coordsorted.bam -o female1_catkin.gtf -A female1_catkin.gene_abund
```

Run StringTie for all samples

```
# Directory containing BAM files
bam_dir="../../hisat"

# Loop over all sorted BAM files in that directory
for bam in ${bam_dir}/*_coordsorted.bam; do
    # Get the base sample name, removing path and suffix
    base=$(basename "$bam" "_coordsorted.bam")
    
    # Run StringTie with output files named by sample
    stringtie "$bam" -o "${base}.gtf" -A "${base}.gene_abund"
done
```

Make list of all gtf files and merge

```
ls *.gtf > gtfs.list

stringtie --merge gtfs.list -o merged.gtf
```

### Obtain read counts 

We will use **[HTSeq]([http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](https://htseq.readthedocs.io/en/release_0.11.1/count.html))**. For faster processing, bam files must be sorted by read name, so paired reads appear adjacent in the file and are counted more quickly (default option in HTSeq). If your BAM is coordinate sorted (sorted by genomic position, the default with samtools sort), you must specify -r pos to tell htseq-count the BAM is position sorted. However, this requires more memory, as mate pairs may be separated and htseq-count buffers them until it sees the mate.

**(DO NOT RUN)**

```
samtools sort -n -T namesorted -O bam -o ../hisat/female1_catkin_namesorted.bam female1_catkin_coordsorted.bam
htseq-count -f bam -r name -s no female1_catkin_namesorted.bam ../stringtie/merged.gtf > female1_catkin_htseqcount.txt
```

Copy the htseq-count output on the full dataset and copy scripts for downstream analyses.

```
cd ../../
cp -r /home/ubuntu/Share/day4/willow/scripts ./
mkdir htseq
cp -r /home/ubuntu/Share/day4/willow/htseq/catkin ./htseq/
cp -r /home/ubuntu/Share/day4/willow/htseq/leaf ./htseq/
```

Merge read count outputs per tissue. Below is an example for catkin. Run the same for the leaf tissue.

```
cd scripts
python3 extract-counts.py ../htseq/catkin ../htseq/catkin/read_counts_catkin.txt
```

### Filter lowly expressed genes

Copy the merged.gtf file based on the full dataset, then create a file with the length of each gene. You only need to run this using one tissue, as all genes are included in both.

```
cd ../stringtie
mkdir fullset
cp /home/ubuntu/Share/day4/willow/stringtie_gtfs/fullset/merged.gtf ./stringtie/fullset/
cd scripts
python3 extract-gene-lengths.py ../htseq/catkin/read_counts_catkin.txt ../stringtie/fullset/merged.gtf ../stringtie/fullset/gene_length.txt
```

Download the read counts file for each tissue and the gene_length txt to your local directory. Use edgeR to convert counts to RPKM.

```
library("edgeR")

data <- read.table("read_counts.txt",stringsAsFactors=F,header=T, row.names=1)
names(data)
dim(data)

expr <- DGEList(counts=data)
gene_length <- read.table("gene_length.txt",stringsAsFactors=F)
head(gene_length)
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE
rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector)
write.table(rpkm, file="rpkm.txt",quote=F, sep="\t")
```

Filter genes that do not have a minimum of 2 RPKM expression in at least half of the individuals of one sex.

```
python filter-expression.py rpkm_catkin.txt read_counts_catkin.txt rpkm_catkin_filtered.txt read_counts_catkin_filtered.txt 'F,F,F,M,M,M'
```

### Normalization of gene expression with edgeR

```
library("edgeR")

data <- read.table("read_counts_catkin.txt",stringsAsFactors=F,header=T, row.names=1)
head(data)
dim(data)
conditions <- factor(c("F","F","F","M","M","M"))
length(conditions)

#Check raw read count data
expr <- DGEList(counts=data,group=conditions)
plotMDS(expr,xlim=c(-6,6))
```

<img width="643" height="623" alt="Screenshot 2025-10-01 at 19 19 32" src="https://github.com/user-attachments/assets/6a38a5da-05f2-4e37-b26a-3f1edbba5d08" />


```
cpm_expr <- cpm(expr)
sample1 <- density(log2(cpm_expr[,1]))
sample2 <- density(log2(cpm_expr[,2]))
sample3 <- density(log2(cpm_expr[,3]))
sample4 <- density(log2(cpm_expr[,4]))
sample5 <- density(log2(cpm_expr[,5]))
sample6 <- density(log2(cpm_expr[,6]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Raw log2 cpm",col="red")
lines(sample2, type="l",lwd=2,col="red")
lines(sample3, type="l",lwd=2,col="red")
lines(sample4, type="l",lwd=2,col="blue")
lines(sample5, type="l",lwd=2,col="blue")
lines(sample6, type="l",lwd=2,col="blue")
```

<img width="617" height="587" alt="Screenshot 2025-10-01 at 19 23 46" src="https://github.com/user-attachments/assets/cd1c03fe-cbcd-4be4-8c98-de6c04a13cc0" />


```
#Check normalised read count data
norm_expr <- calcNormFactors(expr)
plotMDS(norm_expr,xlim=c(-6,6))
```

<img width="547" height="519" alt="Screenshot 2025-10-01 at 19 24 47" src="https://github.com/user-attachments/assets/78bb9a97-8292-494f-8b75-3bbef923fba7" />


```
cpm_norm_expr <- cpm(norm_expr)
sample1 <- density(log2(cpm_norm_expr[,1]))
sample2 <- density(log2(cpm_norm_expr[,2]))
sample3 <- density(log2(cpm_norm_expr[,3]))
sample4 <- density(log2(cpm_norm_expr[,4]))
sample5 <- density(log2(cpm_norm_expr[,5]))
sample6 <- density(log2(cpm_norm_expr[,6]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Raw log2 cpm",col="red")
lines(sample2, type="l",lwd=2,col="red")
lines(sample3, type="l",lwd=2,col="red")
lines(sample4, type="l",lwd=2,col="blue")
lines(sample5, type="l",lwd=2,col="blue")
lines(sample6, type="l",lwd=2,col="blue")
```

<img width="576" height="549" alt="Screenshot 2025-10-01 at 19 25 46" src="https://github.com/user-attachments/assets/42ed8f62-8e49-454f-bbf5-e01bcfd7ad46" />


```
#Normalise and extract rpkm
expr <- DGEList(counts=data)
norm_expr <- calcNormFactors(expr)
gene_length <- read.table("gene_length.txt",stringsAsFactors=F)
head(gene_length)
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE
rpkm_norm <- rpkm(norm_expr, log=FALSE,gene.length=gene_length_vector)
write.table(rpkm_norm, file="",quote=F, sep="\t")
```

Plot the clustering of samples by gene expression information.

```
#Heat maps
library(pheatmap)
library(pvclust)

palette2 <-colorRamps::"matlab.like2"(n=200)
bootstraps = pvclust(log2(rpkm_norm+1), method.hclust="average", method.dist="euclidean")
plot(bootstraps)
```

<img width="651" height="552" alt="Screenshot 2025-10-01 at 19 29 13" src="https://github.com/user-attachments/assets/da90da67-5ab8-4459-8dc8-e724c064be07" />


```
palette2 <-colorRamps::"matlab.like2"(n=200)
pheatmap(log2(rpkm_norm+1), show_colnames=T, show_rownames=F, color = palette2, clustering_distance_cols = "euclidean", clustering_method="average") 
```

<img width="640" height="647" alt="Screenshot 2025-10-01 at 19 36 37" src="https://github.com/user-attachments/assets/296c130a-1563-4355-8d8e-82f0a3e4976d" />


Run the analysis based on the leaf gene expression, and see what differences can you observe.

### Differential gene expression

```
library("edgeR")

data <- read.table("read_counts_catkin.txt",stringsAsFactors=F,header=T, row.names=1)
head(data)
dim(data)
conditions <- factor(c("F","F","F","M","M","M"))
expr <- DGEList(counts=data,group=conditions)
norm_expr <- calcNormFactors(expr)
norm_expr

norm_expr <- estimateCommonDisp(norm_expr)
norm_expr <- estimateTagwiseDisp(norm_expr)
et <- exactTest(norm_expr)
et

p <- et$table$PValue
p_FDR <- p.adjust(p, method = c("fdr"), n = length(p))

de_results_catkin <- et$table
de_results_catkin$Padj <- p_FDR

#How many genes are significant?
sum(de_results_catkin$Padj < 0.05)
# 12753

# How many genes show 2-fold enrichment in males?
sum(de_results_catkin$Padj < 0.05 & de_results_catkin$logFC > 1)
# 6307

# How many genes show 2-fold enrichment in females?
sum(de_results_catkin$Padj < 0.05 & de_results_catkin$logFC < -1)
# 5675

# Create an MA plot
plot(de_results_catkin$logCPM, de_results_catkin$logFC, main="MA plot", xlab="log CPM", ylab="log FC")
sig <- de_results_catkin$Padj < 0.05 & (de_results_catkin$logFC > 1 | de_results_catkin$logFC < -1)
points(de_results_catkin$logCPM[sig], de_results_catkin$logFC[sig], col="orange", pch=19)
abline(h=c(-1,1), lty="dashed", col="grey")
```

<img width="615" height="639" alt="Screenshot 2025-10-01 at 23 16 23" src="https://github.com/user-attachments/assets/29199cbe-792d-4b00-95ea-b7323418b97c" />


Run the differential gene expression on the leaf samples, and see what contrasts can you make to the catkin results.


## Y gene activity decay

```
cd day4
mkdir Y_gene_activity
cd Y_gene_activity
cp /home/ubuntu/Share/day4/guppy/gene_expression/transcriptome ./
```

Obtain read counts using Salmon.....

First, we must index the transcriptome assembly. (Takes a while, so DO NOT RUN)

```
salmon index -t Poecilia_picta_transcripts.fasta -i Poecilia_picta_transcripts
```

Then, align reads to the transcriptome.

```
mkdir salmon_quantification
cd salmon_quantification
salmon quant --numBootstraps 100 --gcBias --seqBias -p 12 -l A -i ../transcriptome/Poecilia_picta_transcripts -1 ../rnaseq_reads/picta/282_1.out.fastq.gz -2 ../rnaseq_reads/picta/282_2.out.fastq.gz -o female1
```





