# Day 4 Practical - 01. Differential gene expression analysis

For this part of the practical, we will use an example for quantifying gene expression in the absence of a GTF file or transcriptome reference. Several tools support genome-guided quantification, where RNA-seq reads are aligned to the genome and expression is estimated from assembled transcript structures. Here, we will focus on a pipeline based on HISAT2 and StringTie, and the model system is the willow (Salix viminalis), a dioecious species for which we have male and female RNA-seq data from both reproductive tissue (catkin) and somatic tissue (leaf).

## 00. Prepare work folder for day 4

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@44.251.209.2
mkdir day4
cd day4
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

## 01. Map RNA-seq reads

Mapping reads can be done with **[HISAT2](https://daehwankimlab.github.io/hisat2/)**, a fast and sensitive alignment program for mapping next-generation sequencing reads. 

Generate genome index.

```
hisat2-build -f genome_assembly.fa genome_assembly
```

```
mkdir differential_gene_expression
cd differential_gene_expression
```

Align paired-end reads to the genome, then sort.

```
cp -r ~/Share/day4/willow/reads/ ./
mkdir hisat
cd hisat

hisat2 ~/Share/day4/willow/genome/genome_assembly_1k -1 ../reads/female1_catkin_R1.fastq -2 ../reads/female1_catkin_R2.fastq -q --no-discordant --no-mixed --no-unal --dta -S female1_catkin.sam

samtools view -bS female1_catkin.sam | samtools sort -o female1_catkin_coordsorted.bam
rm female1_catkin.sam
```

HISAT2 options:

--no-discordant: suppress discordant alignments for paired reads

--no-mixed: suppress unpaired alignments for paired reads

--no-unal: suppresses the output of SAM records for reads that fail to align

--dta/--downstream-transcriptome-assembly: report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

Run hisat2 and sorting for all the samples.

```
reads_dir="../reads"
# Genome index prefix
genome_index="~/Share/day4/willow/genome/genome_assembly_1k"

# Loop over all fastq files
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

If hisat2 is taking too long to run, then copy the outputs to your folder

```
cd ~/day4/differential_gene_expression
rm -r hisat
cp -r ~/Share/test_day4/differential_gene_expression/hisat ./
```

## 02. Extract gene coordinates

**[StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)** - A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

```
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

## 03. Obtain read counts 

We will use **[HTSeq]([http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](https://htseq.readthedocs.io/en/release_0.11.1/count.html))**. For faster processing, bam files must be sorted by read name, so paired reads appear adjacent in the file and are counted more quickly (default option in HTSeq). If your BAM is coordinate sorted (sorted by genomic position, the default with samtools sort), you must specify -r pos to tell htseq-count the BAM is position sorted. However, this requires more memory, as mate pairs may be separated and htseq-count buffers them until it sees the mate.

**(DO NOT RUN)**

```
samtools sort -n -T namesorted -O bam -o ../hisat/female1_catkin_namesorted.bam female1_catkin_coordsorted.bam

htseq-count -f bam -r name -s no female1_catkin_namesorted.bam ../stringtie/merged.gtf > female1_catkin_htseqcount.txt
```

Copy the htseq-count output on the full dataset and copy scripts for downstream analyses.

```
cd ~/day4/differential_gene_expression
cp -r ~/Share/day4/willow/scripts ./
mkdir htseq
cp -r ~/Share/day4/willow/htseq/catkin ./htseq/
cp -r ~/Share/day4/willow/htseq/leaf ./htseq/
```

Merge read count outputs per tissue.

```
cd scripts
python3 extract-counts.py ../htseq/catkin ../htseq/catkin/read_counts_catkin.txt
head ../htseq/catkin/read_counts_catkin.txt

python3 extract-counts.py ../htseq/leaf ../htseq/leaf/read_counts_leaf.txt
head ../htseq/leaf/read_counts_leaf.txt
```

**[Here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)** is a very useful resource on the different normalization methods and the factors to consider during this process.

## 04. Filter lowly expressed genes

Copy the merged.gtf file based on the full dataset, then create a file with the length of each gene. You only need to run this using one tissue, as all genes are included in both.

```
cd ../stringtie
mkdir fullset
cp ~/Share/day4/willow/stringtie_gtfs/fullset/merged.gtf ./fullset/
cd ../scripts
python3 extract-gene-lengths.py ../htseq/catkin/read_counts_catkin.txt ../stringtie/fullset/merged.gtf ../stringtie/fullset/gene_length.txt
head ../stringtie/fullset/gene_length.txt
```

Download the read counts file for each tissue and the gene_length txt to your local directory. First, use edgeR to convert counts to RPKM.

```
library("edgeR")

data <- read.table("read_counts_catkin.txt",stringsAsFactors=F,header=T, row.names=1)
head(data)
names(data)
dim(data)

#Extract RPKM
expr <- DGEList(counts=data)
gene_length <- read.table("gene_length.txt",stringsAsFactors=F)
head(gene_length)
names(gene_length) <- c("gene","length")
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, gene %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$gene),]
gene_length_vector <- c(gene_length$length)
all(gene_length$V1 == rownames(expr))
#should print TRUE
rpkm <- rpkm(expr, log=FALSE,gene.length=gene_length_vector)
write.table(rpkm, file="rpkm_catkin.txt",quote=F, sep="\t")
```

Filter genes that do not have a minimum of 2 RPKM expression in at least half of the individuals of one sex. 

```
awk -F"\t" '
NR==1 { print; next }
{
  female_count = 0
  male_count = 0
  for (i=2; i<=4; i++) if ($i > 2) female_count++
  for (i=5; i<=7; i++) if ($i > 2) male_count++
  if (female_count >= 2 || male_count >= 2) print
}' rpkm_catkin.txt > rpkm_catkin_filtered.txt

awk -F"\t" 'NR==FNR {genes[$1]; next} FNR==1 || $1 in genes' rpkm_catkin_filtered.txt read_counts_catkin.txt > read_counts_catkin_filtered.txt
```

## 05. Normalization of gene expression with edgeR

```
library("edgeR")

data <- read.table("read_counts_catkin_filtered.txt",stringsAsFactors=F,header=T, row.names=1)
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
# Check normalised read count data
# calcNormFactors computes normalization factors, 
# accounting for differences in sequencing depth and RNA composition between samples
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
#Normalise and extract rpkm (with the calcnormfactors() function that uses edgeR's TMM method)
expr <- DGEList(counts=data)
norm_expr <- calcNormFactors(expr)
gene_length <- read.table("gene_length.txt",stringsAsFactors=F)
head(gene_length)
names(gene_length) <- c("gene","length")
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, gene %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$gene),]
gene_length_vector <- c(gene_length$length)
all(gene_length$V1 == rownames(expr))
#should print TRUE
rpkm_norm <- rpkm(norm_expr, log=FALSE,gene.length=gene_length_vector)
write.table(rpkm_norm, file="rpkm_catkin_normalized.txt",quote=F, sep="\t")
```

Plot the clustering of samples by gene expression information.

```
#Heat maps
library(pheatmap)
library(pvclust)

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

## 06. Differential gene expression (with edgeR)

```
library("edgeR")
library("ggplot2")

data <- read.table("read_counts_catkin_filtered.txt",stringsAsFactors=F,header=T, row.names=1)
head(data)
dim(data)
conditions <- factor(c("F","F","F","M","M","M"))
expr <- DGEList(counts=data,group=conditions)
norm_expr <- calcNormFactors(expr)
norm_expr

# Estimates a single common dispersion parameter across all genes. This dispersion parameter represents the biological variability across samples
norm_expr <- estimateCommonDisp(norm_expr)
# Estimates dispersions for each gene (tag) individually, allowing for gene-level variability
norm_expr <- estimateTagwiseDisp(norm_expr)
# Performs the differential expression exact test using the modeled dispersions to find genes with significant expression changes.
et <- exactTest(norm_expr)
et

# Correct for multiple testing
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

# Volcano Plot (apply -log10 for better visualization)
de_results_catkin$negLogFDR <- -log10(de_results_catkin$Padj)

ggplot(de_results_catkin, aes(x=logFC, y=negLogFDR)) +
	geom_point(alpha=0.4, size=1) +
	theme_minimal() +
	xlab("log2 Fold Change") +
	ylab("-log10 Adjusted P-value") +
	ggtitle("Volcano Plot") +
	geom_hline(yintercept = -log10(0.05), col="red", linetype="dashed") +
	geom_vline(xintercept = c(-1,1), col="blue", linetype = "dashed") +
	geom_point(data=subset(de_results_catkin, Padj < 0.05 & abs(logFC) > 1),aes(x=logFC, y=negLogFDR), color="orange", size=1.5)

# Save sex-biased genes
sex_biased <- subset(de_results_catkin, Padj < 0.05 & abs(logFC) > 1)
write.table(sex_biased, file = "sex_biased_genes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
```

<img width="676" height="665" alt="Screenshot 2025-10-08 at 00 08 52" src="https://github.com/user-attachments/assets/56bb457f-53a2-4d04-887c-9c95d105cb13" />

## 07. Differential gene expression (with DESeq2)

```
# install DESeq2 if not there already
BiocManager::install("DESeq2")

library("DESeq2")
library("ggplot2")

# Load your count data matrix (genes x samples)
data <- read.table("read_counts_catkin_filtered.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# Define sample conditions (female 'F' and male 'M')
conditions <- factor(c("F", "F", "F", "M", "M", "M"))

# Create a metadata dataframe for DESeq2
colData <- data.frame(condition = conditions)
rownames(colData) <- colnames(data)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)

# Filter low count genes (we already filtered so can skip)
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq normalization and differential expression
dds <- DESeq(dds)

# Extract results for the condition comparison (default: second level vs first level)
res <- results(dds, contrast=c("condition","M","F"))

# Add an adjusted p-value cutoff to identify significant genes
res$Padj <- p.adjust(res$pvalue, method = "fdr")

# How many genes are significant at FDR < 0.05
sum(res$Padj < 0.05, na.rm=TRUE)
15556 

# How many genes have 2-fold or greater enrichment in males (log2FC > 1)
sum(res$Padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE)
6766

# How many genes have 2-fold or greater enrichment in females (log2FC < -1)
sum(res$Padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE)
5716

# Volcano plot preparation
res$negLogFDR <- -log10(res$Padj)
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0  # Replace NA logFC with 0 for plotting

ggplot(as.data.frame(res), aes(x = log2FoldChange, y = negLogFDR)) +
  geom_point(alpha = 0.4, size = 1) +
  theme_minimal() +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  ggtitle("Volcano Plot DESeq2") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed") +
  geom_point(data = subset(res, Padj < 0.05 & abs(log2FoldChange) > 1),
             aes(x = log2FoldChange, y = negLogFDR), color = "orange", size = 1.5)
```

## 08. Plot male and female expression for genes with different sex-bias thresholds

```
library(ggplot2)
library(dplyr)
library(tidyr)

# Load expression data (genes x samples) and sex-biased gene list (with log2 fold change)
expr <- read.table("rpkm_catkin_filtered.txt", header=TRUE, row.names=1)
sexbias <- read.table("sex_biased_genes.txt", header=TRUE)

# Filter genes present in expression data
sexbias$gene_id <- rownames(sexbias)
sexbias <- sexbias %>% filter(gene_id %in% rownames(expr))

# Separate male and female samples
female_samples <- grep("^female", colnames(expr), value=TRUE)
male_samples <- grep("^male", colnames(expr), value=TRUE)

# Calculate mean expression per gene for male and female samples, plus log2 transform (adding small pseudocount)
expr$mean_female <- rowMeans(expr[,female_samples])
expr$mean_male <- rowMeans(expr[,male_samples])

# Log2 transform with pseudocount 1
expr$log2_female <- log2(expr$mean_female + 1)
expr$log2_male <- log2(expr$mean_male + 1)

# Merge fold change info
expr$gene_id <- rownames(expr)
data <- merge(expr, sexbias, by="gene_id")

# Define fold change groups (absolute value)
data$abs_logFC <- abs(data$logFC)
data$group <- cut(data$abs_logFC, breaks=c(1,3,5,7,Inf), 
                  labels=c("≥1","≥3","≥5","≥7"), right=FALSE)

# Create data for male-biased genes (log2FC > 1) and female-biased genes (log2FC < 1)
male_biased <- subset(data, logFC > 1)
female_biased <- subset(data, logFC < 1)

# Melt data for ggplot
male_long <- male_biased %>% 
  select(gene_id, group, log2_male, log2_female) %>% 
  pivot_longer(cols=c(log2_male, log2_female), 
               names_to="sex_expression", values_to="log2_RPKM") %>% 
  mutate(sex_expression = ifelse(sex_expression=="log2_male", "Male expression", "Female expression"))

female_long <- female_biased %>% 
  select(gene_id, group, log2_male, log2_female) %>% 
  pivot_longer(cols=c(log2_male, log2_female), 
               names_to="sex_expression", values_to="log2_RPKM") %>% 
  mutate(sex_expression = ifelse(sex_expression=="log2_male", "Male expression", "Female expression"))

# Plot for male-biased genes
p1 <- ggplot(male_long, aes(x=group, y=log2_RPKM, fill=sex_expression)) + 
  geom_boxplot(position=position_dodge(width=0.8)) + 
  scale_fill_manual(values=c("firebrick", "dodgerblue")) +
  labs(title="Male-biased genes", x=expression(log[2]*FC), y="Mean log2 RPKM") +
  theme_bw()

# Plot for female-biased genes
p2 <- ggplot(female_long, aes(x=group, y=log2_RPKM, fill=sex_expression)) + 
  geom_boxplot(position=position_dodge(width=0.8)) + 
  scale_fill_manual(values=c("firebrick", "dodgerblue")) +
  labs(title="Female-biased genes", x=expression(log[2]*FC), y="Mean log2 RPKM") +
  theme_bw()

# Combine plots side by side
library(gridExtra)
grid.arrange(p1, p2, ncol=2)
```

### Run the differential gene expression on the leaf samples, and see what contrasts can you make to the catkin results.
