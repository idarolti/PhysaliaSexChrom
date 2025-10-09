
# Day 4 Practical - 03. Dosage compensation

## 00. Prepare work folder

A lot of the steps will be done in R, so prepare a work folder on your local machine.

```
mkdir ~Desktop/physalia/day4/dosage_compensation
cd ~Desktop/physalia/day4/dosage_compensation
```

## 01. Quantify gene expression

Obtain read counts using **[Salmon](https://combine-lab.github.io/salmon/)**. Salmon is a fast quasi-mapping approach that directly maps reads to the transcriptome without full base-to-base alignment, drastically reducing runtime and storage compared to other aligners. It does take a while to run, so you can copy the outputs directly to your working folder.

<img width="760" height="517" alt="Screenshot 2025-10-04 at 14 43 01" src="https://github.com/user-attachments/assets/166c0336-0fd1-428d-8e1b-28e39d942279" />

The two steps in Salmon are indexing the transcriptome, and then aligning with bootstrapping. 

```
salmon index -t Poecilia_picta_transcripts.fasta -i Poecilia_picta_transcripts

salmon quant --numBootstraps 100 --gcBias --seqBias -p 12 -l A -i ~/Share/day4/guppy/transcriptome/Poecilia_picta_transcripts -1 ~/Share/day4/guppy/rnaseq_reads/picta/female1_R1.fastq.gz -2 ~/Share/day4/guppy/rnaseq_reads/picta/female1_R2.fastq.gz -o female1
```

```
scp -i chrsex25.pem ubuntu@44.254.129.237:~/Share/day4/guppy/transcriptome ./
scp -i chrsex25.pem ubuntu@44.254.129.237:~/Share/day4/guppy/salmon_quantification ./
```

## 02. Obtain merged counts data

First, make list of sample names.

```
mkdir expression
cd expression

(echo "sample_id"; ls ../salmon_quantification) > samples_list.txt
```

Import counts data with **[tximport](https://www.bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)** package in R.

```
library(tximport)

# Get file with sample names
picta_samples <- read.table("samples_list.txt", header=T, sep=",")
head(picta_samples)

# Read quantification files
files <- file.path("./salmon_quantification_fullset", picta_samples$sample_id, "quant.sf")
names(files) <- paste0(picta_samples$sample_id)
all(file.exists(files))

# Run tximport
txi <- tximport(files, type = "salmon", , txOut = TRUE, countsFromAbundance = "no")
head(txi$counts)

# Write read counts
df <- data.frame(gene = rownames(txi$counts), txi$counts, check.names = FALSE)
write.table(df, file="./expression/Poecilia_picta_counts.txt", quote = FALSE, sep = ",", row.names = FALSE)
```

## 03. Convert raw counts to RPKM

First, obtain gene lengths.

```
awk '/^>/{if(seqlen){print seqname"\t"seqlen}; seqname=substr($0,2); seqlen=0; next}
     {seqlen += length($0)}
     END{print seqname"\t"seqlen}' ../transcriptome/Poecilia_picta_transcripts.fasta > ../transcriptome/gene_lengths.txt
```

Covert read counts to RPKM.

```
library("edgeR")

data <- read.table("Poecilia_picta_counts.txt",stringsAsFactors=F, header=T, row.names=1, sep=",")
head(data)
dim(data)

#Extract RPKM
expr <- DGEList(counts=data)
gene_length <- read.table("gene_lengths.txt",stringsAsFactors=F)
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

rpkm <- rpkm(expr, log=FALSE, gene.length=gene_length_vector)
write.table(rpkm, file="Poecilia_picta_rpkm.txt",quote=F, sep=",")
```

## 04. Remove lowly-expressed genes

Remove genes that do not have at least 2 RPKM in half of the individuals of one sex.

```
awk -F, '
NR==1 { print; next }
{
  female_count = 0
  male_count = 0
  for (i=2; i<=4; i++) if ($i > 2) female_count++
  for (i=5; i<=7; i++) if ($i > 2) male_count++
  if (female_count >= 2 || male_count >= 2) print
}' Poecilia_picta_rpkm.txt > Poecilia_picta_rpkm_filtered.txt

awk -F, 'NR==FNR {genes[$1]; next} FNR==1 || $1 in genes' Poecilia_picta_rpkm_filtered.txt Poecilia_picta_counts.txt > Poecilia_picta_counts_filtered.txt
```

## 05. Normalize gene expression

```
library("edgeR")

data <- read.table("Poecilia_picta_counts_filtered.txt",stringsAsFactors=F,header=T, row.names=1,sep=",")
head(data)
dim(data)
conditions <- factor(c("F","F","F","M","M","M"))

#Check raw read count data
expr <- DGEList(counts=data,group=conditions)
plotMDS(expr,xlim=c(-6,6))
```

<img width="618" height="580" alt="Screenshot 2025-10-04 at 16 46 47" src="https://github.com/user-attachments/assets/d5994cde-4c8b-4096-b991-946283323db9" />

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

<img width="661" height="601" alt="Screenshot 2025-10-04 at 16 43 21" src="https://github.com/user-attachments/assets/1144aee0-796f-4b08-9999-a8a7f800d463" />

```
#Check normalised read count data
expr <- DGEList(counts=data,group=conditions)
norm_expr <- calcNormFactors(expr)
plotMDS(norm_expr,xlim=c(-2,2))
```

<img width="666" height="578" alt="Screenshot 2025-10-04 at 16 43 55" src="https://github.com/user-attachments/assets/a55d697d-df5b-4c7f-9e27-9630544f981b" />

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

<img width="665" height="602" alt="Screenshot 2025-10-04 at 16 44 18" src="https://github.com/user-attachments/assets/9305286d-290b-40a7-9bbd-6709b3665f28" />

```
#Normalise and extract rpkm
expr <- DGEList(counts=data)
norm_expr <- calcNormFactors(expr)
gene_length <- read.table("./transcriptome/gene_lengths.txt",stringsAsFactors=F)
head(gene_length)
dim(gene_length)
expressed_genes <- rownames(data)
length(expressed_genes)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(rownames(expr),gene_length$V1),]
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(expr))
#should print TRUE

rpkm_norm <- rpkm(norm_expr, log=FALSE, gene.length=gene_length_vector)

rpkm_df <- as.data.frame(rpkm_norm)
rpkm_df$gene <- rownames(rpkm_df)
rpkm_df <- rpkm_df[, c(ncol(rpkm_df), 1:(ncol(rpkm_df)-1))]

write.table(rpkm_df, file="./expression/Poecilia_picta_rpkm_normalized.txt", quote = FALSE, sep = ",", row.names = FALSE)
```

## 06. Obtain chromosomal information for each gene

```
awk '$3 == "CDS" { 
  split($9, a, ";"); 
  for(i in a) { 
    if(a[i] ~ /^ID=/) { 
      split(a[i], b, "="); 
      split(b[2], c, ":"); 
      print c[1], $1 
    } 
  } 
}' ../transcriptome/Poecilia_picta_annotation.gff3 > ../transcriptome/gene_chromosome.txt
```

## 07. Compare expression between males and females

```
library(ggplot2) 
library(tidyr)
library(dplyr)

# Read expression matrix
expr <- read.table("Poecilia_picta_rpkm_normalized.txt", stringsAsFactors = FALSE,sep=",",header=T)
head(expr)

# Read gene-to-chromosome mapping (no header assumed)
genes_chr <- read.table("gene_chromosome.txt", stringsAsFactors = FALSE, col.names = c("gene", "chromosome"))
head(genes_chr)

# Merge by gene
merged <- merge(expr, genes_chr, by = "gene", all.x = TRUE)
head(merged)

# Apply log2 transformation (+1 pseudocount) to expression columns
expr_cols <- setdiff(colnames(merged), c("gene", "chromosome"))
merged[expr_cols] <- log2(merged[expr_cols] + 1)

# Filter rows for sex chromosome
chr12_data <- merged %>% filter(chromosome == "LG12")

# Select relevant columns and reshape to long format
long_data <- chr12_data %>%
  select(gene, female1, female2, female3, male1, male2, male3) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

# Add a sex column based on sample name prefix
long_data$sex <- ifelse(grepl("^female", long_data$sample), "Female", "Male")

# Create boxplot of expression by sex
ggplot(long_data, aes(x = sex, y = expression,fill=sex)) + 
  geom_boxplot(outlier.shape = NA,notch = TRUE) +
  scale_fill_manual(values = c("Female" = "firebrick3", "Male" = "dodgerblue")) +
  labs(title = "Sex chromosome",
       x = "Sex", y = "Log2 RPKM") +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 8))
```

<img width="404" height="531" alt="Screenshot 2025-10-04 at 16 56 45" src="https://github.com/user-attachments/assets/32de4ec2-a51d-40a2-99f9-440ad7d2b38f" />

```
# Create boxplot of expression by sample
ggplot(long_data, aes(x = sample, y = expression,fill=sex)) + 
  geom_boxplot(outlier.shape = NA,notch = TRUE) +
  scale_fill_manual(values = c("Female" = "firebrick3", "Male" = "dodgerblue")) +
  labs(title = "Sex chromosome",
       x = "Sample", y = "log2 RPKM") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img width="488" height="525" alt="Screenshot 2025-10-04 at 17 04 30" src="https://github.com/user-attachments/assets/1497c2fb-8398-4404-b480-14a8c4da798a" />

**Try plotting male and female expression for autosomal data**

autosomal_data <- merged %>% filter(chromosome != "LG12")

## 08. Compare expression between the autosomes and the sex chromosomes

```
sexchr_data <- chr12_data %>%
  select(gene, female1, female2, female3, male1, male2, male3) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")

sexchr_data$sex <- ifelse(grepl("^female", sexchr_data$sample), "Female", "Male")
sexchr_data$chr_type <- "SexChr"
long_data$chr_type <- "Autosome"
combined_data <- bind_rows(sexchr_data, long_data)

combined_males <- combined_data %>% filter(sex == "Male")

ggplot(combined_males, aes(x = chr_type, y = expression, fill = chr_type)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = c("SexChr" = "forestgreen", "Autosome" = "gray")) +
  labs(title = "Males", x = "Category", y = "log2 RPKM") +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 8))
```

<img width="352" height="531" alt="auto_sexchr_males" src="https://github.com/user-attachments/assets/a21bcd00-5b6f-4e78-801c-cb8effa56589" />

### Try plotting the same for the female dataset

## 09. Moving average plot of male and female expression

```
library(dplyr)
library(ggplot2)
library(zoo)

# Read your expression data csv
expr <- read.csv("Poecilia_picta_rpkm_normalized.txt")

# Remove the "-RA" suffix in the gene name
expr$gene <- sub("-RA$", "", expr$gene)

# Log-transform all RPKM columns
expr[, -which(names(expr) == "gene")] <- log2(expr[, -which(names(expr) == "gene")] + 1)

# Read GFF
gff <- read.table("../transcriptome/Poecilia_picta_annotation.gff3", sep="\t", header=FALSE, stringsAsFactors=FALSE)

# Assign meaningful colnames to GFF columns (based on GFF3 specification)
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
head(gff)

# Filter for chromosome LG12 (the sex chromosome)
lg12_genes <- gff %>%
  filter(seqname == "LG12" & feature == "gene")

# Extract gene ID from `attribute` column, assuming format 'ID=GENEID;...'
lg12_genes$gene <- sub(".*ID=([^;]+);.*", "\\1", lg12_genes$attribute)

# Join expression data with gene positions; only keep genes present in expression data
expr_lg12 <- expr %>%
  filter(gene %in% lg12_genes$gene) %>%
  left_join(lg12_genes[, c("gene", "start")], by = "gene")

# Order genes by start position on LG12
expr_lg12 <- expr_lg12 %>%
  arrange(start)

# Calculate average expression per sex per gene
expr_lg12 <- expr_lg12 %>%
  mutate(female_avg = rowMeans(select(., starts_with("female"))),
         male_avg = rowMeans(select(., starts_with("male"))))

# Calculate moving average (window size = e.g. 10 genes)
expr_lg12$female_ma <- zoo::rollmean(expr_lg12$female_avg, k=20, fill=NA)
expr_lg12$male_ma <- zoo::rollmean(expr_lg12$male_avg, k=20, fill=NA)

# Plot moving averages
ggplot(expr_lg12, aes(x=start/1000000)) +
  geom_line(aes(y=female_ma, color="Female")) +
  geom_line(aes(y=male_ma, color="Male")) +
  labs(title="Gene expression moving average on sex chromosome",
       x="Genomic position (bp)",
       y="Log2 RPKM") +
  scale_color_manual(values=c("Female"="firebrick", "Male"="dodgerblue")) +
  theme_classic()+
  coord_cartesian(ylim = c(0, 8))
```

<img width="752" height="436" alt="Screenshot 2025-10-08 at 09 53 02" src="https://github.com/user-attachments/assets/f354ebec-2dfb-4a91-b177-0b0b365aace9" />

## 10. Get sex-biased genes

```
library("edgeR")
library("ggplot2")

data <- read.csv("Poecilia_picta_counts_filtered.txt",row.names=1)
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

de_results <- et$table
de_results$Padj <- p_FDR

#How many genes are significant?
sum(de_results$Padj < 0.05)
# 3131

# How many genes show 2-fold enrichment in males?
sum(de_results$Padj < 0.05 & de_results$logFC > 1)
# 997

# How many genes show 2-fold enrichment in females?
sum(de_results$Padj < 0.05 & de_results$logFC < -1)
# 1147

# Volcano Plot
de_results$negLogFDR <- -log10(de_results$Padj)

ggplot(de_results, aes(x=logFC, y=negLogFDR)) +
	geom_point(alpha=0.4, size=1) +
	theme_minimal() +
	xlab("log2 Fold Change") +
	ylab("-log10 Adjusted P-value") +
	ggtitle("Volcano Plot") +
	geom_hline(yintercept = -log10(0.05), col="red", linetype="dashed") +
	geom_vline(xintercept = c(-1,1), col="blue", linetype = "dashed") +
	geom_point(data=subset(de_results, Padj < 0.05 & abs(logFC) > 1),
		aes(x=logFC, y=negLogFDR), color="orange", size=1.5)

# Output files
de_results$gene <- rownames(de_results)
de_results <- de_results[, c("gene", setdiff(colnames(de_results), "gene"))]
write.table(de_results, file = "all_genes.txt", row.names = FALSE, quote=FALSE)

# Male-biased genes: significant, logFC > 1
male_biased <- subset(de_results, Padj < 0.05 & logFC > 1)
male_biased$gene <- rownames(male_biased)
male_biased <- male_biased[, c("gene", setdiff(colnames(male_biased), "gene"))]
write.table(male_biased, file = "male_biased_genes.txt", row.names = FALSE, quote=FALSE)

# Female-biased genes: significant, logFC < -1
female_biased <- subset(de_results, Padj < 0.05 & logFC < -1)
female_biased$gene <- rownames(female_biased)
female_biased <- female_biased[, c("gene", setdiff(colnames(female_biased), "gene"))]
write.table(female_biased, file = "female_biased_genes.txt", row.names = FALSE, quote=FALSE)
```

## 11. Check for sexualization of gene expression on the sex chromosome

```
# Read GFF file, skipping comment lines
gff <- read.delim("../transcriptome/Poecilia_picta_annotation.gff3", comment.char="#", header=FALSE, stringsAsFactors=FALSE)
colnames(gff) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "phase", "attribute")

# Filter for gene features
genes_gff <- gff[gff$feature == "gene", ]

# Extract gene IDs using base R regex
genes_gff$gene_id <- sub(".*ID=([^;]+);.*", "\\1", genes_gff$attribute)

# Get unique gene-to-chromosome mapping
gene_chr <- unique(genes_gff[, c("gene_id", "chromosome")])

# Read male-biased genes
colnames(male_biased)[1] <- "gene_id"
# Remove the "-RA" suffix in the gene name
male_biased$gene_id <- sub("-RA$", "", male_biased$gene_id)

# Read female-biased genes
colnames(female_biased)[1] <- "gene_id"
# Remove the "-RA" suffix in the gene name
female_biased$gene_id <- sub("-RA$", "", female_biased$gene_id)

# Merge male-biased genes with chromosome info
male_biased_chr <- merge(male_biased, gene_chr, by="gene_id")

# Merge female-biased genes with chromosome info
female_biased_chr <- merge(female_biased, gene_chr, by="gene_id")

# Calculate total genes per chromosome
total_counts <- as.data.frame(table(gene_chr$chromosome))
colnames(total_counts) <- c("chromosome", "total_genes")

# Calculate male-biased genes per chromosome
male_counts <- as.data.frame(table(male_biased_chr$chromosome))
colnames(male_counts) <- c("chromosome", "male_biased_genes")

# Calculate female-biased genes per chromosome
female_counts <- as.data.frame(table(female_biased_chr$chromosome))
colnames(female_counts) <- c("chromosome", "female_biased_genes")

#FOR MALES

# Merge counts and replace NA with 0
chr_summary <- merge(total_counts, male_counts, by="chromosome", all.x=TRUE)
chr_summary$male_biased_genes[is.na(chr_summary$male_biased_genes)] <- 0

# Calculate proportion
chr_summary$proportion_male_biased <- chr_summary$male_biased_genes / chr_summary$total_genes

print(chr_summary)

# Filter chromosomes starting with "LG"
lg_rows <- grep("^LG", chr_summary$chromosome)
chr_summary_lg <- chr_summary[lg_rows, ]

# Extract LG number as integer
lg_nums <- as.integer(sub("LG", "", chr_summary_lg$chromosome))

# Order data frame by LG number
ordered_indices <- order(lg_nums)
chr_summary_lg <- chr_summary_lg[ordered_indices, ]

# Create a barplot with ordered LG chromosomes on x-axis
barplot(height = chr_summary_lg$proportion_male_biased,
        names.arg = chr_summary_lg$chromosome,
        las = 2, # vertical labels
        col = "cornflowerblue",
        ylab = "Proportion of Male-Biased Genes",
        main = "Proportion of Male-Biased Genes per LG Chromosome",
        cex.names = 0.8)


#FOR FEMALES

# Merge counts and replace NA with 0
chr_summary <- merge(total_counts, female_counts, by="chromosome", all.x=TRUE)
chr_summary$female_biased_genes[is.na(chr_summary$female_biased_genes)] <- 0

# Calculate proportion
chr_summary$proportion_female_biased <- chr_summary$female_biased_genes / chr_summary$total_genes

print(chr_summary)

# Filter chromosomes starting with "LG"
lg_rows <- grep("^LG", chr_summary$chromosome)
chr_summary_lg <- chr_summary[lg_rows, ]

# Extract LG number as integer
lg_nums <- as.integer(sub("LG", "", chr_summary_lg$chromosome))

# Order data frame by LG number
ordered_indices <- order(lg_nums)
chr_summary_lg <- chr_summary_lg[ordered_indices, ]

# Create a barplot with ordered LG chromosomes on x-axis
barplot(height = chr_summary_lg$proportion_female_biased,
        names.arg = chr_summary_lg$chromosome,
        las = 2, # vertical labels
        col = "cornflowerblue",
        ylab = "Proportion of Female-Biased Genes",
        main = "Proportion of Female-Biased Genes per LG Chromosome",
        cex.names = 0.8)
```
