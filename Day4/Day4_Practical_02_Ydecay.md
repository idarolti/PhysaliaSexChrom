# Day 4 Practical - 02. Y gene activity decay

In this practical we will explore how to measure Y gene activity decay using estimates of allele-specific expression (ASE). ASE quanitfies the ratio of Y-linked allele expression to X-linked allele expression for genes that still have homologs on both sex chromosomes. A low expression ratio (Y/X) indicates reduced Y gene activity, signaling functional degradation or partial gene silencing on Y. We will use RNA-seq data from Poecilia picta males and females, perform genotyping to identify heterozygous sites, calculate the expression ratio for the two alleles at each site, and compare the distribution of ratios between sites on the autosomes and on the sex chromosomes, in both males and females.

## 00. Prepare work folder for day 4

```
mkdir day4
cd day4
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

## 01. Align RNA-seq reads and call SNPs

(DO NOT RUN)

The **[STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)** aligner is often recommended for mapping RNA-seq reads before calling SNPs because of its high sensitivity and accuracy in spliced alignment, especially for exon-exon junctions, and better handling of complex splicing events and repetitive regions.

```
STAR --runMode genomeGenerate --genomeDir reference_genome --genomeFastaFiles Poecilia_picta.fna --genomeSAindexNbases 12 --runThreadN 12

STAR --genomeDir reference_genome --readFilesCommand zcat --readFilesIn male1_R1.fastq.gz male1_R2.out.fastq.gz --runThreadN 12 --outFilterMultimapNmax 1
```

Sort and convert to BAM

```
samtools sort -o male1.bam male1.sam
samtools index male1.bam
```

Call SNPs

```
samtools mpileup -f Poecilia_picta.fna -b list_of_male_bams.txt -o males.mpileup
VarScan.v2.3.9.jar mpileup2snp males.mpileup --min-coverage 2 --min-ave-qual 20 --min-freq-for-hom 0.90 --p-value 1 --strand-filter 0 -—min-var-free 1e-10 --output-vcf 1 > males_unfiltered.vcf
```

## 02. Filter vcf files

Copy the outputs of mpileup2snp to your directory.

```
mkdir ase
cd ase
cp -r ~/Share/day4/ase/scripts ./
cp -r ~/Share/day4/ase/snp_calling/ ./
cd snp_calling
```

Measurements of allele-specific expression can be biased when RNA-seq reads are aligned to a single reference genome haplotype, which contains only one allele at each polymorphic site. That is because reads carrying the non-reference (alternative) allele tend to have more mismatches when aligned to the reference genome, increasing the likelihood that they are either unmapped or incorrectly mapped to a different locus. Alignment algorithms often favor sequences that match the reference perfectly or with fewer mismatches, leading to a systematic bias where reads with the reference allele are mapped more efficiently than those with alternative alleles. As a result, reads containing the non-reference allele are underrepresented, skewing allele-specific expression estimates and potentially masking true biological signals or creating artificial signals of allelic imbalance. (See this article for more details: [https://doi.org/10.1186/1471-2164-14-536](https://doi.org/10.1186/1471-2164-14-536))

To mitigate this, remove regions with a high density of SNPs.

```
cd ../scripts/
python exclude_snp_clusters.py ../snp_calling/picta_males.vcf -l 75 -m 5
```

Apply additional filters to remove triallelic SNPs and missing information in the vcf file.

```
python additional_filters.py ../snp_calling/picta_males_noclusters.vcf ../snp_calling/picta_males_noclusters_filt.vcf
```

Filter for minimum site coverage of 15 and minimum minor allele coverage of 4

```
python filter_coverage.py ../snp_calling/picta_males_noclusters_filt.vcf ../snp_calling/picta_males_noclusters_filt_coverage.vcf
```

## 03. Extract major allele ratio information

Split vcfs into autosomes and sex chromosomes.

```
mkdir ../major_allele_ratio
python split_autosomes_sexchromo.py ../snp_calling/picta_males_noclusters_filt_coverage.vcf positional_information.txt ../major_allele_ratio/picta_males_snps_autosomes.vcf ../major_allele_ratio/picta_males_snps_sexchromosomes.vcf
```

For each heterozygous site, calculate the ratio between the major and the minor alleles.

```
python average_major_allele_fraction_distribution.py ../major_allele_ratio/picta_males_snps_autosomes.vcf ../major_allele_ratio/picta_males_autosomes_major_allele_ratio.txt

python average_major_allele_fraction_distribution.py ../major_allele_ratio/picta_males_snps_sexchromosomes.vcf ../major_allele_ratio/picta_males_sexchromosomes_major_allele_ratio.txt
```

Prepare input for R plot.

```
python prepare_R_input.py ../major_allele_ratio/picta_males_autosomes_major_allele_ratio.txt ../major_allele_ratio/picta_males_sexchromosomes_major_allele_ratio.txt ../major_allele_ratio/picta_males_allele_specific_expression.txt
```

Transfer the file to your local machine, and use R to plot the major allele ratio density for the autosomes and sex chromosomes.

```
library(ggplot2)

data_m <- read.table("picta_males_allele_specific_expression.txt",stringsAsFactors=F,header=T,sep=",")

plot = ggplot(data_m, aes(x=ratio, fill=category)) + 
		geom_density(alpha=0.7) + 
		scale_fill_manual(values=c("#5F808D", "#E1B91A")) +
		coord_cartesian(xlim=c(0.5,1)) +
		theme(
           text=element_text(size=8),
           panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
           plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
           plot.title =element_text(color="black", size=14),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.text.x = element_text(size=16,color="black"),
           axis.text.y = element_text(size=16,color="black"),
           axis.title.x=element_text(size=16, color="black"), 
           axis.title.y=element_text(size=16, color="black"),
           legend.position="top",
           legend.text = element_text(colour="black",size=12),
           legend.title= element_text(colour="black",size=0)
           ) +
       	scale_x_continuous(breaks=seq(0,1,0.1)) +
       xlab('Major allele ratio') +
       ylab(expression('Density')) + 
       labs(title="P. picta males",x="Major allele ratio",y="Density")
plot
```

Run the analysis on the picta female dataset and compare results.

Run the analysis on the male and female datasets from P. reticulata and compare the results.
