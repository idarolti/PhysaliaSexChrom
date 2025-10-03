# Day 4 Practical - 02. Y gene activity decay and Dosage Compensation

In this practical we will explore how to measure Y gene activity decay using estimates of allele-specific expression (ASE). ASE quanitfies the ratio of Y-linked allele expression to X-linked allele expression for genes that still have homologs on both sex chromosomes. A low expression ratio (Y/X) indicates reduced Y gene activity, signaling functional degradation or partial gene silencing on Y. We will use RNA-seq data from Poecilia picta males and females, perform genotyping to identify heterozygous sites, calculate the expression ratio for the two alleles at each site, and compare the distribution of ratios between sites on the autosomes and on the sex chromosomes, in both males and females.

## 00. Prepare work folder for day 4

```
mkdir day4
cd day4
conda activate sexchr
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
cp -r /home/ubuntu/Share/day4/guppy/ase_scripts ./
cp -r /home/ubuntu/Share/day4/guppy/snp_calling/ ./
cd snp_calling
```

Measurements of allele-specific expression can be biased when RNA-seq reads are aligned to a single reference genome haplotype, which contains only one allele at each polymorphic site. That is because reads carrying the non-reference (alternative) allele tend to have more mismatches when aligned to the reference genome, increasing the likelihood that they are either unmapped or incorrectly mapped to a different locus. Alignment algorithms often favor sequences that match the reference perfectly or with fewer mismatches, leading to a systematic bias where reads with the reference allele are mapped more efficiently than those with alternative alleles. As a result, reads containing the non-reference allele are underrepresented, skewing allele-specific expression estimates and potentially masking true biological signals or creating artificial signals of allelic imbalance. (See this article for more details: [https://doi.org/10.1186/1471-2164-14-536](https://doi.org/10.1186/1471-2164-14-536)

To mitigate this by removing regions with a high density of SNPs.

```
cd ../ase_scripts/
python exclude_snp_clusters.py ../snp_calling/males.vcf
```

Apply additional filters to remove triallelic SNPs and missing information in the vcf file.

```
python additional_filters.py ../snp_calling/males_noclusters.vcf
```

