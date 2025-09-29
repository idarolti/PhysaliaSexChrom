# Day 1 Practical

This practical will cover:

1. Read quality check and trimming
2. Read mapping to genome assembly
3. Variant calling
   

## 01. Read quality check and trimming

* **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - A quality control tool for high throughput sequence data

First, create an output directory for FastQC. Then, run FastQC for all the paired fastq files.

```
mkdir fastqc_output_raw_reads

for f in /home/ubuntu/Share/day1/01.quality_trimming/raw_reads/*fastq; do fastqc $f -o ./fastqc_output_raw_reads; done
```

* **[MultiQC](https://multiqc.info)** - A tool for merging FastQC output reports of individual samples into a single summary report

This software uses as input the fastqc.zip files produced by FastQC. After running, download the .html output file to your computer to visualize the results in a web browser.

```
multiqc ./fastqc_output_raw_reads -o ./fastqc_output_raw_reads
```

* **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** - A read trimming tool for Illumina NGS data

The following command will trim reads to remove adapter sequences, regions where the average Phred score in sliding windows of four bases is <15, reads for which the leading/trailing bases have a Phred score <3, and paired-end reads where either read pair is <50 bp. You can find adapter sequences [here](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm). You can use the command on each pair (forward/R1 + reverse/R2) of fastq files.

```
mkdir trimmed_reads
    
input_dir="/home/ubuntu/Share/day1/01.quality_trimming/raw_reads/"
adapter_dir="/home/ubuntu/Share/day1/01.quality_trimming/"
output_dir=./trimmed_reads
    
trimmomatic PE \
   $input_dir/sample1_R1.fastq \
   $input_dir/sample1_R2.fastq \
   $output_dir/sample1_output_R1_paired.fastq.gz \
   $output_dir/sample1_output_R1_unpaired.fastq.gz \
   $output_dir/sample1_output_R2_paired.fastq.gz \
   $output_dir/sample1_output_R2_unpaired.fastq.gz \
   ILLUMINACLIP:$adapter_dir/adapters.fa:2:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```

For comparison, run FastQC and MultiQC on the trimmed reads

```
mkdir fastqc_output_trimmed_reads
for f in ./trimmed_reads/*_paired.fastq.gz; do fastqc $f -o ./fastqc_output_trimmed_reads; done
multiqc ./fastqc_output_trimmed_reads -o ./fastqc_output_trimmed_reads
```

Try running the trimming and fastqc commands for sample2.

## 02. Read mapping

* **[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)** - A tool for aligning short-read data to a reference genome or genomic sequences

First, build an index for the reference genome. DO NOT RUN!

```
mkdir read_alignments

bowtie2-build Poecilia_picta.fna Poecilia_picta
```

Then, align each pair of reads to the indexed genome using bowtie2 and convert the output alignment sam file into a sorted bam file.

```
bowtie2 -p12 -x /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta \
   -1 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R1_subset.fastq \
   -2 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R2_subset.fastq \
   | samtools view -b -S - | samtools sort - -o ./read_alignments/Poecilia_picta_female1_subset.bam
```

TAKES 15 MIN!
```
bowtie2 -p12 -x /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta \
   -1 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R1_chr12.fastq \
   -2 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R2_chr12.fastq \
   | samtools view -b -S - | samtools sort - -o ./read_alignments/Poecilia_picta_female1_chr12.bam
```

Try running the bowtie mapping for female2 subset.

## 03. Variant calling

* **[GATK](https://gatk.broadinstitute.org/hc/en-us)** - A genomic analysis toolkit focused on variant discovery.

Create sequence dictionary for the reference sequence. DO NOT RUN!

```
gatk CreateSequenceDictionary -R Poecilia_picta.fna -O Poecilia_picta.dict
```

Call SNPs using [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)

```
mkdir snp_calling

picard AddOrReplaceReadGroups I=./read_alignments/Poecilia_picta_female1_subset.bam O=./read_alignments/Poecilia_picta_female1_subset_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=picta_female1

samtools index ./read_alignments/Poecilia_picta_female1_subset_RG.bam

gatk HaplotypeCaller \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -I ./read_alignments/Poecilia_picta_female1_subset_RG.bam \
   -O ./snp_calling/Poecilia_picta_female1_subset.gvcf --emit-ref-confidence GVCF \
   --min-base-quality-score 30 --pcr-indel-model NONE --sample-name picta_female1
```

Perform genotyping of variants using [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs)

```
gatk GenotypeGVCFs \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   --variant /home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_female1_chr12.gvcf \
   -O ./snp_calling/Poecilia_picta_female1_chr12.genotyped.gvcf
```

Filter variants using [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants) and [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)

```
gatk SelectVariants \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./snp_calling/Poecilia_picta_female1_chr12.genotyped.gvcf \
   -O ./snp_calling/Poecilia_picta_female1_chr12.selectvar.gvcf --restrict-alleles-to BIALLELIC --select-type-to-include SNP

gatk VariantFiltration \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./snp_calling/Poecilia_picta_female1_chr12.selectvar.gvcf \
   -O ./snp_calling/Poecilia_picta_female1_chr12.selectvar_filtered.gvcf \
   --filter-expression "QUAL <= 30.0 || DP <= 20" --filter-name "low_qual_or_dp"
```

Try genotyping and variant filters on female1_chr8.
