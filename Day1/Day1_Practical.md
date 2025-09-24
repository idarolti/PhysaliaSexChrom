# Day 1 Practical

## 01. Read trimming and quality check

* **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** - A read trimming tool for Illumina NGS data

    The following command will trim reads to remove adapter sequences, regions where the average Phred score in sliding windows of four bases is <15, reads for which the leading/trailing bases have a Phred score <3, and paired-end reads where either read pair is <50 bp. You can use the command on each pair (forward/R1 + reverse/R2) of fastq.gz files.

    > trimmomatic PE sample_input_R1.fastq.gz sample_input_R2.fastq.gz sample_output_R1_paired.fastq.gz sample_output_R1_unpaired.fastq.gz sample_output_R2_paired.fastq.gz sample_output_R2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

* **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - A quality control tool for high throughput sequence data

    First, create an output directory for FastQC.
    > mkdir fastqc_output

    Then, run FastQC for all the paired fastq files that resulted from trimmomatic.
    > for f in _paired.fastq.gz; do fastqc $f -o ./fastqc_output; done

* **[MultiQC](https://multiqc.info)** - A tool for merging FastQC output reports of individual samples into a single summary report

    This software uses as input the fastqc.zip files produced by FastQC.
    > multiqc ./fastqc_output -o ./fastqc_output

    Download the .html output file to your computer to visualize the results in the web browser.

## 02. Read mapping

* **[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)** - A tool for alining short-read data to a reference genome or genomic sequences

    First, build an index for the reference genome.
    > bowtie2-build reference.fasta reference

    Then, align each pair of reads to the indexed genome using bowtie2 and convert the output alignment sam file in to a sorted bam file.
    > bowtie2 -x reference -1 sample_output_R1_paired.fastq.gz -2 sample_output_R2_paired.fastq.gz -p 12 | samtools view -b -S - | samtools sort - -o sample.bam

## 03. Variant calling

* **[GATK](https://gatk.broadinstitute.org/hc/en-us)** - A genomic analysis toolkit focused on variant discovery.

    Variant calling using GATK runs over 4 steps
  
1. Call SNPs and indels per individual, per chromosome using [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)

    Calling variants separately per individual and per chromosome allows these to be run in parallel, therefore faster.

```

mkdir snpcalling
cd snpcalling

# use a for loop to run over every bam file and every chromosome in the reference

for i in `ls ../mapping/*.bam`
    do
    for j in chr{1..n}
        do
        gatk HaplotypeCaller \
         -R reference \
         -I $i \
         -L $j \
         -O $i.$j.gvcf \
         --emit-ref-confidence GVCF \
         --min-base-quality-score 30 \
         --pcr-indel-model NONE

```

2. Combine all GVCFs per chromosome for multiple individuals using [CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/13832710975771-CombineGVCFs)

```

gatk CombineGVCFs \
-R reference \
--variant sample_1.chr1.gvcf \
--variant sample_2.chr1.gvcf \
--variant sample_n.chr1.gvcf \
-O chr1_combined.g.vcf.gz

# repeat this for every chromosome

```

3. Perform joint genotyping of variants for multiple samples using [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs)

```

gatk GenotypeGVCFs \
-R reference \
--variant chr1_combined.g.vcf.gz \
-O chr1_genotyped.g.vcf.gz

# repeat this for every chromosome

```

4. Merge together GVCFs for all chromosomes using [GatherVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/13832683937435-GatherVcfs-Picard)

```

gatk GatherVcfs \
-R reference \
--variant chr1_genotyped.g.vcf.gz \
--variant chr2_genotyped.g.vcf.gz \
--variant chrn_genotyped.g.vcf.gz \
-O AllSamples_Genotyped.merged.g.vcf.gz

```
