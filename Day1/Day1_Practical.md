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

for f in ./Shared/day1/01.quality_trimming/raw_reads/*fastq; do fastqc $f -o ./fastqc_output_raw_reads; done

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
    
input_dir="./Shared/day1/01.quality_trimming/raw_reads/"
adapter_dir="./Shared/day1/01.quality_trimming/"
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
for f in ./trimmed_reads/*.gz; do fastqc $f -o ./fastqc_output_trimmed_reads; done
multiqc ./fastqc_output_trimmed_reads -o ./fastqc_output_trimmed_reads

```

## 02. Read mapping

* **[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)** - A tool for aligning short-read data to a reference genome or genomic sequences

    First, build an index for the reference genome. DO NOT RUN!

```
mkdir reference_genome
mkdir read_alignments
    
cp ./Shared/day1/02.read_mapping/reference_genome/Poecilia_picta ./reference_genome
    
bowtie2-build ./reference_genome/Poecilia_picta.fna ./reference_genome/Poecilia_picta
```

Then, align each pair of reads to the indexed genome using bowtie2 and convert the output alignment sam file in to a sorted bam file.

```
bowtie2 -x ./reference_genome/Poecilia_picta -1 ./Shared/day1/02.read_mapping/reads/Poecilia_picta_female1_R1_subset.fastq -2 ./Shared/day1/02.read_mapping/reads/Poecilia_picta_female1_R2_subset.fastq -p 12 | samtools view -b -S - | samtools sort - -o ./read_alignments/Poecilia_picta_female1_subset.bam
```

## 03. Variant calling

* **[GATK](https://gatk.broadinstitute.org/hc/en-us)** - A genomic analysis toolkit focused on variant discovery.
  
1. Call SNPs using [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)

```

mkdir snpcalling
cd snpcalling

gatk HaplotypeCaller \
    -R reference \
    -I sample.bam \
    -L chr1 \
    -O sample.chr1.gvcf \
    --emit-ref-confidence GVCF \
    --min-base-quality-score 30 \
    --pcr-indel-model NONE

```

2. Perform genotyping of variants using [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs)

```

gatk GenotypeGVCFs \
-R reference \
--variant sample.chr1.gvcf \
-O sample.chr1.genotyped.g.vcf.gz

```

