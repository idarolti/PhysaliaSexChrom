# Day 1 Practical

## 01. Read trimming and quality check

* **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** - A read trimming tool for Illumina NGS data

    The following command will trim reads to remove adapter sequences, regions where the average Phred score in sliding windows of four bases is <15, reads for which the leading/trailing bases have a Phred score <3, and paired-end reads where either read pair is <50 bp. You can use the command on each pair (forward/R1 + reverse/R2) of fastq.gz files

    > trimmomatic PE sample_input_R1.fastq.gz sample_input_R2.fastq.gz
 
    sample_output_R1_paired.fastq.gz sample_output_R1_unpaired.fastq.gz sample_output_R2_paired.fastq.gz sample_output_R2_unpaired.fastq.gz 

    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50


## 02. Read mapping

## 03. Variant calling
