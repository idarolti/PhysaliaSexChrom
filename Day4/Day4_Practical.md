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

**(DO NOT RUN)**

Mapping reads can be done with **[HISAT2](https://daehwankimlab.github.io/hisat2/)**, a fast and sensitive alignment program for mapping next-generation sequencing reads. 

```
hisat2-build -f genome_assembly_1k.fa genome_assembly_1k

hisat2 genome_assembly_1k -1 female1_catkin_R1.out.fastq -2 female1_catkin_R2.out.fastq -q --no-discordant --no-mixed --no-unal --dta -S female1_catkin.sam

samtools view -Su female1_catkin.sam | samtools sort - female1_catkin_sorted.bam
```

HISAT2 options:

--no-discordant: 

--no-mixed:

--no-unal:

--dta/--downstream-transcriptome-assembly: Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computationa and memory usage.

### Extract gene coordinates

**[StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)** - A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

```
stringtie female1_catkin_sorted.bam -o female1_catkin.gtf -p 12 -A female1_catkin.gene_abund
```

```
stringtie --merge gtfs.list -o merged.gtf
```

### Obtain read counts with **[HTSeq]([http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](https://htseq.readthedocs.io/en/release_0.11.1/count.html))**






