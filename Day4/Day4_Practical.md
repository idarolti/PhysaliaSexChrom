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
mkdir genome
cd genome
cp /home/ubuntu/Share/day4/willow/genome/genome_assembly_1k.fa
hisat2-build -f genome_assembly_1k.fa genome_assembly_1k
```

Align paired-end reads to the genome, then sort.

```
cd ..
cp /home/ubuntu/Share/day4/willow/reads/ ./
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

stringtie ../hisat/female1_catkin_coordsorted.bam -o female1_catkin.gtf -A female1_catkin.gene_abund
```

Run StringTie for all samples

```
# Directory containing BAM files
bam_dir="../hisat"

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

### Obtain read counts (DO NOT RUN)

We will use **[HTSeq]([http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](https://htseq.readthedocs.io/en/release_0.11.1/count.html))**. For faster processing, bam files must be sorted by read name, so paired reads appear adjacent in the file and are counted more quickly (default option in HTSeq). If your BAM is coordinate sorted (sorted by genomic position, the default with samtools sort), you must specify -r pos to tell htseq-count the BAM is position sorted. However, this requires more memory, as mate pairs may be separated and htseq-count buffers them until it sees the mate.

```
cd ../
mkdir htseq
cd htseq

samtools sort -n -T namesorted -O bam -o ../hisat/female1_catkin_namesorted.bam female1_catkin_coordsorted.bam
htseq-count -f bam -r name -s no female1_catkin_namesorted.bam ../stringtie/merged.gtf > female1_catkin_htseqcount.txt
```












