# Day 3 Practical

This practical will cover:

1. K-mer analyses for sex chromosome discovery
2. Estimating sequence divergence between gametologs
3. Measuring signal for Fast-x evolution
   

## 01. K-mer analyses

## 02. Gametologs divergence

First, obtain gametolog sequences and scripts.

```
mkdir ./day3/gametolog_divergence
cd ./day3/gametolog_divergence
cp ~/Server/day3/gametolog_divergence/1.gamtolog_sequences ./
cp ~/Server/day3/gametolog_divergence/scripts
```

Align sequences with **[Prank](http://wasabiapp.org/software/prank/).

```
cd ./scripts
python 01.run-prank.py ../1.gametolog_sequences
```




## 03. Fast-X evolution





* **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - A quality control tool for high throughput sequence data

First, create an output directory for FastQC. Then, run FastQC for all the paired fastq files.
