# Day 3 Practical - 02. Finding gametolog sequences

This part of the practical will cover the steps for identifying sex-linked sequences with **[SEX-DETector](https://pmc.ncbi.nlm.nih.gov/articles/PMC5010906/)**

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

## 02. Genotyping

First, obtain alignment bam files and transcriptome assembly.
```
mkdir sexdetector
cd sexdetector
cp -r ~/Share/day3/sexdetector/transcriptome_assembly/ ./
cp -r ~/Share/day3/sexdetector/bam_files/ ./
cp -r ~/Share/day3/sexdetector/scripts/ ./
```

Genotyping will be done using **[reads2snp](https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools)**. This is preferred when using SEX-DETector with RNA-seq data because it allows for allelic expression biases (important for sex chromosome studies because Y or W alleles may be less expressed).

Make a list of bam files to run reads2snp.

```
cd bam_files
ls -d "$PWD"/* > bam_list.txt
```

Run reads2snp:

-aeb: allows alleles to have different expression levels, which is important for sex chromosome anaylses as the Y copy can be less expressed than the X copy

-min: minimum number of reads to call a genotype

-par: 0 (do not remove SNPs that appear to come from paralogous sequences, avoid overfiltering as X/Y SNPs can look like paralogous SNPs)

-bqt: minimum base quality

-rqt: minimum read mapping quality

```
cd ../scripts
mkdir ../reads2snp
./reads2snp_2.0.64.bin -aeb -min 3 -par 0 -bqt 20 -rqt 10 -bamlist ../bam_files/bam_list.txt -bamref ../transcriptome_assembly/trinity.fasta -out ../reads2snp/reads2snp_output
```

Have a look at the two main reads2snp outputs: .alr and .gen

```
cd ../reads2snp
head reads2snp_output.alr
head reads2snp_output.gen
```

The output .alr file is a tab delimited file that shows each position of each gene on a different line with the major allele (most common allele in number of reads), then M if the position is monomorphic or P if the position shows different alleles, then for each individual the total read number if the position is monomorphic or the total read number followed by the number of reads for each base [A/C/G/T] if the position is polymorphic:

```
>contig_name
	maj	M/P	individual1	individual2
	T	M	8	1
	C	P	17[0/0/0/17]	14[0/6/0/8]
```

The output .gen file is a tab delimited file showing each position of each gene on a separate line with the position number starting from 1, then for each individual the inferred genotype followed by a pipe and the posterior probability of the inferred genotype:

```
>contig_name
  position	individual1	individual2
  1	TT|1	NN|0
  2	TT|1	TT|0.98
```

## 03. SNP segregation analysis (subset)

Identify sex-linked sequences with SEX-DETector. The software can be found **[here](https://gitlab.in2p3.fr/sex-det-family)**. For this practical we will use the original version of the software. Newer versions can also identify the type of sex chromosome system (XY or ZW) if you don't know it a priori.

```
mkdir ../sexdetector_output
cd ../scripts
```

Generate gen_summary file - allows SEX-DETector to run faster. The following options should be use for the command below, -hom string (names of homogametic progeny individuals separated by commas), -het string (the names of the heterogametic progeny individuals separated by commas), -hom_par string (homogametic parent name), -het_par string (heterogametic parent name). The gen_summary file is similar to the gen file except that it shows only one occurence of each possible SNP in the dataset and shows the number of times it happens in the first column instead of the position number of the SNP.

```
./SEX-DETector_prepare_file.pl ../reads2snp/reads2snp_output.gen ../reads2snp/reads2snp_output.gen_summary -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father

head ../reads2snp/reads2snp_output.gen_summary
```

Run SEX-DETector

```
./SEX-DETector.pl

./SEX-DETector.pl -alr ../reads2snp/reads2snp_output.alr -alr_gen ../reads2snp/reads2snp_output.gen -alr_gen_sum ../reads2snp/reads2snp_output.gen_summary -system xy -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father -seq -detail -detail-sex-linked -out ../sexdetector_output/Poecilia_reticulata

cd ../sexdetector_output
```

The main outputs from SEX-DETector are:
- SNPs_detail.txt: information on each SNP of each gene (position in gene, likelihood autosomal/sex-linked, genotypes)
- sex-linked_detail.txt: information on SNP inferred as sex-linked (position, likelihood, expression of each allele, error rates)
- sex-linked_sequences.fasta: X and Y sequences fasta file
- assignment.txt: outputs the probability for each gene to be autosomal, sex-linked (with X/Y alleles), sex-linked (XO hemizygous)
  
Look at the outputs produced by SEX-DETector. How many genes are inferred to be sex-linked?


## 04. SNP segregation analysis (full dataset)

Next, get the SEX-DETector output for the full dataset.

```
cd ~/day3/sexdetector
cp -r ~/Share/day3/sexdetector/sexdetector_output_full ./
cd sexdetector_output_full
```

Find how many genes are inferred as autosomal versus sex-linked.

```
cat Poecilia_reticulata_assignment.txt | grep -w "sex-linked" -c
cat Poecilia_reticulata_assignment.txt | grep -w "autosomal" -c
```

Obtain a single sequence for each gene (to then blast onto the assembly).

```
head Poecilia_reticulata_sex-linked_sequences.fasta

awk '
  # on lines that start with ">"
  /^>/ {
	# save the header line in a variable "gene"
    gene = $0
	# remove the leading ">" symbol, leaving just the gene name
    sub(/^>/, "", gene)
	# remove the suffix after the last underscore in the gene name
    sub(/_[^_]+$/, "", gene)
	# use an array "seen: to check if this gene ID was seen before
    if (seen[gene]++) {
      skip = 1
    } else {
      # print header
      print $0
      skip = 0
    }
  }
  !/^>/ {
    # print sequence only if not skipping
    if (!skip) print $0
  }
' Poecilia_reticulata_sex-linked_sequences.fasta > Poecilia_reticulata_sex-linked_sequences_unique.fasta

grep ">" Poecilia_reticulata_sex-linked_sequences_unique.fasta -c
```

Use Blast so see where on the genome the identifyied sex-linked genes align. A Blast database (with makeblastdb) has been already created in the Shared folder.

```
blastn -db ~/Share/day1/02.read_mapping/reference_genome/Poecilia_reticulata -query Poecilia_reticulata_sex-linked_sequences_unique.fasta -out blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

head blastout
```

Identify top blast hits for each sequence. This script takes a blast output file (format: outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq") and identifies the top blast hit for each query. Top blast hit = minimum 30 pidentity, greatest blast score and greatest pidentity. If a query has two hits with identical blast score and pidentity one is chosen randomly as the tophit. 

```
python ../scripts/blast_tophit.py blastout blastout_tophits
head blastout_tophits
```

Count the number of times a chromosome is assigned a sex-linked gene.

```
awk -F',' '
{
  # chromosome + gene combined key
  key = $2 SUBSEP $1
  # mark unique pairs
  count[key] = 1
}
END {
  for (k in count) {
	# split key back into chromosome and gene
    split(k, parts, SUBSEP)
    chrom = parts[1]
    gene = parts[2]
	# increment count of unique genes per chromosome
    genes_per_chrom[chrom]++
  }
  for (chrom in genes_per_chrom) {
	# output chromosome and gene count
    print chrom, genes_per_chrom[chrom]
  }
}
' blastout_tophits
```

Extract the inferred sex-linked genes that align to the sex chromosome (CM002717.1).

```
awk -F',' '$2 == "CM002717.1"' blastout_tophits > blastout_tophits_sexchromo
```

Transfer this file to your local machine, and plot the distribution of sex-linked genes across the sex chromosome usind R.

```
sexlinked = read.csv("blastout_tophits_sexchromo", header=F)

dim(sexlinked)
head(sexlinked)

names(sexlinked) <- c("Gene","Chromosome","Bitscore","PIdentity","Start","End")
head(sexlinked)

positions <- sexlinked$Start
genes <- sexlinked$Gene

dotchart(sexlinked$Start/1000000,labels=sexlinked$Gene,cex=.7,main="Sex-linked genes",xlab="Sex chromosome position (Mb)",xlim=c(0,26))
```

Run the last part of the analysis (04.SNP segregation analysis) using the SEX-DETector output from another reticulata cross. What differences can you notice in the distribution of sex-linked genes? 

```
cd ~/day3/sexdetector/
cp -r ~/Share/day3/sexdetector/sexdetector_output_full_Ret9 ./
```
