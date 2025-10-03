# Day 3 Practical - 02. Finding gametolog sequences

This part of the practical will cover the steps for identifying sex-linked sequences with **[SEX-DETector](https://pmc.ncbi.nlm.nih.gov/articles/PMC5010906/)**

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate sexchr
```

## 02. Genotyping

First, obtain alignment bam files and transcriptome assembly.

```
mkdir sexdetector
cd sexdetector
cp -r /home/ubuntu/Share/day3/sexdetector/transcriptome_assembly/ ./
cp -r /home/ubuntu/Share/day3/sexdetector/bam_files/ ./
```

Genotyping will be done using **[reads2snp](https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools)**. 

Make a list of bam files to run reads2snp.

```
cd bam_files
ls *.bam > bam_list.txt

(or try: ls -d "$PWD"/* > bam_list.txt)
```

Run reads2snp:

-aeb: allows alleles to have different expression levels, which is important for sex chromosome anaylses as the Y copy can be less expressed than the X copy

-min: minimum number of reads to call a genotype

-par: 0 no paraclean usage (do not clean for paralogous SNPs because X/Y SNPs can look like paralogous SNPs)

-bqt: minimum base quality

-rqt: minimum read mapping quality

```
mkdir ../reads2snp
reads2snp -aeb -min 3 -par 0 -bqt 20 -rqt 10 -bamlist bam_list.txt -bamref ../transcriptome_assembly/trinity.fasta -out ../reads2snp/reads2snp_output
```

Have a look at the two main reads2snp outputs: .alr and .gen

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

## 03. SNP segregation analysis

Identify sex-linked sequences with SEX-DETector. The software can be found **[here](https://gitlab.in2p3.fr/sex-det-family)**. For this practical we will use the original version of the software.

Copy the SEX-DETector scripts to your folder and create a folder for the output files.

```
cd day3/sexdetector/
cp /home/ubuntu/Share/day3/sexdetector/scripts/ ./
mkdir sexdetector_output
cd scripts
```

Generate gen_summary file - allows SEX-DETector to run faster. The following options should be use for the command below, -hom string (names of homogametic progeny individuals separated by commas), -het string (the names of the heterogametic progeny individuals separated by commas), -hom_par string (homogametic parent name), -het_par string (heterogametic parent name). The gen_summary file is similar to the gen file except that it shows only one occurence of each possible SNP in the dataset and shows the number of times it happens in the first column instead of the position number of the SNP.

```
./SEX-DETector_prepare_file.pl ../reads2snp/reads2snp_output.gen ../reads2snp_output.gen_summary -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father
```

Run SEX-DETector

```
./SEX-DETector.pl -alr ../reads2snp/reads2snp_output.alr -alr_gen ../reads2snp/reads2snp_output.gen -alr_gen_sum ../reads2snp/reads2snp_output.gen_summary -system xy -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father -seq -detail -detail-sex-linked -out ../sexdetector_output/Poecilia_reticulata
```

Look at the outputs produced by SEX-DETector and see which genes are inferred to be sex-linked and why.

Next, using the SEX-DETector output for the full dataset, check where on the genome do the sex-linked genes map.

```
cd day3/sexdetector
cp -r /home/ubuntu/Share/day3/sexdetector/sexdetector_output_full ./
cd sexdetector_output_full

grep "^>" Poecilia_reticulata_sex-linked_sequences.fasta | sed 's/_[^_]*$//' | sort | uniq | wc -l

awk '
  /^>/ {
    # Extract gene name by removing trailing underscore and suffix (_X, _Y, etc.)
    gene = $0
    sub(/^>/, "", gene)
    sub(/_[^_]+$/, "", gene)
    if (seen[gene]++) {
      # skip this sequence (do not print header)
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

grep ">" Poecilia_reticulata_sex-linked_sequences_unique.fasta | wc -l

blastn -db /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_reticulata -query Poecilia_reticulata_sex-linked_sequences_unique.fasta -out blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

python2 ../scripts/blast_tophit.py blastout blastout_tophits

awk -F',' '
{
  key = $2 SUBSEP $1   # chromosome + gene combined key
  count[key] = 1      # mark unique pairs
}
END {
  for (k in count) {
    split(k, parts, SUBSEP)
    chrom = parts[1]
    gene = parts[2]
    genes_per_chrom[chrom]++
  }
  for (chrom in genes_per_chrom) {
    print chrom, genes_per_chrom[chrom]
  }
}
' blastout_tophits
```

Extract the inferred sex-linked genes that align to the sex chromosome (CM002717.1).

```
awk -F',' '$2 == "CM002717.1"' blastout_tophits > blastout_tophits_sexchromo
```

Lastly, in R, plot the distribution of sex-linked genes across the sex chromosome.

```
sexlinked = read.csv("blastout_tophits_sexchromo", header=F)

dim(sexlinked)

names(sexlinked)

positions <- sexlinked$V5
genes <- sexlinked$V1

dotchart(positions,labels=genes,cex=.7,main="Sex-linked genes",xlab="Chr12 start position",xlim=c(0,26000000))
```

<img width="1069" height="628" alt="sex_linked_genes_distribution" src="https://github.com/user-attachments/assets/8c47fc86-ecfc-45c9-939a-10b8f119938c" />
