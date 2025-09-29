# Day 3 Practical

This practical will cover:

1. K-mer analyses for sex chromosome discovery
2. Estimating sequence divergence between gametologs
3. Measuring signal for Fast-x evolution

## 00. Prepare work folder for day 3
mkdir day3
cd day3
conda activate sexchr

## 01. K-mer analyses

## 02. Identify sex-linked sequences with **[SEX-DETector](https://pmc.ncbi.nlm.nih.gov/articles/PMC5010906/)**

## 03. Gametologs divergence

First, obtain gametolog sequences and scripts.

```
mkdir ./gametolog_divergence
cd ./gametolog_divergence
cp -r /home/ubuntu/Share/day3/gametologs_divergence/1.gametolog_sequences ./
cp -r /home/ubuntu/Share/day3/gametologs_divergence/scripts ./
```

Align sequences with **[Prank](http://wasabiapp.org/software/prank/)**.

```
cd ./scripts
python 01.run-prank.py ../1.gametolog_sequences
```

Remove gaps in alignments and short sequences.

```
python 02.remove-gaps.py ../1.gametolog_sequences ../invalid_gametologs -cutoff 300
```

Convert fasta file to **[phylip](https://www.phylo.org/index.php/help/phylip)** format, which is required by PAML. PRANK includes a built-in feature for format conversion using the -convert option along with the -f flag to specify the output format.

```
cp -r ../1.gametolog_sequences ../2.gametolog_sequences_phylip
python 03.convert-fasta-phylip.py ../2.gametolog_sequences_phylip gapsrm.fa
```

**[PAML](https://snoweye.github.io/phyclust/document/pamlDOC.pdf)** is a suite of programs for phylogenetic analyses of DNA or protein sequences using maximum likelihood (ML). The **[yn00]()** module is a method for estimating synonymous and nonsynonymous substitution rates in pairwise comparison of protein-coding DNA sequences. 

We must first create a paml control file that specifies input alignment and output files, plus options like the genetic code and analyses to perform. Then run pmal yn00.

```
python 04.make-paml-control-file.py /home/ubuntu/Share/test_day3/2.gametolog_sequences_phylip
python 05.run-paml-yn00.py /home/ubuntu/Share/test_day3/2.gametolog_sequences_phylip
```

The pairwise dS values can be found in the 2YN.dS files. We can extract the dS values for all gametologs using:

```
cd ../
mkdir plot
cp -r /home/ubuntu/Share/day3/gametologs_divergence/5.paml ./
cd 5.paml

for d in ./Gametologs_ENSORLT000000*; do
   folder=$(basename "$d")
   colname=${folder#Gametologs_}
   lastval=$(grep -E '[0-9]+\.[0-9]+' "$d/2YN.dS" | tail -n 1 | grep -Eo '[0-9]+\.[0-9]+' | tail -n 1)
   echo -e "${colname}\t${lastval}" >> ../plot/gametologs_dS.txt
done

cd ../plot
head gametologs_dS.txt
```

Merge the file with dS values with the file with positional information on the sex chromosome for each gametolog.

```
cp /home/ubuntu/Share/day3/gametologs_divergence/gametologs_position.txt ./
join -t $'\t' -1 1 -2 1 <(sort gametologs_dS.txt) <(sort gametologs_position.txt) > gametologs_dS_position.txt
head gametologs_dS_position.txt

scp -i chrsex25.pem ubuntu@44.249.25.243:/path/gametologs_dS_position.txt ~/Desktop
```

Visualize results in **[R](https://www.r-project.org/)**.

```
library(ggplot2)

gametologs <- read.csv("gametologs_dS_position.txt", row.names=1,header=F,sep="\t")

ggplot(gametologs, aes(x=V3, y=V2)) + 
	geom_point(size=2, alpha=0.7) +
	labs(title="Gametologs divergence",
			x="Chromosomal position (Mb)",
			y="Pairwise divergence dSxy") +
	theme_minimal()
```
