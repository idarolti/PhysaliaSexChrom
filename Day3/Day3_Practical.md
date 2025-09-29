# Day 3 Practical

This practical will cover:

1. K-mer analyses for sex chromosome discovery
2. Identifying sex-linked sequences
3. Estimating sequence divergence between gametologs

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate sexchr
```

## 01. K-mer analyses

Set up folder

```
mkdir kmersGWAS
cd kmersGWAS
mkdir Ppicta
cd Ppicta
cp /bin/Ppicta_dirlist.txt .
cp /bin/Ppicta_phenotype.txt .
cp /bin/Poecilia_picta* .
```

Generate kmers.table

This is a table of all kmers and their presence/absence between all individuals

```
bash run_kmersGWAS_step1_Ppicta.sh
```

Generate kinship table, in case useful for future analysis

```
bin/emma_kinship_kmers -t kmers_table -k 31 --maf 0.05 > kmers_table.kinship
```

Testing for associating with phenotype using plink

```
bin/kmers_table_to_bed -t kmers_table -k 31 -p Ppicta_phenotype.txt --maf 0.05 --mac 2 -b 1000000000 -o kmerGWAS_plink
```

Run plink

```
plink --noweb --bfile kmerGWAS_plink.0 --allow-no-sex --assoc --out kmers
```

Edit format of output file

```
awk -v OFS='\t' '{ $1=$1; print }' kmers.assoc > kmers.assoc.tab
```

Check the output file for the most significant p-value, and filter for only kmers with this value

```
cut -f 9 kmers.assoc.tab | grep -v 'P' | sort -g | head -1
```

Edit the command for the p-value output from the step above

```
awk '$9 < P' kmers.assoc.tab > most_significant_assoc.txt
```

Convert plink file

```
python plink_to_abyss_kmers.py most_significant_assoc.tab plink_abyss_input.txt
```

Assemble small contigs with ABYSS

```
ABYSS -k25 -c0 -e0 plink_abyss_input.txt -o plink_abyss_output.txt
```

Blast contigs against a reference genome

```
REF_FASTA=Poecilia_picta.fna

makeblastdb -in $REF_FASTA -dbtype nucl

blastn -query plink_abyss_output.txt -db $REF_FASTA -outfmt 6 -out kmers_blast.out
```

Download kmers_blast.out file
In R

```
library(tidyverse)

# Read in data
kmerblast_out <- read_tsv(paste0(species, "/kmers_blast.out"), col_names = F) %>%
  rename("query_sid" = 1,
         "scaf" = 2,
         "ident.match_p" = 3,
         "alignmt_len" = 4,
         "mismatch_n" = 5,
         "gapopen_n" = 6,
         "query_start" = 7,
         "query_end" = 8,
         "ref_start" = 9,
         "ref_end" = 10,
         "exp_value" = 11,
         "bitscore" = 12)

# query' is unique for each k-mer, so we want to select the best hit for each one and plot those:
kmerblast_filter <- 
  kmerblast_out %>% 
  group_by(query_sid) %>% 
  top_n(n = 1, wt = -exp_value)

# Count hits per chromosome

kmerblast_chrom <-
  kmerblast_filter %>% 
  ungroup() %>% 
  count(scaf)

(chr <-
    ggplot(kmerblast_chrom,
           aes(x=scaf,
               y=n)) + 
    geom_col()

### kmer positions

(points <-
	filter(kmerblast_filter, scaf == "") %>%
    ggplot(aes(x = ref_start, y = ident.match_p)) +
    geom_point(
      size = 1)

```

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
