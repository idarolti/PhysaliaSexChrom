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
kmerblast_out <- read_tsv("kmers_blast.out") %>%
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

# Read in reference chromosome list

ref <- read_tsv("chrlength_Ppicta.tsv") %>%
  rename("scaf" = 1,
         "length" = 2)

# query' is unique for each k-mer, so we want to select the best hit for each one and plot those:
kmerblast_filter <- 
  filter(kmerblast_out, !grepl('JAVY', scaf)) %>% 
  group_by(query_sid) %>% 
  top_n(n = 1, wt = -exp_value)

# Count hits per chromosome

kmerblast_chrom <-
  kmerblast_filter %>% 
  ungroup() %>% 
  count(scaf)

kmer_chr <- left_join(ref, kmerblast_chrom) %>% arrange(scaf)

kmer_chr$n_prop <- kmer_chr$n / kmer_chr$length

(chr <-
    ggplot(kmer_chr,
           aes(x=scaf,
               y=n_prop)) + 
    geom_col() +
    theme_bw() +
    labs(title = bquote("Blasted kmerGWAS results for P. picta"),
         subtitle = "Counts per chromosome, proportional to chromosome length") +
    xlab("Chromosome") +
    ylab("ABySS contig blastn hit proportion") +
    theme_bw() + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(size = 15, vjust = .5,  angle = 90),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(size=15)
    )
)
  
### kmer positions
  
  (points <-
      filter(kmerblast_filter, scaf == "CM065354.1") %>%
      ggplot(aes(x = ref_start, y = ident.match_p)) +
      geom_point(
        size = 1) +
      labs(title = bquote("Blasted kmerGWAS results for P. picta"),
           subtitle = "All kmers")) +
      labs(x = "Position (bp)", y = "Blast similarity best hit per contig (%)") +
      theme_bw() + 
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(size = 15, vjust = .5,  angle = 90),
            axis.title.x = element_text(size = 15),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            title = element_text(size=15)
      )
```

## 02. Identify sex-linked sequences with **[SEX-DETector](https://pmc.ncbi.nlm.nih.gov/articles/PMC5010906/)**

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

dotchart(positions,labels=genes,cex=.7,main="Sex-linked genes",xlab="Chr12 start position",gcolor="black")
```

## 03. Gametologs divergence

First, obtain gametolog sequences and scripts.

```
mkdir gametologs_divergence
cd gametologs_divergence
cp -r /home/ubuntu/Share/day3/gametologs_divergence/1.gametolog_sequences ./
cp -r /home/ubuntu/Share/day3/gametologs_divergence/scripts ./
```

Align sequences with **[Prank](http://wasabiapp.org/software/prank/)**.

```
cd scripts
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
