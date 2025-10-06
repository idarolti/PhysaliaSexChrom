# Day 3 Practical - 01. K-mer analyses

This first part of the practical will cover the steps for K-mer analyses for sex chromosome discovery.

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate /home/ubuntu/miniconda3/envs/sexchr
```
Download all .sh scripts from this day's GitHub folder and upload to your day3 on the server

## 01. Setup for K-mer analyses

Within day3, set up folder for K-mer analyses

```
mkdir kmersGWAS
cd kmersGWAS
mkdir Ppicta
cd Ppicta
cp ~/Share/day3/kmersGWAS/picta/picta_kmers.table* .
cp ~/Share/day3/kmersGWAS/picta/Ppicta_phenotype.txt .
cd ..
```

## 02. Generate K-mer counts

This first step generates a list of all kmers and their presence/absence across all individuals. It can be run with the script below, but takes a very long time, so we won't run it today. Have a look at the file using the code below, and see if you can understand what it does.

```
cat ~/Share/day3/kmersGWAS/run_kmersGWAS_step1_Ppicta.sh
```

The output of this step is a file called kmers.table.table. This is a long list of every kmer found in the samples and their presence/absence in each individual. The file is already saved in the server, and we can now use this file for the later steps.  

```
cp ~/Share/day3/kmersGWAS/picta/picta_kmers.table* .
```

## 03. Generate kinship table  

Generate kinship table, in case useful for future analysis:

```
~/bin/emma_kinship_kmers -t kmers_table -k 31 --maf 0.05 > kmers_table.kinship
```

## 04. Test for association of K-mers and phenotype    
Testing for association with phenotype using the software **[PLINK](https://www.cog-genomics.org/plink/)**   
First generate PLINK compatible input file and filter for allele frequency

```
~/bin/kmers_table_to_bed -t kmers_table -k 31 -p Ppicta_phenotype.txt --maf 0.05 --mac 2 -b 1000000000 -o kmerGWAS_plink
```

Then run association test with PLINK

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

## 05. Assemble contigs from significant kmers    

Convert PLINK association file to a fasta input file for assembly of potentially sex-specific contigs

```
python plink_to_abyss_kmers.py most_significant_assoc.tab plink_abyss_input.txt
```

Assemble small contigs with ABYSS

```
ABYSS -k25 -c0 -e0 plink_abyss_input.txt -o plink_abyss_output.txt
```

## 06. Locate the assembled contigs in reference genome    

Use **[BlastN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch)**  to locate contigs in the reference genome, for this first generate a blast reference database for the reference genome and then run blast

```
REF_FASTA=~/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna

makeblastdb -in $REF_FASTA -dbtype nucl

blastn -query plink_abyss_output.txt -db $REF_FASTA -outfmt 6 -out kmers_blast.out
```

## 07. Visualize genomic locations of the assembled contigs in R

Download the outout file of the BLAST above kmers_blast.out file to your local machine   
Open **[R](https://cran.r-project.org)** or **[RStudio](https://posit.co/products/open-source/rstudio/?sid=1)**   

In R execute the code below

```
install.packages("tidyverse")
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

#####ADD EXAMPLE OF PLOTS AS IMAGE FILE###
