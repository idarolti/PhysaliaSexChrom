# Day 3 Practical - 03. Gametologs divergence

This part of the practical will cover the steps for estimating sequence divergence between gametologs.

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

## 01. Align gametolog sequences

First, obtain gametolog sequences and scripts.

```
mkdir gametologs_divergence
cd gametologs_divergence
cp -r ~/Share/day3/gametologs_divergence/1.gametolog_sequences ./
cp -r ~/Share/day3/gametologs_divergence/scripts ./
```

Simplify gene names in the fasta files.

```
root_dir="1.gametolog_sequences"

# Loop through all subfolders
for subdir in "$root_dir"/Gametologs_*; do
    # Define the expected fasta filename inside the subdir
    fasta_file=$(find "$subdir" -maxdepth 1 -type f -name "*.stripped.cdna.fa")
    
    # Skip if no such fasta file found
    if [[ -z "$fasta_file" ]]; then
        echo "No .stripped.cdna.fa file found in $subdir"
        continue
    fi
    
    echo "Processing $fasta_file ..."
    
    # Use sed to simplify headers in place
    sed -E '/^>/ s/(_.*)//' "$fasta_file" > "${fasta_file}.tmp" && mv "${fasta_file}.tmp" "$fasta_file"
done
```

Align sequences with **[Prank](http://wasabiapp.org/software/prank/)**. Aligning sequences is important because dS estimation depends on the correct placement of codons and identifying homologous nucleotide positions. This part takes a few seconds to run per gametolog pair, so we can start the command and then quit if it takes too long.

```
cd scripts
python 01.run-prank.py ../1.gametolog_sequences
rm -r ../1.gametolog_sequences
cp -r ~/Share/day3/gametologs_divergence/1.gametolog_sequences_prank/ ../
```

Remove gaps in alignments and short sequences. Gaps represent insertions or deletions (indels) that can disrupt the reading frame and homology established at the codon level. For statistical reasons, you can also filter for a minimum gene length (shorter sequences provide very few codons for substitution rate estimation, which can lead to unstable and unreliable substitution estimates due to insufficient mutation counts and sampling noise).

```
python 02.remove-gaps.py ../1.gametolog_sequences_prank ../invalid_gametologs -cutoff 300
```

## 02. Prepare files for PAML

Convert fasta file to **[phylip](https://www.phylo.org/index.php/help/phylip)** format, which is required by PAML. PRANK includes a built-in feature for format conversion using the -convert option along with the -f flag to specify the output format.

This next step will delete all the files that don't end in "gapsrm.fa" and the convert that file to phylip format. So if we want to still have a copy of the previous files, we can duplicate the 1.gametolog_sequences folder.

```
cp -r ../1.gametolog_sequences_prank ../2.gametolog_sequences_phylip
python 03.convert-fasta-phylip.py ../2.gametolog_sequences_phylip gapsrm.fa
```

**[PAML](https://snoweye.github.io/phyclust/document/pamlDOC.pdf)** is a suite of programs for phylogenetic analyses of DNA or protein sequences using maximum likelihood (ML). The **[yn00]()** module is a method for estimating synonymous and nonsynonymous substitution rates in pairwise comparison of protein-coding DNA sequences. 

We must first create a paml control file that specifies input alignment and output files, plus options like the genetic code and the analyses to be performed. Then run pmal yn00.

```
cd ../2.gametolog_sequences_phylip
for d in Gametologs_*; do
    if [ -d "$d" ]; then
        base="${d#Gametologs_}"
        ctl_path="$d/$base.ctl"
        cat > "$ctl_path" <<EOF
seqfile = $base.phy * sequence data filename
outfile = $base.txt * main result file name

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
weighting = 0
commonf3x4 = 0
ndata = 1
EOF
        echo "Wrote $ctl_path"
    fi
done
```

## 03. Run PAML and extract pairwise divergence estimates

```
for d in Gametologs_*; do
    if [ -d "$d" ]; then
        echo "Running yn00 in $d"
        (
            cd "$d"
            yn00 *.ctl
        )
    fi
done
```

The pairwise dS values can be found in the 2YN.dS files. We can extract the dS values for all gametologs using:

```
cd ../
mkdir plot

cd 2.gametolog_sequences_phylip
for d in ./Gametologs_ENSORLT000000*; do
   folder=$(basename "$d")
   colname=${folder#Gametologs_}
   lastval=$(grep -E '[0-9]+\.[0-9]+' "$d/2YN.dS" | tail -n 1 | grep -Eo '[0-9]+\.[0-9]+' | tail -n 1)
   echo -e "${colname}\t${lastval}" >> ../plot/gametologs_dS.txt
done

cd ../plot
head gametologs_dS.txt
```

## 04. Plot dSxy estimates

Merge the file with dS values with the file with positional information on the sex chromosome for each gametolog.

```
cp ~/Share/day3/gametologs_divergence/gametologs_position.txt ./
head gametologs_position.txt

join -t $'\t' -1 1 -2 1 <(sort gametologs_dS.txt) <(sort gametologs_position.txt) > gametologs_dS_position.txt
head gametologs_dS_position.txt
```

Download the gametologs_dS_position.txt file and visualize results in **[R](https://www.r-project.org/)**.

```
scp -i chrsex25.pem ubuntu@44.249.25.243:/path/gametologs_dS_position.txt ~/Desktop
```

```
library(ggplot2)

gametologs <- read.csv("gametologs_dS_position.txt", row.names=1,header=F,sep="\t")
head(gametologs)
names(gametologs) <- c("dS","Position")

ggplot(gametologs, aes(x=Position, y=dS)) + 
	geom_point(size=2, alpha=0.7) +
	labs(title="Gametologs divergence",
			x="Sex chromosome position (Mb)",
			y="Pairwise divergence dSxy") +
	theme_minimal()
```
