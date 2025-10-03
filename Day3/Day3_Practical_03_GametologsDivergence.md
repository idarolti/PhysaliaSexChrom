# Day 3 Practical - 03. Gametologs divergence

This part of the practical will cover the steps for estimating sequence divergence between gametologs.

## 00. Prepare work folder for day 3

```
mkdir day3
cd day3
conda activate sexchr
```

## 01. Align gametolog sequences

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

## 02. Prepare files for PAML.

Convert fasta file to **[phylip](https://www.phylo.org/index.php/help/phylip)** format, which is required by PAML. PRANK includes a built-in feature for format conversion using the -convert option along with the -f flag to specify the output format.

```
cp -r ../1.gametolog_sequences ../2.gametolog_sequences_phylip
python 03.convert-fasta-phylip.py ../2.gametolog_sequences_phylip gapsrm.fa
```

**[PAML](https://snoweye.github.io/phyclust/document/pamlDOC.pdf)** is a suite of programs for phylogenetic analyses of DNA or protein sequences using maximum likelihood (ML). The **[yn00]()** module is a method for estimating synonymous and nonsynonymous substitution rates in pairwise comparison of protein-coding DNA sequences. 

We must first create a paml control file that specifies input alignment and output files, plus options like the genetic code and analyses to perform. Then run pmal yn00.

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

## 03. Run PAML and extract pairwise divergence estimates.

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

## 04. Plot dSxy estimates.

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

<img width="801" height="437" alt="gametologs divergence plot" src="https://github.com/user-attachments/assets/0ab6471f-418d-44ae-a060-199e4f6c44c0" />
