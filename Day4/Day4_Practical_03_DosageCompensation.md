
# Day 4 Practical - 03. Dosage compensation

```
cd day4
mkdir dc
cd dc
cp -r /home/ubuntu/Share/day4/guppy/transcriptome ./
```

Obtain read counts using Salmon.....

First, we must index the transcriptome assembly. (Takes a while, so DO NOT RUN)

```
salmon index -t Poecilia_picta_transcripts.fasta -i Poecilia_picta_transcripts
```

Then, align reads to the transcriptome.

```
mkdir salmon_quantification
cd salmon_quantification

salmon quant --numBootstraps 100 --gcBias --seqBias -p 12 -l A -i ../transcriptome/Poecilia_picta_transcripts -1 /home/ubuntu/Share/day4/guppy/rnaseq_reads/picta/female1_R1.fastq.gz -2 /home/ubuntu/Share/day4/guppy/rnaseq_reads/picta/female1_R2.fastq.gz -o female1
```

This step takes a few minutes to run for each sample, so you can copy the salmon outputs to your folder.

```
cd ../
cp -r /home/ubuntu/Share/day4/guppy/salmon_quantification_fullset ./
```



