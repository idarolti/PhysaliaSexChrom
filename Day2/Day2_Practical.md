# Day 2 Practical Sex chromosome discovery - coverage and SNP based methods

This practical will cover methods on sex chromosome discovery:

1. Coverage based analysis
2. SNP-based analyses
    
We will base our analysis on the **[SexFindR pipeline](https://sexfindr.readthedocs.io/en/latest/#)** published by **[Grayson et al](https://doi.org/10.1101/2022.02.21.481346)**


## 00. Prepare work folder for day 2

```
mkdir day2
cd day2
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

## 01. Coverage based analysis
We will follow the workflow of **[DifCover](https://github.com/timnat/DifCover)** to generate a coverage file for females and males, we start with the coordinate sorted BAM files generated yesterday       
Set up directories and files

```
mkdir coverage
cd coverage
mkdir picta
mkdir reticulata
cp ~/Share/day2/Poecilia_picta_*male*_chr8_chr11_chr12.bam* picta/
cp ~/Share/day2/Poecilia_reticulata_*male*_chr8_chr11_chr12.bam* reticulata/
```

### Step 1. Create BED file

Using **[BEDtools](https://bedtools.readthedocs.io/en/latest/)** we transform the coverage information of the BAM files into a unionbed file   
This step takes a long time, so start running to check it works, then cancel (ctrl+c).

```
cd picta
~/bin/dif_cover_scripts/from_bams_to_unionbed.sh Poecilia_picta_female1_chr8_chr11_chr12.bam Poecilia_picta_male1_chr8_chr11_chr12.bam
```

### Step 2. Calculate the coverage ratio per window

In order to trun this command, copy the premade output of day2 to your working directory.
We will calculate average coverage of valid bases across all merged bed intervals for the female and the male file

```
cp ~/Share/day2/coverage/picta_sample1_sample2.unionbedcv .
```
Inspect the coverage statistics from the covstats file located in Day2/coverage to determine the parameter values for the command below   

```
cat ~/Share/day2/covstats.tab
```

Use these values to set the parameters for the following command

a = minimum coverage sample1  
A = maximum coverage sample1  
b = minimum coverage sample2  
B = maximum coverage sample2  
v = target number of valid bases in the window (set to 10000)  
l = minimum size of window to output (set to 1000)   
Then run the command below to generate the unionbed coverage file

```
~/bin/dif_cover_scripts/from_unionbed_to_ratio_per_window_CC0 -a a -A A -b b -B B -v v -l l picta_sample1_sample2.unionbedcv
```

### Step 3. Adjust coverage based on the bam ratio in covstats.tab and generate DNAcopy output file
Run the following command which will XXXXXX   
```
~/bin/dif_cover_scripts/from_ratio_per_window_to_prepare_for_DNAcopy_output.sh sample1_sample2.ratio_per_w_CC0_* bam_ratio
```

### Step 4. create plots in R

```
Rscript ~/bin/dif_cover_scripts/run_DNAcopy_from_bash.R sample1_sample2.ratio_per_w_CC0_*.log2adj_*
```
Check the output pdf file

### Step 5. filter only genomic regions with enrichment scores > p.

The script extracts from file *.DNAcopyout fragments with enrichment scores ≥ p and stores them in *.DNAcopyout{p}, (i.e. fragments where read coverage in sample1 is higher than sample2 ), and *.DNAcopyout.down{-p} fragments with enrichment scores ≤-p, (i.e. fragments where coverage in sample2 is higher than sample1 ).

Set p to 2, to filter sites with double coverage in one sample compared to the other.

```
~/bin/dif_cover_scripts/from_DNAcopyout_to_p_fragments.sh sample1_sample2.*.DNAcopyout" 2
```

### Step 6. generate histograms for further inspection

```
~/bin/dif_cover_scripts/get_DNAcopyout_with_length_of_intervals.sh *.DNAcopyout ref.length.Vk1s_sorted

echo "Generate rough histogram with given precision order"
~/bin/dif_cover_scripts/generate_DNAcopyout_len_histogram.sh *.DNAcopyout.len 1

echo "Generate histogram with bins centered at value X reporting scores from [X-0.25 to X+0.25)"
~/bin/dif_cover_scripts/generate_DNAcopyout_len_vs_scores_histogram_bin0.5.sh *.DNAcopyout.len
```

### Now run this again for P. reticulata  
### If you finish both species, try adjusting the window size in step 2, and see how this changes the results. Try smaller windows, e.g. 1,000 or larger, e.g. 50,000

## 02. SNP based analyses

Set up directories and files

```
cd day2
mkdir Fst
mkdir GWAS
mkdir SNPden
cp ~/Share/day2/SNPbased/picta_female.list Fst/
cp ~/Share/day2/SNPbased/picta_male.list Fst/
cp ~/Share/day2/SNPbased/picta_sex.list GWAS/
cp ~/Share/day2/SNPbased/picta_female.list SNPden/
cp ~/Share/day2/SNPbased/picta_male.list SNPden/
cp ~/Share/day2/vcf_files/Poecilia_picta_allchromo_merged.vcf.gz .
```

## 03. Calculate intersex Fst 
We will use **[VCFtools](https://vcftools.github.io)** to calculate Fst (fixation index, a measure of genetic differentiation between populations) between males and female in 10kb windows based on the variant file you learned how to generate yesterday   

```
cd Fst
vcftools --gzvcf ../Poecilia_picta_allchromo_merged.vcf.gz --weir-fst-pop picta_FEMALE.list --weir-fst-pop picta_MALE.list --fst-window-size 10000 --out picta.fst.10kb
cd ..
```

### Plot in R  
Download the output file picta.fst.10kb.windowed.weir.fst and the file Poecilia_picta.fna.fai and run the following commands in R. Remember to set your working directory to where the files are saved.  
Download the picta_Fst_plots.R script from day2/SNPbased and run in R.

## 04. Run association test for sex with GEMMA   

We will first generate a filtered input file that we will then further format with PLINK and then use as input for **[GEMMA](https://github.com/genetics-statistics/GEMMA)** which runs linear mixed models to test for an association between each SNP and sex.   

```
cd GWAS
vcftools --gzvcf ../Poecilia_picta_allchromo_merged.vcf.gz --plink --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --out Poecilia_picta.filt
plink --file Poecilia_picta.filt --pheno picta_sex.list --make-bed --out Poecilia_picta.plink --noweb --allow-no-sex
gemma -bfile Poecilia_picta.plink -lm 2 -o Poecilia_picta.gemma
cd output
```

Result is in the output/ folder. Download the Poecilia_picta.gemma.assoc file and the picta_GWAS_plots.R script from day2/SNPbased and run in R.

### Calculate male and female SNP density
We will format the SNP variant file with **[bcftools](https://samtools.github.io/bcftools/bcftools.html)** to subset female and male entries into separate files  
We will then calculate SNP density in 10kb windows with VCFtools
```
gunzip Poecilia_picta_allchromo_merged.vcf.gz
mv Poecilia_picta_allchromo_merged.vcf ./SNPden
cd SNPden

for SAMPLE in `cat picta_FEMALE.list`;
  do
  bcftools view -a -s ${SAMPLE} -o FEMALE_${SAMPLE}.vcf Poecilia_picta_allchromo_merged.vcf
  bgzip -c FEMALE_${SAMPLE}.vcf > FEMALE_${SAMPLE}.vcf.gz
  bcftools index FEMALE_${SAMPLE}.vcf.gz
  vcftools --gzvcf FEMALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out FEMALE_${SAMPLE}.snpden
done

for SAMPLE in `cat picta_MALE.list`;
  do
  bcftools view -a -s ${SAMPLE} -o MALE_${SAMPLE}.vcf Poecilia_picta_allchromo_merged.vcf
  bgzip -c MALE_${SAMPLE}.vcf > MALE_${SAMPLE}.vcf.gz
  bcftools index MALE_${SAMPLE}.vcf.gz
  vcftools --gzvcf MALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out MALE_${SAMPLE}.snpden
done
```

## Now do the same for Poecilia reticulata
