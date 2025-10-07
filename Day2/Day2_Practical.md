# Day 2 Practical Sex chromosome discovery - coverage and SNP based methods

This practical will cover methods on sex chromosome discovery:

1. Coverage based analysis
2. SNP-based analyses
    
We will base our analysis on the **[SexFindR pipeline](https://sexfindr.readthedocs.io/en/latest/#)** published by **[Grayson et al](https://doi.org/10.1101/2022.02.21.481346)**


## 00. Prepare work folder for day 2
Open your terminal application as yesterday, connect to the server as yesterday, replace "25" with the number of your user account and use today's IP adress

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@44.251.209.2
mkdir day2
cd day2
conda activate /home/ubuntu/miniconda3/envs/sexchr
```
Don't forget to open a Filezilla connection, and change the IP to today's address  

## 01. Coverage based analysis
We will follow the workflow of **[DifCover](https://github.com/timnat/DifCover)** to generate a coverage file for females and males, we start with the coordinate sorted BAM files generated yesterday.       
Set up directories and files

```
mkdir coverage
cd coverage
mkdir picta
mkdir reticulata
mkdir dif_cover_scripts
cp ~/Share/day2/coverage/dif_cover_scripts/* ./dif_cover_scripts
```

### Step 1. Create BED file

Using **[BEDtools](https://bedtools.readthedocs.io/en/latest/)** we transform the coverage information of the BAM files into a unionbed file   
This step takes a long time, so start running to check it works, then cancel (ctrl+c).

```
cd picta
../dif_cover_scripts/from_bams_to_unionbed.sh ~/Share/day2/bam_files/Poecilia_picta_female1_subset.bam ~/Share/day2/bam_files/Poecilia_picta_male2_subset.bam
```

### Step 2. Calculate the coverage ratio per window

In order to run this command, use the premade output of day2.
We will calculate average coverage of valid bases across all merged bed intervals for the female and the male file

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
from_unionbed_to_ratio_per_window_CC0 -a a -A A -b b -B B -v v -l l ~/Share/day2/coverage/picta/sample1_sample2.unionbedcv
```

### Step 3. Adjust coverage based on the bam ratio in covstats.tab and generate DNAcopy output file
The covstats.tab file also contains a value 'bam_ratio', which is the ratio of the sizes of the two .bam files. We use this to adjust the ratio per window, to account for the fact that one sample might have been sequenced to a higher coverage than the other. This then adjusts the ratio calculated in the previous step by this value. Change the value of 'bam_ratio' in the following command to the value in covstats.tab, and run.  

```
../dif_cover_scripts/from_ratio_per_window_to_prepare_for_DNAcopy_output.sh sample1_sample2.ratio_per_w_CC0_* bam_ratio
```

### Step 4. Create coverage plots in R

This script uses the R package **[DNAcopy[(https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)** to generate output plots of the coverage ratio. DNAcopy uses both the raw coverage ratio, and computes the average coverage ratio of adjacent windows to create the plot. Run this on the server
```
Rscript ../dif_cover_scripts/run_DNAcopy_from_bash.R sample1_sample2.ratio_per_w_CC0_*.log2adj_*
```
Dowload the pdf file to your machine using Filezilla (or equivalent) and inspect it.

### Step 5. Filter only genomic regions with enrichment scores > p.

The script extracts from file *.DNAcopyout fragments with enrichment scores ≥ p and stores them in *.DNAcopyout{p}, (i.e. fragments where read coverage in sample1 is higher than sample2 ), and *.DNAcopyout.down{-p} fragments with enrichment scores ≤-p, (i.e. fragments where coverage in sample2 is higher than sample1 ).

Set p to 2, to filter sites with double coverage in one sample compared to the other.

```
../dif_cover_scripts/from_DNAcopyout_to_p_fragments.sh sample1_sample2.*.DNAcopyout" 2
```

### Step 6. generate histograms for further inspection

```
../dif_cover_scripts/get_DNAcopyout_with_length_of_intervals.sh *.DNAcopyout ref.length.Vk1s_sorted

echo "Generate rough histogram with given precision order"
../dif_cover_scripts/generate_DNAcopyout_len_histogram.sh *.DNAcopyout.len 1

echo "Generate histogram with bins centered at value X reporting scores from [X-0.25 to X+0.25)"
../dif_cover_scripts/generate_DNAcopyout_len_vs_scores_histogram_bin0.5.sh *.DNAcopyout.len
```

### Now run this again for P. reticulata  
### If you finish both species, try adjusting the window size in step 2, and see how this changes the results. Try smaller windows, e.g. 1,000 or larger, e.g. 50,000


   
## 02. SNP based analyses

If you logged out, remember to reactivate the conda environment

```
conda activate /home/ubuntu/miniconda_envs/sexchr
```

On the server in your home set up directories and files

```
cd 
cd day2
mkdir Fst
mkdir GWAS
mkdir SNPden
cp ~/Share/day2/SNPbased/picta_FEMALE.list Fst/
cp ~/Share/day2/SNPbased/picta_MALE.list Fst/
cp ~/Share/day2/SNPbased/picta_sex.list GWAS/
cp ~/Share/day2/SNPbased/picta_FEMALE.list SNPden/
cp ~/Share/day2/SNPbased/picta_MALE.list SNPden/
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
gunzip Poecilia_picta_allchromo_merged.vcf.gz > SNPden/Poecilia_picta_allchromo_merged.vcf
cd SNPden

for SAMPLE in `cat picta_FEMALE.list`;
  do
  ~/Share/day2/SNPbased/bcftools-1.22/bcftools view -a -s ${SAMPLE} -o FEMALE_${SAMPLE}.vcf Poecilia_picta_allchromo_merged.vcf
  bgzip -c FEMALE_${SAMPLE}.vcf > FEMALE_${SAMPLE}.vcf.gz
  ~/Share/day2/SNPbased/bcftools-1.22/bcftools index FEMALE_${SAMPLE}.vcf.gz
  vcftools --gzvcf FEMALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out FEMALE_${SAMPLE}.snpden
done

for SAMPLE in `cat picta_MALE.list`;
  do
  ~/Share/day2/SNPbased/bcftools-1.22/bcftools view -a -s ${SAMPLE} -o MALE_${SAMPLE}.vcf Poecilia_picta_allchromo_merged.vcf
  bgzip -c MALE_${SAMPLE}.vcf > MALE_${SAMPLE}.vcf.gz
   ~/Share/day2/SNPbased/bcftools-1.22/bcftools index MALE_${SAMPLE}.vcf.gz
  vcftools --gzvcf MALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out MALE_${SAMPLE}.snpden
done
```

## Now do the same for Poecilia reticulata
