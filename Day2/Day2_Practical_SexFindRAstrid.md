# Day 2 Practical Sex chromosome discovery - coverage and SNP based methods

This practical will cover methods on sex chromosome discovery:

1. Coverage based analysis
2. SNP-based analyses

We will base our analysis on the **[SexFindR pipeline](https://sexfindr.readthedocs.io/en/latest/#)** published by **[Grayson et al](https://doi.org/10.1101/2022.02.21.481346)**

## 00. Prepare work folder for day 2

```
mkdir day2
cd day2
conda activate sexchr
```

## 01. Coverage based analysis

## 02. Set up working directories and input files 

Set up directories and files

```
mkdir Fst
mkdir GWAS
mkdir SNPden
cp /home/ubuntu/Share/day2/Ppicta_female.list Fst/
cp /home/ubuntu/Share/day2/Ppicta_male.list Fst/
cp /home/ubuntu/Share/day2/Ppicta_sex.list GWAS/
cp /home/ubuntu/Share/day2/Ppicta_female.list SNPden/
cp /home/ubuntu/Share/day2/Ppicta_male.list SNPden/
cp /home/ubuntu/Share/day2/Poecilia_picta.g.vcf.gz .
```

## 03. Calculate intersex Fst 
We will use **[VCFtools](https://vcftools.github.io)** to calculate Fst (fixation index, a measure of genetic differentiation between populations) between males and female in 10kb windows based on the variant file you learned how to generate yesterday   

```
cd Fst
vcftools --gzvcf ../Poecilia_picta.g.vcf.gz --weir-fst-pop Ppicta_female.list --weir-fst-pop Ppicta_male.list --fst-window-size 10000 --out Ppicta.fst.10kb
cd ..
```

## 04. Run association test for sex with GEMMA   

We will first generate a filtered input file that we will then further format with PLINK and then use as input for **[GEMMA](https://github.com/genetics-statistics/GEMMA)** which runs linear mixed models to test for an association between each SNP and sex.   

```
cd GWAS
vcftools --gzvcf ../Poecilia_picta.g.vcf.gz --plink --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --out Poecilia_picta.filt
plink --file Poecilia_picta.filt --pheno Ppicta_sex.list --make-bed --out Poecilia_picta.plink --noweb --allow-no-sex
gemma -bfile Poecilia_picta.plink -lm 2 -o Poecilia_picta.gemma
cd ..
```

### Calculate male and female SNP density
We will format the SNP variant file with **[bcftools](https://samtools.github.io/bcftools/bcftools.html)** to subset female and male entries into separate files  
We will then calculate SNP density in 10kb windows with VCFtools
```
gunzip Poecilia_picta.g.vcf.gz
mv Poecilia_picta.g.vcf SNPden
cd SNPden

for SAMPLE in `cat Ppicta_female.list`;
  do
  bcftools view -a -s ${SAMPLE} -o FEMALE_${SAMPLE}.vcf Poecilia_picta.g.vcf
  bgzip -c FEMALE_${SAMPLE}.vcf > FEMALE_${SAMPLE}.vcf.gz
  bcftools index FEMALE_${SAMPLE}.vcf.gz
  vcftools --vcf FEMALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out FEMALE_${SAMPLE}.snpden
done

for SAMPLE in `cat Ppicta_male.list`;
  do
  bcftools view -a -s ${SAMPLE} -o MALE_${SAMPLE}.vcf Poecilia_picta.g.vcf
  bgzip -c MALE_${SAMPLE}.vcf > MALE_${SAMPLE}.vcf.gz
  bcftools index MALE_${SAMPLE}.vcf.gz
  vcftools --vcf MALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out MALE_${SAMPLE}.snpden
done
```
