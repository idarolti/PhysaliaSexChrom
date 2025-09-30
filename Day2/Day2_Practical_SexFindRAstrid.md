# Day 2 Practical

This practical will cover methods on sex chromosome discovery:

1. Coverage based analysis
2. SNP-based analyses

## 00. Prepare work folder for day 2

```
mkdir day2
cd day2
conda activate sexchr
```

## 01. Coverage based analysis

## 02. SNP based analysis

Set up directories and files

```
mkdir Fst
mkdir GWAS
mkdir SNPden
cp /bin/Ppicta_female.list Fst/
cp /bin/Ppicta_male.list Fst/
cp /bin/Ppicta_sex.list GWAS
cp /bin/Ppicta_female.list SNPden
cp /bin/Ppicta_male.list SNPden
cp Poecilia_picta.g.vcf.gz .
```

### Fst

```
cd Fst
vcftools --gzvcf ../Poecilia_picta.g.vcf.gz --weir-fst-pop Ppicta_female.list --weir-fst-pop Ppicta_male.list --fst-window-size 10000 --out Ppicta.fst.10kb
```

### GWAS

```
cd GWAS
vcftools --gzvcf ../Poecilia_picta.g.vcf.gz --plink --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --out Poecilia_picta.filt
plink --file Poecilia_picta.filt --pheno Ppicta_sex.list --make-bed --out Poecilia_picta.plink --noweb --allow-no-sex
gemma -bfile Poecilia_picta.plink -lm 2 -o Poecilia_picta.gemma
```

### SNP density

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
