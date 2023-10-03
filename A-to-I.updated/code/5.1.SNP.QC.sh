#!/usr/bin/sh

## step1

/public/home/zhanjiamin/plink --bfile PPMI.All.SNP --keep sample.ID --make-bed --out PPMI.RNAed.SNP 

## sex check
plink --bfile  PPMI.RNAed.SNP --check-sex 

## impute-sex
plink --bfile  PPMI.RNAed.SNP --impute-sex --make-bed --out  PPMI.RNAed.SNP

## keep Chr1-22
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' PPMI.RNAed.SNP.bim > snp_1_22.txt
plink --bfile PPMI.RNAed.SNP --extract snp_1_22.txt --make-bed --out PPMI.RNAed.SNP_1-22 --noweb
## step2
plink --bfile PPMI.RNAed.SNP_1-22 --mind 0.05 --make-bed --out PPMI.RNAed.SNP_step2 --noweb


## step3
plink --bfile  PPMI.RNAed.SNP_step2 --geno 0.05 --make-bed --out  PPMI.RNAed.SNP_step3 --noweb


## step4
/public/home/zhanjiamin/plink --bfile PPMI.RNAed.SNP_step3 --hwe 1e-6 --make-bed --out PPMI.RNAed.SNP_step4 --noweb


## step5
/public/home/zhanjiamin/plink --bfile PPMI.RNAed.SNP_step4 --maf 0.05 --make-bed --out PPMI.RNAed.SNP.filter --noweb


## PCA
/public/home/zhanjiamin/plink --bfile PPMI.RNAed.SNP.filter --pca 10 --out PPMI.RNAed.SNP.10pca --noweb

## to vcf
/public/home/zhanjiamin/plink --bfile PPMI.RNAed.SNP.filter --recode vcf --out PPMI.RNAed.SNP.filter --noweb
bgzip -f PPMI.RNAed.SNP.filter.vcf
tabix -p vcf  PPMI.RNAed.SNP.filter.vcf.gz
