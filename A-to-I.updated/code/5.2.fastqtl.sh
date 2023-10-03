#!usr/bin/bash

fastQTL --vcf SNP_genoype.vcf.gz \
    --region 1:50000000-60000000 \
    --bed A-to-I.edd.level.qqnorm.gz \
    --cov covarite.txt \
    --permute 1000 10000 \
    --seed 2022 \
    --window 100000 \
    --out eQTL.output.txt
