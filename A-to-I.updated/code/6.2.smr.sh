#SMR
#!/usr/bin/sh
#Official documentation: https://yanglab.westlake.edu.cn/software/smr/#Overview


## fastqtl format: gene, snp, distance, pvalue, slope
# ENSG00000237438.1 indel:4D_22_16518157 -999303 0.542909 -0.0510761
awk -F "," '{print$1"\t"$6"\t"$7"\t"$8"\t"$9}' eQTL-result > eQTL-result.nominalP.txt


smr-1.3.1 --eqtl-summary eQTL-result.nominalP.txt --fastqtl-nominal-format --make-besd --out eQTL-result.nominalP 
## Next, it is necessary to organize the esi file and epi file based on the site annotation as shown in the following example:
## esi format: chr, rsid, 0, position, A1, A2, maf_A1
# 1    rs1001  0   744055  A   G   0.23
# 1    rs1002  0   765522  C   G   0.06
# 1    rs1003  0   995669  T   C   0.11

## epi format: chr, probe, 0, position, geneSymbol, strand
# 1    probe1001   0   924243  Gene01  +
# 1    probe1002   0   939564  Gene02  -
# 1    probe1003   0   1130681 Gene03  -

#gwas summary data: PD.gwas.ma
#example
SNP	A1	A2	freq	b	se	p	N
rs11020170	T	C	0.9931	0.1575	0.1783	0.3771	12517
rs116406626	A	G	0.9336	0.0605	0.0456	0.1846	468692
...
...

################################################################
smr-1.3.1 --bfile SNP_genotype \
--gwas-summary PD.gwas.ma \
--beqtl-summary eQTL-result.nominalP \
--maf 0.05 \
--out smr_result \
--thread-num 10 \
--heidi-mtd 1 \
###################################################################
