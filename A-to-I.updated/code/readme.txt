
#When you start using the following code, please make sure read the readme.



############################################################################################################################################
#1. raw fastq file QC and alignment
#1.1 QC 
# Tools: fastp (version 0.20.0)
#Official documentation: https://github.com/OpenGene/fastp
# input file : raw fastq file
# output file : clean fastq
bash 1.1.fastq.QC.sh
#Parameter: 
#raw_data.R1.fastq.gz/raw_data.R2.fastq.gz: raw RNA-seq fastq file
#clean_data.R1.fastq.gz/clean_data.R2.fastq.gz: after quality control RNA-seq fastq file

#1.2 alignment
# Tools: STAR (version 2.7.5a)
#        samtools (version 1.9)
#Official documentation: https://github.com/alexdobin/STAR
#Official documentation: https://github.com/samtools/samtools
# input file : clean fastq file; hg19 genome fastq file
# output file : bam file
bash 1.2.STAR.alignment.sh
#Parameter: --outFileNamePrefix XXX.Aligned.sortedByCoord.out.bam : Output the aligned BAM file.
#           --readFilesIn clean_data.R1.fastq.gz clean_data.R2.fastq.gz : Input RNA-seq fastq files for alignment
#Create an index file for the aligned bam file
# input file : bam file
# output file : index file
bash 1.3.index.sh



############################################################################################################################################
#2.Quantification 
#2.1Detection of A-to-I editing sites levels using RNA-seq data from the PPMI cohort
# Tools: Query_Editing_Level.GRCh37.20161110.pl (version 20161110)
#        perl (version v5.26.2)
# input file : bam file for each sample
# output file : A-to-I editing level files for each sample output
bash 2.1.Quantification.sh
#Parameter: 
#Query_Editing_Level.GRCh37.20161110.pl: Perl code that quantifies the level of A-to-I editing sites.
#rediportal_hg19.txt: rediportal database annotation data (http://srv00.recas.ba.infn.it/atlas/).
#XXX.Aligned.sortedByCoord.out.bam: Input the name of the BAM file.


#2.2 A-to-I editing sites levels QC and normalization
# Tools: python (version v3.7.4)
# input file : raw A-to-I editing sites levels matrix file
# output file : Quality-controlled and normalized A-to-I editing sites levels matrix file 
python3 2.2.QC.normalization.py



############################################################################################################################################
#3.Correlation between A-to-I editing levels and expression of ADARs
# Tools: R (version v4.0.3)
#        prcomp() 
# input file : normalized A-to-I editing sites levels matrix file; normalized A-to-I editing sites levels matrix PC1;ADARs TPM
# output file : Correlation between A-to-I editing levels and expression of ADARs
Rscript 3.Correlation.R



############################################################################################################################################
#4.Differential A-to-I RNA editing and prediction model analysis
# Tools: R (version v4.0.3)
# input file : Quality-controlled A-to-I editing sites levels matrix file; covariate file
# output file : A-to-I editing levels different analysis Results
Rscript 4.limma_findDAG.PDvsHC.R



############################################################################################################################################
#5.Cis-RNA editing quantitative trait loci analysis
#5.1 snp genotype QC
# Tools: vcftools (version v0.1.13)
#       plink (version v1.90b6.20)
#Official documentation: https://vcftools.github.io/
#Official documentation: https://zzz.bwh.harvard.edu/plink/
# input file : raw SNP genotype plink file
# output file : SNP genotype plink file (Quality control)
bash 5.1.SNP.QC.sh

#5.2 eQTL
# Tools: fastQTL (version v2.0)
# input file : normalized A-to-I editing sites levels matrix file; covariate file; Gene expression TPM file; SNP genotype (vcf file)
# output file : A-to-I editing eQTL Results
#Official documentation: https://github.com/francois-a/fastqtl
#Parameter: 
#--vcf SNP Genotype VCF file
#--region set region
#--bed A-To-I Editing level data
bash 5.2.fastqtl.sh



############################################################################################################################################
#6.Mendelian randomization analysis
#6.1.Two-sample Mendelian randomization analysis
# Tools:  R (version v4.0.3)
#        TwoSampleMR (R package; version v0.5.7)
#        ieugwasr (R package; version v0.1.5;API public: http://gwas-api.mrcieu.ac.uk/)
#        plink (version v1.90b6.20)
#Official documentation: https://github.com/MRCIEU/TwoSampleMR
# input file : A-to-I editing eQTL Results; 1KGP EUR SNP genotype(plink file);PD GWAS summary data(access for ieugwasr "ieu-b-7")
# output file : Two-sample Mendelian randomization analysis Results
Rscript 6.1.TMR.R

#6.2.summary data-based mendelian randomization(SMR) analysis
#Official documentation: https://yanglab.westlake.edu.cn/software/smr/#Overview
# Tools: smr (version v1.3.1)
# input file : A-to-I editing eQTL Results;PD GWAS summary data(access for ieugwasr "ieu-b-7");SNP genotype (plink file)
# output file : summary data-based mendelian randomization(SMR) analysis Results
bash 6.2.smr.sh



############################################################################################################################################
#7.Longitudinal A-to-I editing changes associated with the progression of cognitive decline in Parkinson’s disease
# Tools:  R (version v4.0.3)
# input file : normalized A-to-I editing sites levels matrix file (BL,V04,V06,V08); PPMI cohort Longitudinal clinical data (MocA, montreal cognitive assessment);covariate file
# output file : Result of A-to-I editing sites associated with cognitive decline
Rscript 7.longitudinal.R














