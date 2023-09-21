
#######################################################################################################################
############1.Quantification Detection of A-to-I editing siteslevels using RNA-seq data from the PPMI cohort###########
#######################################################################################################################

#rediportal_hg19.txt: rediportal database (http://srv00.recas.ba.infn.it/atlas/)
perl Query_Editing_Level.GRCh37.20161110.pl rediportal_hg19.txt PPMI_Aligned.sortedByCoord.out.bam outfile

##Query_Editing_Level.GRCh37.20161110.pl

#!/bin/perl
############################################################
#perl script that queries editing level of known sites in a BAM file

use warnings;
use strict;
require "parse_pileup_query.pl"; #NEED PARSE PILEUP LIBRARY

if (@ARGV != 3) {
	die "need to provide 3 input:Edit Site list, INDEXED BAM alignment file and output file name\n";
}
my ($inputfile, $bamfile, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2]);

#GLOBAL VARIABLES - PLEASE MODIFY THESE

my $minbasequal = 20; # MINIMUM BASE QUALITY SCORE
my $minmapqual = 1; # MINIMUM READ MAPPING QUALITY SCORE
my $sampath = "samtools"; #PATH TO THE SAMTOOLS EXECUTABLE
my $genomepath = "GRCh37.primary_assembly.genome.fa"; #PATH TO REFERENCE GENOME

my $offset = 33; #BASE QUALITY SCORE OFFSET - 33 FOR SANGER SCALE, 64 FOR ILLUMINA SCALE

##END GLOBAL VARIABLES

my $bedtemp = join '', $outputfile, '.bed';
system("awk \'\$1\!\=\"chromosome\"\{print \$1\"\t\"\$2-1\"\t\"\$2\}\' $inputfile \> $bedtemp");
my $piletemp = join '', $outputfile, '.pileup';
system("$sampath mpileup -A -B -d 1000000 -q $minmapqual -Q $minbasequal -f $genomepath -l $bedtemp $bamfile \> $piletemp");

my %sitehash;
open (my $PILEUP, "<", $piletemp);
while(<$PILEUP>) {
	chomp;
	my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
	my $location = join '_', $chr, $position;
	my ($refnuccount, $acount, $tcount, $ccount, $gcount) = &parse_pileup($_, $minbasequal, $offset);# parse each line of pileup
	my $counts = join ',', $refnuccount, $ccount, $gcount;
	$sitehash{$location} = $counts;
}
system("rm $bedtemp");
system("rm $piletemp");

open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);
print $OUTPUT "#chrom\tposition\tgene\tstrand\tannot1\tannot2\tcoverage\teditedreads\teditlevel\n";

while (<$INPUT>) { #READ IN LIST OF KNOWN EDITED SITES AND QUERY EDITING STATUS
	chomp;
	my @fields = split;
	next if ($fields[0] eq 'chromosome');
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $location = join '_', $chr, $position;
	my ($gene, $strand, $annot1, $annot2) = ($fields[2], $fields[3],$fields[4], $fields[5]);

	if ($sitehash{$location}) { #PRINT OUT RESULT
		my ($refcount, $ccount, $gcount) = split(/\,/,$sitehash{$location});
		my ($newcov, $newmismatch) = (0,0);
		if ($strand eq '+') {
			$newmismatch = $gcount;
		} else {
			$newmismatch = $ccount;
		}
		$newcov = $refcount + $newmismatch;
		if ($newcov) {		
			my $varfreq = 0;
			$varfreq = sprintf("%.3f", $newmismatch/$newcov);
			print $OUTPUT "$fields[0]\t$fields[1]\t$newcov\t$newmismatch\t$varfreq\n";
		} else {
			print $OUTPUT "$fields[0]\t$fields[1]\t0\t0\tN/A\n";
		}
	} else {
		print $OUTPUT "$fields[0]\t$fields[1]\t0\t0\tN/A\n";
	}
}
close $INPUT;	
close $OUTPUT;


###################################################################################
###########2.Correlation between A-to-I editing levels and expression of ADARs#####
###################################################################################

#/public/apps/bioinformatics/anaconda3/bin/Rscrip
#AtoI:  A-to-I editing levels ()
PCA <- prcomp(AtoI, center=T, scale=T)
var2 <- PC5$x[,1:2] %>% as.data.frame()
#ADARs gene
cor.test(log2(ADARs_PCA$ENSG00000160710.11+1), ADARs_PCA$PC1)
cor.test(log2(ADARs_PCA$ENSG00000197381.11+1), ADARs_PCA$PC1)
cor.test(log2(ADARs_PCA$ENSG00000185736.11+1), ADARs_PCA$PC1)


###################################################################################
################3.Differential A-to-I RNA editing and prediction model analysis####
###################################################################################

#/public/apps/bioinformatics/anaconda3/bin/Rscript
#prediction model
#data: A-to-I editing levels
fit <- lm(data[,11] ~ data$age + data$gen + data$rin) #+ covariate$batch
summary(fit)
data$residuals <- fit$residuals
roc(data$DX, data$AtoIsites, ci=T, auc=T, plot=T, print.auc=T, direction=">")


####################################################################################
#####################4.The quality control steps on genotype data###################
####################################################################################

plink --bfile data --hwe/--keep/--gnom/--maf --make-bed --out
#--hwe/--keep/--gnom/--maf: QC step


####################################################################################
#####################5.Cis-RNA editing quantitative trait loci analysis#############
####################################################################################
#/bin/bash
#fastQTL (Official documentation: https://github.com/francois-a/fastqtl)
#--vcf SNP Genotype VCF file; --region region;--bed A-To-I Editing level data
fastQTL --vcf PPMI.RNAed.SNP.filter.vcf.gz \
    --region 1:50000000-60000000 \
    --bed PPMI.edd.qqnorm.gz \
    --cov PPMI.cov.txt \
    --permute 1000 10000 \
    --seed 2022 \
    --window 100000 \
    --out permutations.1.6.txt



######################################################################################
#####################6.Mendelian randomization analysis###############################
######################################################################################
#/public/apps/bioinformatics/anaconda3/bin/Rscrip
#Two sample MR
exp_dat <- format_aries_mqtl(dat, type = "exposure")
kg_rsid <- read.table("/public1/data/zhanjiamin/PPMI_Methylation/PPMI_120_Methylation/1kg.v3/EUR.bim", header=F, sep = "\t")
exp_dat_clu <- exp_dat[,c(1,2,12)]
names(exp_dat_clu) <- c("rsid", "pval", "id")
exp_dat_clu  <- exp_dat_clu %>% filter(rsid %in% kg_rsid$V2)
ld_clump_result <- ld_clump(dplyr::tibble(rsid=exp_dat_clu$rsid, pval=exp_dat_clu$pval, id=exp_dat_clu$id), plink_bin = "/public/home/zhanjiamin/plink",
                            bfile = "/public1/data/zhanjiamin/PPMI_Methylation/PPMI_120_Methylation/1kg.v3/EUR")
exp_dat_ld <- merge(exp_dat, ld_clump_result, by=c("SNP", "pval.exposure", "id.exposure"))
#read GWAS summary data
ao <- gwasinfo()
head(subset(ao, select=c(trait, id)))
ao[grepl("Parkinson", ao$trait), ]
PD_out_dat <- extract_outcome_data(snps = exp_dat_ld$SNP,outcomes = "ieu-b-7", proxies = TRUE, 
                                   maf_threshold = 0.01,access_token = NULL) #r2>0.8
ha_dat <- harmonise_data(exposure_dat = exp_dat_ld, outcome_dat = PD_out_dat)
res <- mr(ha_dat, method_list=c("mr_wald_ratio","mr_ivw"))
res$p.adjust <- p.adjust(res$pval, method = "bonferroni")

#SMR
smr-1.3.1 --eqtl-summary PPMI.FastedQTL.chr1-22.clinical-cov.nominalP.txt --fastqtl-nominal-format --make-besd --out PPMI.FastedQTL.chr1-22.clinical-cov.nominalP 

$smr-1.3.1 --bfile PPMI.RNAed.SNP.filter \
--gwas-summary $nallsEtAl2019.rsid.ma \
--beqtl-summary $PPMI.FastedQTL.chr1-22.clinical-cov.nominalP \
--maf 0.05 \
--out PPMI.FastedQTLsmr \
--thread-num 10 \
--heidi-mtd 1 \


#######################################################################################################################
###########7.Analysis for longitudinal changes of A-to-I editing and cognitive decline in patients with PD#############
#######################################################################################################################

#/public/apps/bioinformatics/anaconda3/bin/Rscrip
site<-read.table("PPMI.HC.PD.edMat.20cov.400samps.NA.0.4.Sites.txt",header=T)
phenodata <- read.table("PPMI.edQTL.longvisitsID.LMM.cov.txt", header=T, check.name=F, sep="\t")
eddata <- data.table::fread("PPMI.HC.PD.edMat.20cov.400samps.BL_V04_V06_V08.new.txt",header=T)
eddata<-data.frame(eddata,check.names=F)
eddata$ID<-gsub(":","_",eddata$ID)
row.names(eddata)<-eddata$ID
eddata<-eddata[,-1]
eddata<-t(eddata)
eddata<-data.frame(eddata,check.names=F)
eddata$ID<-row.names(eddata)
lmdata <- merge(phenodata, eddata, by="ID")
dim(lmdata)

covdata <- lmdata[1:11]
edexp <- lmdata[12:ncol(lmdata)]
edexp <- scale(edexp)
data2 <- data.frame(covdata,edexp)

varname <- colnames(data2)[12:ncol(data2)]

result1<-c()
result2<-c()
for (m in varname)
	{
		print(paste0("start site:",which(varname==m),", all site: 243094"))
    lmdata<-data2[,c(m,'age','YearsafM0','Sex','sample','Durbase','RIN', "Status")]
		lmdata$Durbase <- as.numeric(lmdata$Durbase)
		fitlme<-lmer(lmdata[,1] ~ 1 + Status * YearsafM0 + age + RIN * YearsafM0 + Sex + Durbase + (1 | sample), na.action = na.omit, data = lmdata) 
		res <- summary(fitlme)
		pvalue<-res$coefficients['StatusPD:YearsafM0','Pr(>|t|)']
		beta<-res$coefficients['StatusPD:YearsafM0','Estimate']
		agid<-m
		result1<-data.frame(agid, beta, pvalue)
    result2<-rbind(result1,result2)
	}
write.csv(result2,"PPMI.eddata.LMM.result.all.csv", quote=F, row.names=F)



