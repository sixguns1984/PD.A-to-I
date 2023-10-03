## https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)


## step1
dat <- read.csv("A-to-I editing eQTL Results", header = T)

# remove FDR>=0.05; norsid
dat <- dat %>% filter(FDR<0.05)
dat <- dat %>% filter(grepl("^rs",SNP))
dat$se = sqrt(((dat$beta)^2)/qchisq(dat$pval,1,lower.tail=F)) ##se
exp_dat <- format_aries_mqtl(dat, type = "exposure")

exp_dat_clu <- exp_dat[,c(1,2,12)]
names(exp_dat_clu) <- c("rsid", "pval", "id")

## step2
##SNP clumping
kg_rsid <- read.table("1KGP.v3.EUR.bim", header=F, sep = "\t")

exp_dat_clu  <- exp_dat_clu %>% filter(rsid %in% kg_rsid$V2)

ld_clump_result <- ld_clump(dplyr::tibble(rsid=exp_dat_clu$rsid, pval=exp_dat_clu$pval, id=exp_dat_clu$id), plink_bin = "plink",
                            bfile = "1KGP.v3.EUR")

colnames(ld_clump_result) <- c("SNP", "pval.exposure", "id.exposure")

## step3
exp_dat_ld <- merge(exp_dat, ld_clump_result, by=c("SNP", "pval.exposure", "id.exposure"))

ao <- gwasinfo() #access PD GWAS summary data
head(subset(ao, select=c(trait, id)))
ao[grepl("Parkinson", ao$trait), ]
PD_out_dat <- extract_outcome_data(snps = exp_dat_ld$SNP,outcomes = "ieu-b-7", proxies = TRUE, 
                                   maf_threshold = 0.01,access_token = NULL) #r2>0.8

## step4
ha_dat <- harmonise_data(exposure_dat = exp_dat_ld, outcome_dat = PD_out_dat)
res <- mr(ha_dat, method_list=c("mr_wald_ratio","mr_ivw"))
res_ <- res
res_$p.adjust <- p.adjust(res_$pval, method = "bonferroni")
res_ <- merge(exp_dat_ld[,c(1,2,3,8)], res_, by="id.exposure")
write.csv(res_, file = "TMR_result.full.csv", quote = F, row.names = F)
res_bonf <- res_ %>% filter(p.adjust < 0.05)
write.csv(res_bonf, file = "TMR_result.sig.csv", quote = F, row.names = F)
