#/use/bin/Rscrip

#PC1 and ADARs expression
AtoI<-data.table::fread("A-to-I editing levels matrix",header=T)
ADARs_exp <- read.table("ADARs TPM matrix", header = T, sep = "\t")
AtoI<-data.frame(AtoI)
PCA <- prcomp(AtoI, center=T, scale=T)
PCA1_2 <- PCA$x[,1:2] %>% as.data.frame()
PCA1_2$ID <- row.names(PCA1_2)
ADARs_PCA <- merge(ADARs_exp, PCA1_2, by="ID")
cor.test(log2(ADARs_PCA$ENSG00000160710.11+1), ADARs_PCA$PC1)
cor.test(log2(ADARs_PCA$ENSG00000197381.11+1), ADARs_PCA$PC1)
cor.test(log2(ADARs_PCA$ENSG00000185736.11+1), ADARs_PCA$PC1)
