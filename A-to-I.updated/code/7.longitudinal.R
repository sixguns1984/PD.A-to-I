#/usr/bin/Rscrip

library(dplyr)
library(reshape2)
library(do)

BL_input <- data.frame(data.table::fread("edd.qqnorm.BL.v3.txt", header=T, check.names=F, sep="\t"))
V04_input <- data.frame(data.table::fread("edd.qqnorm.V04.v3.txt", header=T, check.names=F, sep="\t"))
V06_input <- data.frame(data.table::fread("edd.qqnorm.V06.v3.txt", header=T, check.names=F, sep="\t"))
V08_input <- data.frame(data.table::fread("edd.qqnorm.V08.v3.txt", header=T, check.names=F, sep="\t"))
covariate<-reasd.table("covariate",header=T)

visits_data1 <- merge(BL_input, V04_input,  by="ID")
visits_data2 <- merge(V06_input, V08_input, by="ID")
visits_data <- merge(visits_data1, visits_data2, by="ID")
cov_BL_V02_V04_V06_V08_MoCA<-merge(visits_data,covariate)

cov_BL_V02_V04_V06_V08_MoCA <- cov_BL_V02_V04_V06_V08_MoCA %>% dplyr::filter(V08BLMoCA!="NA")


## MoCA group 
cov_BL_V02_V04_V06_V08_MoCA$status[cov_BL_V02_V04_V06_V08_MoCA$V08BLMoCA< 0] <- "Down"
cov_BL_V02_V04_V06_V08_MoCA$status[cov_BL_V02_V04_V06_V08_MoCA$V08BLMoCA==0] <- "Unchanged"
cov_BL_V02_V04_V06_V08_MoCA$status[cov_BL_V02_V04_V06_V08_MoCA$V08BLMoCA>0] <- "Up"
# table(cov_BL_V02_V04_V06_V08_MoCA$HYstatus)
#table(cov_BL_V02_V04_V06_V08_MoCA$status)
#########################################################################


MoCA_down <- cov_BL_V02_V04_V06_V08_MoCA %>% dplyr::filter(status=="Down")
MoCA_down <- MoCA_down[,c("BL", "V08BLMoCA")]
mean(MoCA_down$BL)
mean(MoCA_down$V08BLMoCA)


MoCA_up <- cov_BL_V02_V04_V06_V08_MoCA %>% dplyr::filter(status=="Up")
MoCA_up <- MoCA_up[,c("BL", "V08BLMoCA")]
mean(MoCA_up$BL)
mean(MoCA_up$V08BLMoCA)

MoCA_unchanged <- cov_BL_V02_V04_V06_V08_MoCA %>% dplyr::filter(status=="Unchanged")
MoCA_unchanged <- MoCA_unchanged[,c("BL", "V08BLMoCA")]
mean(MoCA_unchanged$BL)
mean(MoCA_unchanged$V08BLMoCA)

## Up and unchanged means non progressing patients
non_progressPatient <- cov_BL_V02_V04_V06_V08_MoCA %>% dplyr::filter(status!="Down") ## filter(status=="Unchanged")


#start LMM analysis
for (i in seq(1,22)) {

AI_id <- c()
cor_p_proPD <- c()
cor_rho_proPD <- c()
cor_p_nonproPD <- c()
cor_rho_nonproPD <- c()

# i=22  ##Calculate separately on different chromosomes and merge the results in the end.
PD_edit <- data.table::fread(paste0("visits.qqnorm.data.chr", i, ".txt"), header = T,sep = "\t", stringsAsFactors = F, check.names = F)
PD_edit<-data.frame(PD_edit,check.names=F)
colfile <- read.table("colnames_new.txt", header=T, sep="\t", check.names=F)
PD_edit <- cbind(colfile, PD_edit)
PD_edit <- PD_edit[order(PD_edit$ID1, PD_edit$Visits),]
PD_edit$ID<-Replace(PD_edit$ID,"_BL_",".")
PD_edit$ID<-Replace(PD_edit$ID,"_V04_",".")
PD_edit$ID<-Replace(PD_edit$ID,"_V06_",".")
PD_edit$ID<-Replace(PD_edit$ID,"_V08_",".")
colnames(PD_edit)[1]<-"ID_raw"
colnames(PD_edit)[3]<-"ID1"

Progress_PD_edit <- PD_edit %>% dplyr::filter(ID1 %in% progressPatient$ID)
Progress_PD_edit <- merge(cov_BL_V02_V04_V06_V08_MoCA_change, Progress_PD_edit, by=c("ID1","Visits"))
Progress_PD_edit$Visits <- as.character(Progress_PD_edit$Visits)

non_Progress_PD_edit <- PD_edit %>% dplyr::filter(ID1 %in% non_progressPatient$ID)
non_Progress_PD_edit <- merge(cov_BL_V02_V04_V06_V08_MoCA_change, non_Progress_PD_edit, by=c("ID1","Visits"))
non_Progress_PD_edit$Visits <- as.character(non_Progress_PD_edit$Visits)

Progress_PD_edit$Visits[Progress_PD_edit$Visits=="BL"] <- 0
Progress_PD_edit$Visits[Progress_PD_edit$Visits=="V02"] <- 0.5
Progress_PD_edit$Visits[Progress_PD_edit$Visits=="V04"] <- 1
Progress_PD_edit$Visits[Progress_PD_edit$Visits=="V06"] <- 2
Progress_PD_edit$Visits[Progress_PD_edit$Visits=="V08"] <- 3
Progress_PD_edit$Visits <- as.numeric(Progress_PD_edit$Visits)

non_Progress_PD_edit$Visits[non_Progress_PD_edit$Visits=="BL"] <- 0
non_Progress_PD_edit$Visits[non_Progress_PD_edit$Visits=="V02"] <- 0.5
non_Progress_PD_edit$Visits[non_Progress_PD_edit$Visits=="V04"] <- 1
non_Progress_PD_edit$Visits[non_Progress_PD_edit$Visits=="V06"] <- 2
non_Progress_PD_edit$Visits[non_Progress_PD_edit$Visits=="V08"] <- 3
non_Progress_PD_edit$Visits <- as.numeric(non_Progress_PD_edit$Visits)

AIedit_sites_proPD <- colnames(Progress_PD_edit)[6:ncol(Progress_PD_edit)]
length(AIedit_sites_proPD)
AIedit_sites_nonproPD <- colnames(non_Progress_PD_edit)[6:ncol(non_Progress_PD_edit)]
length(AIedit_sites_nonproPD)

for (j in seq(1, length(AIedit_sites_proPD))){
  AIsite <- AIedit_sites_proPD[j]
  
  AI_id <- c(AI_id, AIsite) ## AI_id
  cor_gene_proPD <- cor.test(~ get(AIsite) + MoCA, method="spearm", exact=FALSE, data = Progress_PD_edit, na.action = na.omit, alternative="t")
  cor_gene_nonproPD <- cor.test(~ get(AIsite) + MoCA, method="spearm", exact=FALSE, data = non_Progress_PD_edit, na.action = na.omit, alternative="t")
  
  cor_p_proPD <- c(cor_p_proPD, cor_gene_proPD$p.value)
  cor_rho_proPD <- c(cor_rho_proPD, cor_gene_proPD$estimate)
  
  cor_p_nonproPD <- c(cor_p_nonproPD, cor_gene_nonproPD$p.value)
  cor_rho_nonproPD <- c(cor_rho_nonproPD, cor_gene_nonproPD$estimate)  
}
result_cor_df <- data.frame("AI_id"=AI_id, "cor_p_proPD"=cor_p_proPD, "cor_rho_proPD"=cor_rho_proPD,
                             "cor_p_nonproPD"=cor_p_nonproPD,
                             "cor_rho_nonproPD"=cor_rho_nonproPD)

cor_p_proPD.adjust <- p.adjust(result_cor_df$cor_p_proPD, method = "BH")
result_cor_df$cor_p_proPD.adjust <- cor_p_proPD.adjust
cor_p_nonproPD.adjust <- p.adjust(result_cor_df$cor_p_nonproPD, method = "BH")
result_cor_df$cor_p_nonproPD.adjust <- cor_p_nonproPD.adjust
write.csv(result_cor_df, file=paste0("result_corMoCA_df_chr", i, ".csv"), quote=F, row.names=F)

}
