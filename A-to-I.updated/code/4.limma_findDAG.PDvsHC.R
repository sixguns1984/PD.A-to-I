
library(dplyr)
library(tidyverse) #42; base for data/table manipulation
library(tidyr) #for pivot tables
library(edgeR) #expression normalization
library(limma) #more expression normalization
library(corrplot) #correlation plot matrix


outdir <- "/xx/xx"
dir.create(outdir)

################################################################################## load data

## count input file
Counts_input = " A-to-I editing levels matrix"

## all genes
input_data <- read.delim(Counts_input, row.names=1, check.names = F)

## variance not zero
input_data_ <- input_data[apply(input_data,1,sd)>0.06,] ##0.04
## covariate input file
Cov_input <- "covariate file"

input_cov <- read.delim(Cov_input, check.names = F)
input_cov <- input_cov[,c(1:5)] ## keep PC1-PC5
head(input_cov)
#                  ID CNO DX APPRDX sex EDUCYRS   ENROLLDT    BIRTHDT
# PPMI.3000 PPMI.3000   1         HC      2      1      18  2/01/2011 12/01/1941
#################################################################################################

## match samples in countdata and covdata
samples <- intersect(colnames(input_data), input_cov$ID)
input_data_ <- input_data_[, samples]
input_cov <- input_cov %>% filter(ID %in% samples)

ID_order <- match(colnames(input_data_), input_cov$ID)
input_cov <- input_cov[ID_order,] 
rownames(input_cov) <- input_cov$ID

## check match
identical(colnames(input_data_), rownames(input_cov))


## covariants as factor
input_cov$APPRDX <- as.factor(input_cov$APPRDX)
input_cov$DX <- ifelse(input_cov$APPRDX==1, "PD", "HC")
input_cov$DX <- as.factor(input_cov$DX)
input_cov$gen <- as.factor(input_cov$gen)

###############################################################################################################
sex <- input_cov$gen
RIN <- input_cov$rin
AGE <- input_cov$age
DX <- input_cov$DX

##################################################################
## the library sizes are quite variable between samples, then the voom approach is theoretically more powerful than limma-trend.

design <- model.matrix(~ 0 + sex + AGE + RIN + DX) ## 0 + DX + sex + AGE + RIN
colnames(design) <- gsub("DX", "", colnames(design))  ## gsub("DX", "", colnames(design))
colnames(design)

#########################################################
## Method1 limma-voom
## normalization and filtering in limma
keep <- filterByExpr(input_data_, group = DX,
                     min.count = 0, min.total.count = 0, large.n = 10) ## group = DX; scorstatus
dge <- DGEList(input_data_[keep,])
dge  <- calcNormFactors(dge)

# contr.matrix <- makeContrasts(
#   PD = PD - HC,  ## PD = PD - HC
#   levels = colnames(design))

# contr.matrix
v <- voom(dge, design, plot = F) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)
head(coef(vfit))
# vfit <- contrasts.fit(vfit, contr.matrix)
vfit <- contrasts.fit(vfit, coef=5)

##########################################################
## Method2 don't correct libsize
fit <- lmFit(input_data_, design) 
head(coef(fit))
fit <- contrasts.fit(fit, coef=5)

########################
vfit <- vfit
efit <- eBayes(vfit)
summary(decideTests(efit))
dim(efit)

top.table <- topTable(efit, adjust="BH", sort.by = "P", n = Inf)
head(top.table)
write.csv(top.table, file = paste0(outdir,"A-to-I editing levels different analysis output"), row.names = T, quote = F)

