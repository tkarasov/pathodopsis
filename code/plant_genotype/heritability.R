library(HLMdiag)
library(lme4)
library(tidyr)
library(gaston)
library(coxme)

varcomp_gemma=function(phenfile,n, Kfile){
  system(paste("gemma -vc 2 -n ",n," -p ", phenfile, " -k ", Kfile, " -o output_var", sep=""))
  con=file("./output/output_var.log.txt")
  line=readLines(con) 
  (line[grep("## pve estimate", line)] %>% strsplit("="))[[1]][2] %>%gsub(" ", "", .)%>% as.numeric()->pve
  (line[grep('## se\\(pve\\)', line)]%>% strsplit("="))[[1]][2]%>%gsub(" ", "", .)%>% as.numeric()->se
  return(c(pve=pve, se=se))
}

# #################################
# Read in metadata Feb. 2020
# #################################
setwd("/ebio/abt6_projects9/pathodopsis_microbiomes/data/")
genot=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/genot.txt")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

#Perform MDS on otu_table
my.responseorig <- data.frame(sqrt(plant_clim$otu_table/1000)) %>%
  dist() %>% cmdscale(eig = F) %>% data.frame()
colnames(my.responseorig) = c("MDS1", "MDS2")
col1 = plant_clim$clim_data$Tour_ID

# Correlation of MDS1 and MDS2 with lat
cor.test(plant_clim$clim_data$Lat, my.responseorig$MDS1) # R = 0.75

cor.test(plant_clim$clim_data$Lat, my.responseorig$MDS2) # R = -0.24

#################################
# Read phenotype information
#################################
my.response1 <- my.responseorig$MDS1 #my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))]
my.response2 <- make.names(as.factor(plant_clim$clim_data$cluster))
my.response3 <- my.responseorig$MDS2
tot = data.frame(sample = plant_clim$clim_data$Sequence_ID, MDS1 = my.response1, MDS2 = my.response3, cluster = my.response2 )
tot_genot = data.frame(sample = t(genot), phen = NA)
tot_genot <- merge(tot_genot, tot, all.x = TRUE )
tot_genot <- tot_genot[,!c("phen")]
phen = tot_genot[,c(3,4,5)]
write.table(phen,"/ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype.txt", 
            row.names = FALSE, col.names = FALSE) 
##rownames(my.plantID) <-my.total.matrix$Plant_ID

# 
# # #################################
# # Read in kinship matrix and perform mixed model regression with kinship matrix when selecting SNPs with 70% genotyping rate
# # #################################
 #x <- read.bed.matrix(
   "/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.LDprune100_10_5.bed" )
x <- read.bed.matrix("/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3.bed")

# # Compute Genetic Relationship Matrix
standardize(x) <- "p"
 K <- GRM(x, autosome.only = TRUE)
#write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.LDprune100_10_5.kinship", row.names = FALSE,col.names = FALSE)
write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink.kinship", row.names = FALSE,col.names = FALSE)
# need an ID variable
dat <- data.frame(MDS1=my.responseorig$MDS1, plant_clim$clim_data, id = as.factor(rownames(my.responseorig)), MDS2 =my.responseorig$MDS2)

#subet myresponses to those in K
keep_rownames <- rownames(my.responseorig)[which(rownames(my.responseorig) %in% rownames(K))]
keep_K <- K[keep_rownames, keep_rownames]
keep_dat <- (dat[which(rownames(my.responseorig) %in% rownames(keep_K)),])
rownames(keep_dat) <- keep_dat$id

gfit1 <- lmekin(scale(MDS1) ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")

gfit2 <- lmekin(MDS2 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")


##############

# # #################################
# # Read in kinship matrix and perform mixed model regression with kinship matrix  when selecting SNPs with 80% genotyping rate
# # #################################
x <- read.bed.matrix("/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.2.bed")

# # Compute Genetic Relationship Matrix
standardize(x) <- "p"
K <- GRM(x, autosome.only = TRUE)
#write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.LDprune100_10_5.kinship", row.names = FALSE,col.names = FALSE)
write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink.kinship", row.names = FALSE,col.names = FALSE)
# need an ID variable
dat <- data.frame(MDS1=my.responseorig$MDS1, plant_clim$clim_data, id = as.factor(rownames(my.responseorig)), MDS2 =my.responseorig$MDS2)

#subet myresponses to those in K
keep_rownames <- rownames(my.responseorig)[which(rownames(my.responseorig) %in% rownames(K))]
keep_K <- K[keep_rownames, keep_rownames]
keep_dat <- (dat[which(rownames(my.responseorig) %in% rownames(keep_K)),])
rownames(keep_dat) <- keep_dat$id

gfit1 <- lmekin(scale(MDS1) ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")

gfit2 <- lmekin(MDS2 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")



# # #################################
# # Read in kinship matrix and perform mixed model regression with kinship matrix  when selecting SNPs with 80% genotyping rate
# # #################################
x <- read.bed.matrix("/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3_mind0.7.bed")

# # Compute Genetic Relationship Matrix
standardize(x) <- "p"
K <- GRM(x, autosome.only = TRUE)
#write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.LDprune100_10_5.kinship", row.names = FALSE,col.names = FALSE)
write.table(K, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink.kinship", row.names = FALSE,col.names = FALSE)
# need an ID variable
dat <- data.frame(MDS1=my.responseorig$MDS1, plant_clim$clim_data, id = as.factor(rownames(my.responseorig)), MDS2 =my.responseorig$MDS2)
rownames(dat) <- dat$Sequence_ID
dat <- dat[which(rownames(dat) %in% colnames(K)),]
#subet myresponses to those in K
keep_rownames <- rownames(dat) #rownames(my.responseorig)[which(rownames(my.responseorig) %in% rownames(K))]
keep_K <- K[keep_rownames, keep_rownames]
keep_dat <- dat #(dat[which(rownames(my.responseorig) %in% rownames(keep_K)),])
rownames(keep_dat) <- keep_dat$id

gfit1 <- lmekin(MDS1 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")

gfit2 <- lmekin(MDS2 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")



















############# Using the bed output by Gautam
x2 <-read.bed.matrix("/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.bed")
#standardize(x2) <- "p"
K2<- GRM(x2, autosome.only = TRUE)


######## Plink-generated matrix by gautam
K2 = read.table("/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.rel")
K2_names = read.table("/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.fam")
rownames(K2) = K2_names[,1]
colnames(K2) = K2_names[,1]
keep_K2 <- as.matrix(K2[keep_rownames, keep_rownames])
gfit2 <- lmekin(MDS1 ~ (1|id), data=keep_dat, varlist=keep_K2, method = "REML")

######## Plink-generated matrix by me
normalize_kinmat <- function(kinmat){
  #normalize kinship so that Kij \in [0,1]
  tmp=kinmat - min(kinmat)
  tmp=tmp/max(tmp)
  tmp[1:454,1:454]
  #fix eigenvalues to positive
  diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
  tmp[1:454,1:454]  
  return(tmp)
}


K2 = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink2.rel")
K2_names = read.table("/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.fam")
rownames(K2) = K2_names[,1]
colnames(K2) = K2_names[,1]
keep_K2 <- as.matrix(K2[keep_rownames, keep_rownames])

#keep_K2 <- normalize_kinmat(keep_K2)
gfit2 <- lmekin(MDS1 ~ (1|id), data=keep_dat, varlist=scale(keep_K2), method = "REML")

# It's saying that it is non-negative definite

######## plink kinship
K3 = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink.kinship")
K3_names = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plink.fam")
rownames(K3) = K3_names[,1]
colnames(K3) = K3_names[,1]
keep_K3 <- as.matrix(K3[keep_rownames, keep_rownames])
gfit3 <- lmekin(MDS1 ~ (1|id), data=keep_dat, varlist=keep_K3, method = "ML")
# Variance = 0.071, Residual Error = 0.0011






# # Fit model only with kinship using REML. The residuals look bad, showing a strong negative correlation between predicted value and residual.
gfit1 <- lmekin(scale(MDS1) ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")

gfit2 <- lmekin(MDS2 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")



gfit1_pred <- keep_K %*% gfit1$coefficients$random$id
resid <- keep_dat$MDS1 - gfit1_pred
plot(gfit1_pred, resid)
fixef(gfit1)
ranef(gfit1) 



# keep_dat <- scale(keep_dat)
# 
# #h2=apply(keep_dat,2, H2nokin, id=keep_K)
# 
# 
# # write a quick function for predicted values from lmekin
# pred_val <- function(keep_K, gfit){
#   keep_K_red <- keep_K[names(gfit$coefficients$random$id), names(gfit$coefficients$random$id)]
#   keep_dat_red <- keep_dat[names(gfit$coefficients$random$id), names(gfit$coefficients$random$id)]
#   return(keep_K_red, keep_dat_red)
# }
# 
# 
# y = keep_dat$MDS1
# 

# 
# data = data.frame(Y = keep_dat$MDS1, id = colnames(keep_K))
# data$id = as.factor(data$id)
# model = lmekin(Y ~ (1|id), data, varlist = keep_K, method = "REML")
# g='id'
# varID=model$vcoef[g]
# varexp=model$vcoef
# varexp=varexp[names(varexp)!=g]
# H2=sum(unlist(varID))/sum(unlist(varID), unlist(varexp),model$sigma^2)
# # why is H2 0.73 for gfit1 and 
# 
# 
# r.squaredLR(gfit1)