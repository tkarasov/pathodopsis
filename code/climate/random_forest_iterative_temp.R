# This function does feature selection using the caret and mbench packages then does random forest modeling
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
library(caret)
library(mlbench)
library(tidyr)
library(Hmisc)
library(dplyr)
#library(randomForest)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(phyloseq)
library(gplots)
library(ggthemes)
library(MLeval)
library(plotROC)
library(GROAN)
library(GENESIS)
library(coxme)
library(gaston)
library(MuMIn)
library(pROC)
source("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/climate/generic_random_forest.R")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/color_pal.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/theme_pubs.rds")


# #################################
# Read in metadata Feb. 2020
# #################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

#Perform MDS on otu_table
my.responseorig <- data.frame(sqrt(plant_clim$otu_table/1000)) %>% 
  dist() %>% cmdscale(eig = F) %>% data.frame()
colnames(my.responseorig) = c("MDS1", "MDS2")
col1 = plant_clim$clim_data$Tour_ID

# Correlation of MDS1 and MDS2 with lat
cor.test(plant_clim$clim_datal$Lat, my.responseorig$MDS1) # R = 0.75

cor.test(plant_clim$clim_datal$Lat, my.responseorig$MDS1) # R = -0.24

#################################
# Set response variable as MDS1 or as cluster
#################################
my.response1 <- my.responseorig$MDS1 #my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))] 
my.response2 <- make.names(as.factor(plant_clim$clim_data$cluster))
my.response3 <- my.responseorig$MDS2
##rownames(my.plantID) <-my.total.matrix$Plant_ID

#################################
# Random forest on cluster
#################################
my.total.matrix.num <- generate_my_total_matrix(my.response2, plant_clim)
preprocess1 <- preprocess_data(my.total.matrix.num)
x_full <- preprocess1[1][[1]]
x_train <- preprocess1[2][[1]]
x_test <- preprocess1[3][[1]]
train_ind <- preprocess1[4][[1]]

# Subset training and test response variable
y = as.factor(my.total.matrix.num$response)
y_train <- y[train_ind]
y_test_cluster <- y[-train_ind]

#Select features from full data
subsets <- c(1:20,25, 30, 33)
rfProfile_class = my_feat_selec(x_full, y, subsets)
my.predictors = predictors(rfProfile_class)
x_train = x_train[,my.predictors]
x_test = x_test[,my.predictors]

# Random Forest model (caret) on training data
rf_train.output_cluster <-my_random_forest(my.predictors, x=x_train, y=y_train, classification = TRUE)
rf_test.output_cluster<-predict(rf_train.output_cluster, newdata = x_test, type = "prob")

importance <- varImp(rf_train.output_cluster, scale = TRUE) #varImp(rf_train.output_cluster$finalModel, scale=TRUE)
import_class <- plot(importance)
import_1 <- importance$importance %>% filter(Overall > 1) %>% arrange(desc(Overall))
import_2 <- data.frame(import_1[1:15,])
rownames(import_2) <- rownames(import_1)[1:15]
import_2$Variable <- rownames(import_2)
colnames(import_2)[1] <- "Importance"

# plot relative importance
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/feature_importance.pdf",
    useDingbats = FALSE, fonts = "ArialMT", width = 3.5, height = 2)
importance <- ggplot(import_2, aes(x = reorder(Variable, Importance), 
                        y = Importance)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance"
    #title = "Feature Importance: <Model>"
  )
#dev.off()
# 82% prediction accuracy. PDSI and vapor pressure, site type were the best predictors, mtry = 31

PDSI_1 <- ggplot(data = my.total.matrix.num, aes(x = response, y = PDSI)) +
  #geom_line() +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, aes(color = response)) + 
  scale_colour_manual(name = "Cluster", values = c('#d95f02', '#1b9e77')) +
  geom_jitter(data = my.total.matrix.num, aes(x = response, y = PDSI, color = response), 
              cex = 1, alpha = 0.5, width = 0.1, height = 0.1, show.legend = FALSE) +
  theme_pubs +
  xlab("") + 
  ylab("PDSI")

plot_grid(importance, PDSI_1)

dev.off()
#
glm(plant_clim$clim_data$cluster ~ plant_clim$clim_data$PDSI + plant_clim$clim_data$Lat, family = "binomial")

model1 <- lm(my.response1 ~ plant_clim$clim_data$PDSI + plant_clim$clim_data$Lat)
model2 <- lm(my.response1 ~ plant_clim$clim_data$Lat)

#################################
# Correlation between variables
#################################
nums <- unlist(lapply(my.total.matrix.num, is.numeric))  
correlationMatrix <- cor(my.total.matrix.num[,nums], use = "pairwise.complete.obs")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/env_heatmap.pdf",
    useDingbats = FALSE, fonts = "ArialMT")

heatmap.2(correlationMatrix, scale = "none", density.info="none", 
          trace="none", dendrogram = "none" )#, col = Colors)

dev.off()

# summarize the correlation matrix
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)

#################################
# Random forest on MDS1
#################################
my.total.matrix.num <- generate_my_total_matrix(response = my.response1, plant_clim)
preprocess1 <- preprocess_data(my.total.matrix.num)
x_full <- preprocess1[1][[1]]
x_train <- preprocess1[2][[1]]
x_test <- preprocess1[3][[1]]
train_ind <- preprocess1[4][[1]]

# Subset training and test response variable
y = my.total.matrix.num$response
y_train <- y[train_ind]
y_test <- y[-train_ind]

#Select features from full data
subsets <- c(1:20,25, 30, 33)
rfProfile_reg = my_feat_selec(x_full, y, subsets)
my.predictors = predictors(rfProfile_reg)
x_train = x_train[,my.predictors]
x_test = x_test[,my.predictors]

# Random Forest model (caret) on training data
rf_train.output <-my_random_forest(my.predictors, x=x_train, y=y_train, classification = FALSE)
rf_test.output<-predict(rf_train.output, newdata = x_test)
importance <- varImp(rf_train.output$finalModel, scale=TRUE)

import_class <- plot(importance)
import_1 <- importance$Overall %>% filter(Overall > .2) %>% arrange(desc(Overall))
import_2 <- data.frame(import_1[1:15,])
rownames(import_2) <- rownames(import_1)[1:15]
import_2$Variable <- rownames(import_2)
colnames(import_2)[1] <- "Importance"

# plot relative importance
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/feature_importance.pdf",
    useDingbats = FALSE, fonts = "ArialMT", width = 3.5, height = 3.5)
importance <- ggplot(import_2, aes(x = reorder(Variable, Importance), 
                                   y = Importance)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance"
    #title = "Feature Importance: <Model>"
  )


yup <- plant_clim$clim_data
ggplot(aes(x = my.response1, y=my.response3, col = plant_clim$clim_data$cluster))

#Rsquared = 0.58 with mtry = 49 and 45 predictor variables

# #################################
# Read in kinship matrix and perform mixed model regression with kinship matrix
# #################################
x <- read.bed.matrix(
    "/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.LDprune100_10_5.bed" )

# Compute Genetic Relationship Matrix
K <- GRM(x, autosome.only = FALSE)
dim(K)


# load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.vcf.kin0')
# file.king <- "/ebio/abt6_projects9/pathodopsis_microbiomes/data/poolsGVCF.filtered_snps_final.PASS.bi.vcf.kin0"
# KINGmat <- kingToMatrix(file.king, estimator="Kinship")


# need an ID variable
dat <- data.frame(MDS1=my.responseorig$MDS1, plant_clim$clim_data, id = rownames(my.responseorig))

#subet myresponses to those in K
keep_rownames <- rownames(my.responseorig)[which(rownames(my.responseorig) %in% rownames(K))]
keep_K <- K[keep_rownames, keep_rownames]
keep_dat <- dat[which(rownames(my.responseorig) %in% rownames(keep_K)),]
rownames(keep_dat) <- keep_dat$id
# keep_dat <- scale(keep_dat)

# write a quick function for predicted values from lmekin
pred_val <- function(keep_K, gfit){
  keep_K_red <- keep_K[names(gfit$coefficients$random$id), names(gfit$coefficients$random$id)]
  keep_dat_red <- keep_dat[names(gfit$coefficients$random$id), names(gfit$coefficients$random$id)]
  return(keep_K_red, keep_dat_red)
}

#only kinship
gfit1 <- lmekin(MDS1 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")
gfit1_pred <- keep_K %*% gfit1$coefficients$random$id
resid <- keep_dat$MDS1 - gfit1_pred
plot(gfit1_pred, resid)
fixef(gfit1)
ranef(gfit1) 
r.squaredLR(gfit1)

#kinship and drought
gfit2 <- lmekin(MDS1 ~ (1|id) + scale(PDSI), data=keep_dat, method = "REML")
red <- pred_val(keep_K, gfit2)
gfit2_pred <- keep_K[names(gfit2$coefficients$random$id), names(gfit2$coefficients$random$id)] %*% gfit2$coefficients$random$id + scale(keep_dat[names(gfit2$coefficients$random$id),]$PDSI) %*% gfit2$coefficients$fixed[2]
resid <- keep_dat[names(gfit2$coefficients$random$id),]$MDS1 - gfit2_pred
plot(gfit2_pred, resid, pch = 20)
r.squaredLR(gfit2, null.RE = TRUE)

#kinship and Lat
gfit3 <- lmekin(MDS1 ~  (1|id) + scale(Lat), data=keep_dat, varlist=keep_K, method = "REML")
r.squaredLR(gfit3, null.RE = TRUE)

#kinship and Lat and Drought
gfit4 <- lmekin(MDS1 ~  (1|id) + scale(Lat) + scale(PDSI), data=keep_dat, varlist=keep_K, method = "REML")
r.squaredLR(gfit4, null.RE = TRUE)

#################################
# Plot residuals
#################################


# 
# ###
# # openmx
# Sys.setenv(OMP_NUM_THREADS=32)
# library(OpenMx)
# 
# #greml
# library(qgg)
# fitG <- greml(MDS1 ~ (1|id), data=keep_dat, varlist=keep_K, method = "REML")
# 
# #################################
# ASV x ASV association with drought
#################################
asv_relation <- function(asv){
  #keep_otu <- plant_clim$otu_table[keep_dat$id,]
  rownames(keep_otu) <- keep_dat$id
  #keep_otu <- keep_otu[keep_dat[which(keep_dat$Lat<45 & keep_dat$Lat>30),]$id,]
  #keep_otu <- keep_otu[keep_dat$id,]
  model1 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat$Lat)
  model2 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat$PDSI)
  model3 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat$PDSI + keep_dat$Lat)
  pval_lat <- summary(model1)$coefficients[2,4]
  print("poo")
  pval_pds <- summary(model2)$coefficients[2,4]
  pval_lat_tog <- summary(model3)$coefficients[3,4]
  pval_pds_tog <- summary(model3)$coefficients[2,4]
  return(c(p.adjust(pval_lat, "BH"), p.adjust(pval_pds, "BH"), p.adjust(pval_lat_tog, "BH"), p.adjust(pval_pds_tog, "BH")))
}
N = length(colnames(plant_clim$otu_table))
otu_sig <- data.frame(pval_lat = numeric(N), 
    pval_pds = numeric(N),
    pds_tog = numeric(N),
    lat_tog = numeric(N),
    lat_corr = numeric(N))
rownames(otu_sig) = colnames(plant_clim$otu_table)



keep_dat <- dat[which(rownames(my.responseorig) %in% rownames(keep_K)),]
rownames(keep_dat) <- keep_dat$id
keep_otu <- plant_clim$otu_table[keep_dat$id,]

for(asv in colnames(plant_clim$otu_table)){
  print(asv)
  pvals <- asv_relation(asv)
  otu_sig[asv,c(1:4)] = pvals
  otu_sig[asv, 5] = cor.test(as.numeric(keep_otu[,asv]),  keep_dat[rownames(keep_otu),]$Lat)$estimate
}

#what percentage are correlated with latitude? 
lat_only <- length(which(otu_sig$pval_pds<0.01)) / dim(otu_sig)[1] #33%
pdsI_lat_correct <- length(which(otu_sig$pds_tog<0.01)) / dim(otu_sig)[1] #10%


#plot correlation of ASVs with latitude or PDSI

keep_otu_lat_pos <- keep_otu[,which(otu_sig$pval_lat<0.01 & otu_sig$lat_corr>0)]/1000
keep_otu_lat_neg <- keep_otu[,which(otu_sig$pval_lat<0.01 & otu_sig$lat_corr<0)]/1000
neg <- colnames(keep_otu[,which(otu_sig$pval_lat<0.01 & otu_sig$lat_corr<0)])
pos <- colnames(keep_otu[,which(otu_sig$pval_lat<0.01 & otu_sig$lat_corr>0)])
all_tax <- data.frame(plant_clim$tax_table)$Family
pos_tax <- subset_taxa(plant_clim, pos)
#data.frame(plant_clim$tax_table)$Family[pos,]

#################################
# ROC curves for cluster assignment
#################################
all_ROC = cbind(y_test_cluster, rf_test.output_cluster)
y_test = addLevel(all_ROC$y_test, newlevel="X4")
#g1_3 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X1"))
#g3_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X3"))
g1_2 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X1"))
g2_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X2"))
#g2_3 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X2")) 
#g3_2 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X3")) 

full_comp <- ggroc(list("Class1vs2" = g1_2, "Class2vs1" = g2_1, 
                        "Class2vs3" = g2_3, "Class3vs2" = g3_2,
                        "Class1vs3" = g1_3, "Class3vs1" = g3_1), legacy.axes=TRUE) +
  scale_color_viridis_d() +
  geom_abline(intercept = 0, slope = 1,
              color = "darkgrey", linetype = "dashed") +
  theme_pubs +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 20)))

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/ASV_class_ROC.pdf", 
  useDingbats = FALSE, fonts = "ArialMT")
full_comp
dev.off()

auc(g1_3)
#0.9343
auc(g3_1)
#0.9594
auc(g1_2)
#0.8269
auc(g2_1)
# 0.5898
auc(g3_2)
# 0.7385
auc(g2_3)
#0.6639

#################################
# Logistic regression for cluster assignment
#################################
# logistic regression of 1 vs 3
tog <- data.frame(my.total.matrix)
test_tog <- tog[-train_ind,] %>% filter(cluster%in%c("1", "2"))
train_tog <- tog[train_ind,] %>% filter(cluster%in%c("1", "2"))
tog <- tog %>% filter(cluster%in%c("1", "3"))
tog$recode = c(0)
tog[tog$cluster=="1",]$recode = 0
tog[tog$cluster!="3",]$recode = 1
mylogit.PDSI_vpd <- glm(recode ~ PDSI + vpd + srad + Elevation + tmax + def +  ws + Site_type + pet, data = tog, family = "binomial")
mylogit.latitude <- glm(recode ~ Lat, data = tog, family = "binomial")
newdat <- data.frame(PDSI=seq(min(tog$PDSI, na.rm = TRUE), max(tog$PDSI, na.rm = TRUE)),
                              vpd = seq(min(tog$vpd, na.rm = TRUE), max(tog$vpd, na.rm= TRUE),
                                        len=100))
newdat$vs = predict(mylogit.PDSI_vpd , newdata=newdat, type="response")

AIC(mylogit.PDSI_vpd, mylogit.latitude)

logit_test <- glm(data = train_tog,
                  recode ~ PDSI + vpd + srad ,
                    family = "binomial")

logit_Lat_test <- glm(data = train_tog, 
                  recode ~ Lat,
                  family = "binomial")


test_tog$pred <- predict(logit_test, newdata=test_tog, type = "response")
test_tog$pred_lat <- predict(logit_Lat_test, newdata=test_tog, type = "response")
logit_roc <- roc(test_tog, response = cluster, predictor=pred)
lat_roc <- roc(test_tog, response = cluster, predictor=pred_lat)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/logit_RF_ROC_comparison.pdf")
ggroc(list("Random Forest" = g1_3, "Full logit" = logit_roc, "logit Latitude Only" = lat_roc), legacy.axes = TRUE) +
  scale_color_viridis_d() +
  geom_abline(intercept = 0, slope = 1,
              color = "darkgrey", linetype = "dashed") +
  theme_pubs +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 20))) +
  theme(legend.key=element_blank())

dev.off()



# #################################
# # Step 5: Random forest across OTUs
# ################################
# 
# 
# 
# 
# 
# #################################
# # calculate correlation matrix. Only possible with numeric predictors
# #################################
# 
# nums <- unlist(lapply(my.total.matrix.num, is.numeric))  
# correlationMatrix <- cor(my.total.matrix.num[,nums], use = "pairwise.complete.obs")
# 
# pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/env_heatmap.pdf",
#     useDingbats = FALSE, fonts = "ArialMT")
# 
# heatmap.2(correlationMatrix, scale = "none", density.info="none", 
#           trace="none", dendrogram = "none", col = Colors)
# 
# dev.off()
# # summarize the correlation matrix
# # find attributes that are highly corrected (ideally >0.75)
# highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)
# 
# 
# 
# 
# 
# 
# 
# 
# # combine all data sets
# 
# col1 = x$Tour_ID#as.factor(sample_data(only_ath)$Clim)
# col2 = fin_predictors$PDSI
# 
# p1 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color = col1), cex = 3) +
#   scale_color_viridis_d() +
#   theme_bw()
# 
# p2 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color = col2), cex = 3) +
#   scale_colour_gradient2() +
#   theme_bw()
# 
# 
# plot_grid(import_MDS1, p2)

