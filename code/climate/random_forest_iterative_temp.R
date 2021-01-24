# This function does feature selection using the caret and mbench packages then does random forest modeling
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
library(caret)
library(mlbench)
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
source("/ebio/abt6_projects9/pathodopsis_microbiomes/generic_random_forest.R")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/color_pal.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/theme_pubs.rds")


# #################################
# Read in metadata Feb. 2020
# #################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")


#Perform MDS on otu_table
my.responseorig <- data.frame(sqrt(plant_clim$otu_table)) %>% 
  dist() %>% cmdscale(eig = F) %>% data.frame()
colnames(my.responseorig) = c("MDS1", "MDS2")
col1 = plant_clim$clim_data$Tour_ID

#################################
# Set response variable as MDS1 or as cluster
#################################
my.response1 <- my.responseorig$MDS1 #my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))] 
my.response2 <- make.names(as.factor(plant_clim$clim_data$cluster))
my.response3 <- make.names(as.factor(plant_clim$clim_data$MDS2))
##rownames(my.plantID) <-my.total.matrix$Plant_ID

#################################
# Random forest on cluster
#################################
my.total.matrix.num <- generate_my_total_matrix(response = my.response1)
preprocess1 <- preprocess_data(my.total.matrix.num)
x_full <- preprocess1[1]
x_train <- preprocess1[2]
x_test <- preprocess1[3]
train_ind <- preprocess1[4]

# Subset training and test response variable
y = make.names(as.factor(my.total.matrix.num$cluster))
y_train <- y[train_ind]
y_test <- y[-train_ind]

#Select features from full data
subsets <- c(1:20,25, 30, 33)
rfProfile = my_feat_selec(x_full, y, subsets)
my.predictors = predictors(rfProfile)
x_train = x_train[,my.predictors]
x_test = x_test[,my.predictors]

# Random Forest model (caret) on training data
rf_train.output <-my_random_forest(my.predictors, x=x_train, y=y_train)
rf_test.output<-predict(rf_train.output, newdata = x_test, type = "prob")

importance <- varImp(rf_train.output$finalModel, scale=TRUE)
import_class <- plot(importance)

#################################
# Random forest on MDS1
#################################
my.total.matrix.num <- generate_my_total_matrix(response = my.response2)
preprocess1 <- preprocess_data(my.total.matrix.num)
x_full <- preprocess1[1]
x_train <- preprocess1[2]
x_test <- preprocess1[3]
train_ind <- preprocess[4]

# Subset training and test response variable
y = make.names(as.factor(my.total.matrix.num$cluster))
y_train <- y[train_ind]
y_test <- y[-train_ind]

#Select features from full data
subsets <- c(1:20,25, 30, 33)
rfProfile = my_feat_selec(x_full, y, subsets)
my.predictors = predictors(rfProfile)
x_train = x_train[,my.predictors]
x_test = x_test[,my.predictors]

# Random Forest model (caret) on training data
rf_train.output <-my_random_forest(my.predictors, x=x_train, y=y_train)
rf_test.output<-predict(rf_train.output, newdata = x_test, type = "prob")
importance <- varImp(rf_train.output$finalModel, scale=TRUE)
import_class <- plot(importance)

#################################
# Random forest on MDS2
#################################




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
dat <- data.frame(MDS1=my.responseorig$MDS1, my.total.matrix.num, id = rownames(my.responseorig))

#subet myresponses to those in K
keep_rownames <- rownames(my.responseorig)[which(rownames(my.responseorig) %in% rownames(keep_K))]
keep_K <- K[keep_rownames, keep_rownames]
keep_dat <- dat[which(rownames(my.responseorig) %in% rownames(keep_K)),]

# provide the kinship matrix 
gfit1 <- lmekin(MDS1 ~ PDSI + (1|id), data=keep_dat, varlist=keep_K)
fixef(gfit1)
ranef(gfit1) 
r.squaredLR(gfit1)
gfit2 <- lmekin(MDS1 ~ PDSI + (1|id), data=keep_dat)
r.squaredLR(gfit2)
gfit3 <- lmekin(MDS1 ~  (1|id), data=keep_dat, varlist=keep_K)
r.squaredLR(gfit3)

#################################
# ROC curves for cluster assignment
#################################
all_ROC = cbind(y_test, rf_test.output)
y_test = addLevel(all_ROC$y_test, newlevel="X4")
g1_3 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X1"))
g3_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X3"))
g1_2 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X1"))
g2_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X2"))
g2_3 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X2")) 
g3_2 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X3")) 

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
test_tog <- tog[-train_ind,] %>% filter(cluster%in%c("1", "3"))
train_tog <- tog[train_ind,] %>% filter(cluster%in%c("1", "3"))
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


PDSI_1 <- ggplot(data = newdat, aes(x=PDSI, y = vs)) +
  geom_line() +
  geom_jitter(data = tog, aes(x = PDSI, y = recode), cex = 2, alpha = 0.5, width = 0.01, height = 0.1) +
  theme_pubs +
  xlab("PDSI") + 
  ylab("Class 1 vs Class 3")

PDSI_2 <- ggplot(data = tog, aes(x = as.factor(recode), y = PDSI)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  xlab(c("Class")) +
  theme_pubs

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/PDSI_class.pdf", 
    useDingbats = FALSE, fonts = "ArialMT")
plot_grid(PDSI_1, PDSI_2)
dev.off()

#################################
# Step 5: Random forest across OTUs
################################





#################################
# calculate correlation matrix. Only possible with numeric predictors
#################################

nums <- unlist(lapply(my.total.matrix.num, is.numeric))  
correlationMatrix <- cor(my.total.matrix.num[,nums], use = "pairwise.complete.obs")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/env_heatmap.pdf",
    useDingbats = FALSE, fonts = "ArialMT")

heatmap.2(correlationMatrix, scale = "none", density.info="none", 
          trace="none", dendrogram = "none", col = Colors)

dev.off()
# summarize the correlation matrix
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.50)








# combine all data sets

col1 = x$Tour_ID#as.factor(sample_data(only_ath)$Clim)
col2 = fin_predictors$PDSI

p1 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col1), cex = 3) +
  scale_color_viridis_d() +
  theme_bw()

p2 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col2), cex = 3) +
  scale_colour_gradient2() +
  theme_bw()


plot_grid(import_MDS1, p2)

