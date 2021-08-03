# DEPRECATED SCRIPT


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
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/color_pal.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/theme_pubs.rds")

# #################################
# # Step 2: Read in response variable and bind to metadata
# #################################
# load(file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")

# #################################
# # Step 1: Read in metadata Feb. 2020
# #################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

my.responseorig <- data.frame(sqrt(plant_clim$otu_table)) %>% dist() %>% cmdscale(eig = F) %>% data.frame()

colnames(my.responseorig) = c("MDS1", "MDS2")

col1 = plant_clim$clim_data$Tour_ID

#################################
# Set response variable as MDS1 or as cluster
#################################
my.response1 <- my.responseorig$MDS1 #my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))] 

my.response2 <- make.names(as.factor(plant_clim$clim_data$cluster))

my.response3 <- make.names(as.factor(plant_clim$clim_data$MDS2))





my.plantID <- my.total.matrix %>% dplyr::select(c(PDSI, Tour_ID, Plant_ID))
rownames(my.plantID) <-my.total.matrix$Plant_ID

my.total.matrix <- filter(my.total.matrix, is.na(response) == FALSE) %>% dplyr::select (-c(Plant_ID, Sequence_ID, Site_ID))

#################################
# Step 2: Preprocess predictor data into dummy variables
#################################

#need to get the contrast
# Tutorial on RFE in caret: http://topepo.github.io/caret/recursive-feature-elimination.html

#center and scale the predictors
#my.total.matrix$Lat = as.numeric(as.character(my.total.matrix$Lat))
#my.total.matrix$Long = as.numeric(as.character(my.total.matrix$Long))
my.total.matrix[which(is.na(my.total.matrix$Land_cov)),]$Land_cov = "5000"
my.total.matrix[my.total.matrix=="n/a"]<-NA
list_facs <- c("Lat", "Long", "Soil_temp", 
               "Soil_hum", "Humidity_ground.1", 
               "Air_temp", "Air_hum", "R_diameter")

my.total.matrix[,list_facs] <- lapply(my.total.matrix[, list_facs], 
                                      function(x) as.numeric(as.character(x))) 

#normalize the predictors and impute missing values. preProcess won't work with factors. 
my.total.matrix.num <- my.total.matrix %>% 
  dplyr::select(-c(Albugo, Tour_ID, Necrosis, Strata_trees, Strata_shrubs, 
            Strata_road, Strata_wall_tree, Strata_water, 
            Lat, Long, Date, Strata_wall_shrub, 
            Site_Name, Site_name, Growing_on, 
            Rosette_color, Ground_type, Heterogeneity, cluster, hc_cuttree2))


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


#################################
# Functions for feature selection and random forest
#################################
my_feat_selec <- function(x, y, subsets){
  set.seed(16)
  # Create a control object to set the details of feature selection that will occur next
  ctrl = rfeControl(functions = rfFuncs,
                    #method = "metric",
                    repeats = 1,
                    saveDetails = TRUE,
                    verbose = TRUE)
#I was getting a weird error with lmfuncs. Not clear if the problem was me or a bug. The following error was also found in this blog:https://community.rstudio.com/t/rfe-error-logistic-classification-undefined-columns-selected-and-match-requires-vector-arguments/23988
  prettySeq <- function(x) paste("Resample", gsub(" ", "0", format(seq(along = x))), sep = "")
  rfProfile <- rfe(x = x, y = y,
                   sizes = subsets,
                   rfeControl = ctrl, #rfeControl(functions = caretFuncs))
                   na.action = na.omit)

  my.predictors = predictors(rfProfile)[1:2]

  
  trellis.par.set(caretTheme())
  plot(rfProfile, type = c("g", "o"))
  return(rfProfile)
}

# Subset predictors to those chosen in Step 1
my_random_forest<-function(my.predictors, x, y,classification=TRUE){
  set.seed(116)
  if(classification==TRUE){
    metric="Accuracy"
  }
  else{
    metric <- "RMSE"
  }
  fin_predictors <- x %>% select(my.predictors)
  fin_predictors <- cbind(fin_predictors, response = y)
  
  control <- trainControl(method="repeatedcv", 
                          number=3, 
                          repeats = 3,
                          verboseIter = TRUE,
                          classProbs = TRUE,
                          savePredictions = TRUE)
  
  rf.output <- caret::train(response~., 
                            data = fin_predictors, 
                            method="rf", 
                            metric=metric, 
                            trControl = control,
                            verbose = TRUE,
                            ntree = 50)
  
  return(rf.output)
}


my_pairwise_ROC<-function(my.predictors, x, y){
  fin_predictors <- x %>% select(my.predictors)  
  df = list()
  for(level in levels(y)){
    ynew = as.numeric(y)
    hm = which(y!=level)
    ynew[hm] <- 0
    ynew[!hm] <- 1
    ynew <- as.factor(ynew)
    fin_predictors2 <- cbind(fin_predictors, response = ynew)
    rf = my_random_forest(my.predictors, x=x, y=ynew)
    mtry = rf$bestTune
    rf.mytry = rf$pred %>% dplyr::filter(mtry == mtry)
    df1 = data.frame(truth = rf.mytry$obs, pred = rf.mytry$pred)
    df[[level]] = df1
    df[[level]]$"level" = level
  }

  return(df)
}
  
  
#################################
# Step 3: Preprocess data
#################################

normalization <- preProcess(my.total.matrix.num %>% dplyr::select(-c(cluster,response)),
                            method = c("center", "scale", "knnImpute"),
                            na.remove = TRUE)

smp_size <- floor(0.75 * nrow(my.total.matrix.num))
set.seed(123)
train_ind <- sample(seq_len(nrow(my.total.matrix.num)), size = smp_size)
x_full <- predict(normalization, my.total.matrix.num %>% dplyr::select(-c(cluster,response)))
x_train <- as.data.frame(x_full[train_ind,])
x_test <- as.data.frame(x_full[-train_ind,])
#c("knnImpute", "center", "scale"),
#normalization <- preProcess(my.total.matrix, method = c("center", "scale"), na.action = na.omit)
x_train <- x_train[,colSums(is.na(x_train))==0,]   #%>% select(-c(Site_ID))
x_test <- x_test[,colSums(is.na(x_test))==0,]

#################################
# Step 4: Random forest on classification to cluster
################################
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

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/ASV_class_ROC.pdf", useDingbats = FALSE, fonts = "ArialMT")
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



#Notes: rfe is a simple backwards selection (recursive feature elimination algorith)

# problems with feature selection: https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-validation


# #################################
# # Step 1: Read in metadata
# #################################
# # load the metadata. How was this built?
# metadata = read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_1_2020_reads.csv", header=T, fill =TRUE) %>% select(-c(X, X.1, X.2))
# #choose only A. thaliana
# 
# all_metadata = filter(metadata, Host_Species == "Ath")
# 
# 
# #convert factors to factors
# factor_var = c("Albugo", "Necrosis", "Land_cov")
# all_metadata[,factor_var] <- apply(all_metadata[,factor_var], 2, as.factor)
# 
# #remove weird columns
# all_metadata <- all_metadata %>% select(-c("Time"))
# 
# #combine columns over averages
# many <- c("tmax", "tmin", "vap", "ppt", "srad", "soil", "ws", "aet", "def", "PDSI", "vpd", "pet")
# 
# for(val in many){
#   rel = colnames(all_metadata)[startsWith(colnames(all_metadata), val)]
#   mean_val = rowMeans(all_metadata[,rel], na.rm =T)
#   all_metadata <- all_metadata %>% select(-c(rel))
#   all_metadata[,val] = mean_val
#   
# }
# 
# #remove samples for which missing more than 20 predictors
# all_metadata<- all_metadata[-which(apply(all_metadata, 1, function(x) sum(is.na(x)))>20),]
# data_frame_predictors <- all_metadata
# 
# #too many factors
# facs <- unlist(lapply(ictors, is.factor)) 
# 
# 
