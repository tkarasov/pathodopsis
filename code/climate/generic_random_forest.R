


#################################
# Functions for feature selection and random forest
#################################
generate_my_total_matrix <- function(response_choice, mat_clim){
  my.total.matrix <- cbind(mat_clim$clim_data, response = response_choice)#cbind(data_frame_predictors, "otu" = my.response1)
  my.plantID <- my.total.matrix %>% dplyr::select(c(PDSI, Tour_ID, Plant_ID))

  # There were some duplicate plants Let's remove from consideration
  dup <- which(!duplicated(my.plantID[,3]))
  my.total.matrix <- my.total.matrix[c(dup),]
  my.plantID <- my.plantID[c(dup),]

  rownames(my.plantID) <-my.total.matrix$Plant_ID
  my.total.matrix <- filter(my.total.matrix, is.na(response) == FALSE) %>% dplyr::select (-c(Plant_ID, Sequence_ID, Site_ID))
  #################################
  # Step 2: Preprocess predictor data into dummy variables
  #################################
  my.total.matrix$Land_cov <- replace_na(my.total.matrix$Land_cov, "5000")
  my.total.matrix[my.total.matrix=="n/a"]<-NA
  list_facs <- c("Lat", "Long", "Soil_temp", 
               "Soil_hum", "Humidity_ground.1", 
               "Air_temp", "Air_hum", "R_diameter")
  my.total.matrix[,list_facs] <- lapply(my.total.matrix[, list_facs], 
                                      function(x) as.numeric(as.character(x))) 

  #normalize the predictors and impute missing values. preProcess won't work with factors. 
  col_remove <- c("Albugo", "Tour_ID", "Necrosis", "Strata_trees", "Strata_shrubs", 
            "Strata_road", "Strata_wall_tree", "Strata_water", 
            "Date", "Strata_wall_shrub", "Land_cov", "Long", "Lat",
            "Site_Name", "Site_name", "Growing_on", 
            "Rosette_color", "Ground_type", "Heterogeneity", "cluster", "hc_cuttree2")
  my.total.matrix.num <- my.total.matrix %>% dplyr::select(-one_of(col_remove))

  return(my.total.matrix.num)
}
  
#################################
# Step 3: Preprocess data, Preprocess predictor data into dummy variables
#################################
preprocess_data <- function(my.total.matrix.num){
  normalization <- preProcess(my.total.matrix.num %>% dplyr::select(-c(response)),
                            method = c("center", "scale", "knnImpute"),
                            na.remove = TRUE)

  # Create a training set of 75% of the data
  smp_size <- floor(0.75 * nrow(my.total.matrix.num))
  set.seed(123)
  train_ind <- sample(seq_len(nrow(my.total.matrix.num)), 
    size = smp_size)

  # Now perform predictions
  x_full <- predict(normalization, my.total.matrix.num %>% dplyr::select(-c(response)))
  x_train <- as.data.frame(x_full[train_ind,])
  x_test <- as.data.frame(x_full[-train_ind,])
  #c("knnImpute", "center", "scale"),
  #normalization <- preProcess(my.total.matrix, method = c("center", "scale"), na.action = na.omit)
  x_train <- x_train[,colSums(is.na(x_train))==0,]   #%>% select(-c(Site_ID))
  x_test <- x_test[,colSums(is.na(x_test))==0,]

# Return data
  return(list(x_full, x_train, x_test, train_ind))
}

my_feat_selec <- function(x, y, subsets){
  # Create a control object to set the details of feature selection that will occur next
  set.seed(16)
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

my_random_forest<-function(my.predictors, x, y, classification=TRUE){
  # Subset predictors to those chosen in Step 1
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
  
  

# #Select features from full data
# subsets <- c(1:20,25, 30, 33)
# rfProfile = my_feat_selec(x_full, y, subsets)
# my.predictors = predictors(rfProfile)
# x_train = x_train[,my.predictors]
# x_test = x_test[,my.predictors]

# # Random Forest model (caret) on training data
# rf_train.output <-my_random_forest(my.predictors, x=x_train, y=y_train)
# rf_test.output<-predict(rf_train.output, newdata = x_test, type = "prob")

# importance <- varImp(rf_train.output$finalModel, scale=TRUE)
# import_class <- plot(importance)
# all_ROC = cbind(y_test, rf_test.output)

# y_test = addLevel(all_ROC$y_test, newlevel="X4")
# g1_3 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X1"))
# g3_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X3"), predictor = "X3"))
# g1_2 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X1"))
# g2_1 = (roc(all_ROC, response = "y_test", levels = c("X1","X2"), predictor = "X2"))
# g2_3 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X2")) 
# g3_2 = (roc(all_ROC, response = "y_test", levels = c("X2","X3"), predictor = "X3")) 


# full_comp <- ggroc(list("Class1vs2" = g1_2, "Class2vs1" = g2_1, 
#                         "Class2vs3" = g2_3, "Class3vs2" = g3_2,
#                         "Class1vs3" = g1_3, "Class3vs1" = g3_1), legacy.axes=TRUE) +
#   scale_color_viridis_d() +
#   geom_abline(intercept = 0, slope = 1,
#               color = "darkgrey", linetype = "dashed") +
#   theme_pubs +
#   xlab("False Positive Rate") +
#   ylab("True Positive Rate") +
#   theme(legend.position = c(0.8, 0.2),
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(override.aes = list(shape = 20)))

# pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/ASV_class_ROC.pdf", useDingbats = FALSE, fonts = "ArialMT")
# full_comp
# dev.off()

# auc(g1_3)
# #0.9343

# auc(g3_1)
# #0.9594

# auc(g1_2)
# #0.8269

# auc(g2_1)
# # 0.5898

# auc(g3_2)
# # 0.7385

# auc(g2_3)
# #0.6639

# # logistic regression of 1 vs 3
# tog <- data.frame(my.total.matrix)
# test_tog <- tog[-train_ind,] %>% filter(cluster%in%c("1", "3"))
# train_tog <- tog[train_ind,] %>% filter(cluster%in%c("1", "3"))
# tog <- tog %>% filter(cluster%in%c("1", "3"))
# tog$recode = c(0)
# tog[tog$cluster=="1",]$recode = 0
# tog[tog$cluster!="3",]$recode = 1
# mylogit.PDSI_vpd <- glm(recode ~ PDSI + vpd + srad + Elevation + tmax + def +  ws + Site_type + pet, data = tog, family = "binomial")
# mylogit.latitude <- glm(recode ~ Lat, data = tog, family = "binomial")
# newdat <- data.frame(PDSI=seq(min(tog$PDSI, na.rm = TRUE), max(tog$PDSI, na.rm = TRUE)),
#                               vpd = seq(min(tog$vpd, na.rm = TRUE), max(tog$vpd, na.rm= TRUE),
#                                         len=100))
# newdat$vs = predict(mylogit.PDSI_vpd , newdata=newdat, type="response")

# AIC(mylogit.PDSI_vpd, mylogit.latitude)

# logit_test <- glm(data = train_tog,
#                   recode ~ PDSI + vpd + srad ,
#                     family = "binomial")

# logit_Lat_test <- glm(data = train_tog, 
#                   recode ~ Lat,
#                   family = "binomial")


# test_tog$pred <- predict(logit_test, newdata=test_tog, type = "response")
# test_tog$pred_lat <- predict(logit_Lat_test, newdata=test_tog, type = "response")
# logit_roc <- roc(test_tog, response = cluster, predictor=pred)
# lat_roc <- roc(test_tog, response = cluster, predictor=pred_lat)

# pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/logit_RF_ROC_comparison.pdf")
# ggroc(list("Random Forest" = g1_3, "Full logit" = logit_roc, "logit Latitude Only" = lat_roc), legacy.axes = TRUE) +
#   scale_color_viridis_d() +
#   geom_abline(intercept = 0, slope = 1,
#               color = "darkgrey", linetype = "dashed") +
#   theme_pubs +
#   xlab("False Positive Rate") +
#   ylab("True Positive Rate") +
#   theme(legend.position = c(0.8, 0.2),
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(override.aes = list(shape = 20))) +
#   theme(legend.key=element_blank())

# dev.off()


# PDSI_1 <- ggplot(data = newdat, aes(x=PDSI, y = vs)) +
#   geom_line() +
#   geom_jitter(data = tog, aes(x = PDSI, y = recode), cex = 2, alpha = 0.5, width = 0.01, height = 0.1) +
#   theme_pubs +
#   xlab("PDSI") + 
#   ylab("Class 1 vs Class 3")

# PDSI_2 <- ggplot(data = tog, aes(x = as.factor(recode), y = PDSI)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.1, alpha = 0.5) +
#   xlab(c("Class")) +
#   theme_pubs

# pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/PDSI_class.pdf", 
#     useDingbats = FALSE, fonts = "ArialMT")
# plot_grid(PDSI_1, PDSI_2)
# dev.off()

# #################################
# # Step 5: Random forest across OTUs
# ################################




# # combine all data sets

# col1 = x$Tour_ID#as.factor(sample_data(only_ath)$Clim)
# col2 = fin_predictors$PDSI

# p1 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color = col1), cex = 3) +
#   scale_color_viridis_d() +
#   theme_bw()

# p2 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color = col2), cex = 3) +
#   scale_colour_gradient2() +
#   theme_bw()


# plot_grid(import_MDS1, p2)



# #Notes: rfe is a simple backwards selection (recursive feature elimination algorith)

# # problems with feature selection: https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-validation


# # #################################
# # # Step 1: Read in metadata
# # #################################
# # # load the metadata. How was this built?
# # metadata = read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_1_2020_reads.csv", header=T, fill =TRUE) %>% select(-c(X, X.1, X.2))
# # #choose only A. thaliana
# # 
# # all_metadata = filter(metadata, Host_Species == "Ath")
# # 
# # 
# # #convert factors to factors
# # factor_var = c("Albugo", "Necrosis", "Land_cov")
# # all_metadata[,factor_var] <- apply(all_metadata[,factor_var], 2, as.factor)
# # 
# # #remove weird columns
# # all_metadata <- all_metadata %>% select(-c("Time"))
# # 
# # #combine columns over averages
# # many <- c("tmax", "tmin", "vap", "ppt", "srad", "soil", "ws", "aet", "def", "PDSI", "vpd", "pet")
# # 
# # for(val in many){
# #   rel = colnames(all_metadata)[startsWith(colnames(all_metadata), val)]
# #   mean_val = rowMeans(all_metadata[,rel], na.rm =T)
# #   all_metadata <- all_metadata %>% select(-c(rel))
# #   all_metadata[,val] = mean_val
# #   
# # }
# 
# #remove samples for which missing more than 20 predictors
# all_metadata<- all_metadata[-which(apply(all_metadata, 1, function(x) sum(is.na(x)))>20),]
# data_frame_predictors <- all_metadata
# 
# #too many factors
# facs <- unlist(lapply(ictors, is.factor)) 
# 
# 
