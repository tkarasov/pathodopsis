#!/usr/bin/env Rscript


# The goal of this script is to take any variable and the given climate variables and determine their relationship
library(randomForest)
library(dplyr)

#######################################################################
# Read in Climate data
#######################################################################

#terraclim data from Grey
# tmax:Maximum temperature, 
# tmin: minimum temperature, 
# vp:vapor pressure, 
# ppt:precipitation accumulation, 
# srad:downward surface shortwave radiation, 
# ws: wind-speed
# pet:Reference evapotranspiration (ASCE Penman-Montieth), 
# q:Runoff*, 
# aet:Actual Evapotranspiration*, 
# def:Climate Water Deficit*, 
# soil:Soil Moisture*, 
# swe:Snow Water Equivalent*, 
# PDSI:Palmer Drought Severity Index, 
# vpd:Vapor pressure deficit

clim_gen = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/terraclim_1_2020.csv", sep = " ", header = T)

#my metadata
metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course3.txt", header=T, sep="\t")

# Let's do load as a trial. First 14 columns are variables
met_load = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_9_2018_reads.txt", sep ="\t", header = T)
met_load$Total_load = rowSums(met_load[,c(15:dim(met_load)[2])], na.rm = T)

#merge all metadata
all_metadata = metadata %>% full_join(clim_gen, by = c("Site_ID", "Tour_ID")) %>% full_join(met_load, by.x = "Site_ID", by.y = "Site_ID")

#Make variable that is only first letter of Koppen-Geiger
all_metadata$First_Kop = substr(all_metadata$ClimateZ, 1,1)

#choose only A. thaliana
all_metadata = filter(all_metadata, Host_Species == "Ath")

#######################################################################
# Read in variable of interest
#######################################################################
my_var = "Total_load"

# Limit to only non-soil
data1 = filter(all_metadata, Used_for == "M")
is_na = is.na(data1[,my_var])
#data1 = data1[!is_na,]

# Relevant variables
rel = c("Tour_ID", "Latitude", "Longitude", "vpd.6", "ppt.6", "PDSI.6", "ClimateZ", "soil.6", "pet.6", "def.6", "tmax.6", "tmin.6", "First_Kop", "aet.6", "srad.6", my_var) #, my_var, "Peronosporaceae",
data1 = data1[,rel]
#remove samples with any NAs
data1 = data1[complete.cases(data1),]

#Make all variables factors if characters
data1_factor <- data1 %>%
  mutate_if(sapply(data1, is.character), as.factor)

data1 = data1_factor
#data1$Latitude = as.numeric(as.character(data1$Lat.x))
#data1$Lat.y = as.numeric(as.character(data1$Lat.y))

#######################################################################
# Perform random forest analysis of variables
#######################################################################
# When I performed random forests with some variables I got negative variance explained. 
# http://developmentaldatascience.org/post/29-01-18_metaforest_no_effect/


# Split data into training and validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(116)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]
ValidSet <- data1[-train,]
summary(TrainSet)
summary(ValidSet)

# Create a Random Forest model with default parameters
my.formula = as.formula(paste(paste(paste("log10(", my_var, sep=""), ")", sep = ""),"~ .", sep = ""))
model1 <- randomForest(my.formula, data = TrainSet, importance = TRUE, ntree = 1000)
model1
imp = importance(model1)
varImpPlot(model1)

# Predicting on train set
predTrain <- predict(model1, TrainSet, type = "class")

# Checking classification accuracy
table(predTrain, TrainSet[,my_var])  

# Predicting on Validation set
predValid <- predict(model1, ValidSet, type = "class")

# Checking classification accuracy
mean(predValid == ValidSet$Condition)                    
table(predValid,ValidSet$Peronosporaceae)

# # Using For loop to identify the right mtry for model
# a=c()
# i=5
# for (i in 3:8) {
#   model3 <- randomForest(my.formula, data = TrainSet, ntree = 500, mtry = i, importance = TRUE)
#   predValid <- predict(model3, ValidSet, type = "class")
#   a[i-2] = mean(predValid == ValidSet$Condition)
#   print(model3)
# }
# 
# a
# 
# plot(3:8,a)

mymodel.rfP <- rfPermute(
  my.formula, data = TrainSet, ntree = 1000, 
  na.action = na.omit, nrep = 100, num.cores = 1
)

# Plot the null distributions and observed values.
plotNull(mymodel.rfP) 

# Plot the unscaled importance distributions and highlight significant predictors
plot(rp.importance(mymodel.rfP, scale = FALSE))





