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


#################################
# Step 1: Read in metadata
#################################
# load the data
metadata = read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_1_2020_reads.csv", header=T, fill =TRUE) %>% select(-c(X, X.1, X.2))

#choose only A. thaliana
all_metadata = filter(metadata, Host_Species == "Ath")


#convert factors to factors
factor_var = c("Albugo", "Necrosis", "Land_cov")
all_metadata[,factor_var] <- apply(all_metadata[,factor_var], 2, as.factor)

#remove weird columns
all_metadata <- all_metadata %>% select(-c("Time"))

#combine columns over averages
many <- c("tmax", "tmin", "vap", "ppt", "srad", "soil", "ws", "aet", "def", "PDSI", "vpd", "pet")

for(val in many){
  rel = colnames(all_metadata)[startsWith(colnames(all_metadata), val)]
  mean_val = rowMeans(all_metadata[,rel], na.rm =T)
  all_metadata <- all_metadata %>% select(-c(rel))
  all_metadata[,val] = mean_val
  
}

#remove samples for which missing more than 20 predictors
all_metadata<- all_metadata[-which(apply(all_metadata, 1, function(x) sum(is.na(x)))>20),]
data_frame_predictors <- all_metadata

#too many factors
facs <- unlist(lapply(ictors, is.factor)) 


#################################
# Step 2: Read in response variable and bind to metadata
#################################
load(file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")
otu_name = "seq_1"
#my.response <- otu_table(GP_at15_all)[,otu_name]
only_ath <- subset_samples(GP_at15_all, Subject %in% all_metadata$Sequence_ID)
my.responseorig <- data.frame(sqrt(sqrt(otu_table(only_ath)+1)) %>% dist() %>% cmdscale())
#colnames(my.response) <- "otu"
colnames(my.responseorig) = c("MDS1", "MDS2")
col1 = as.factor(sample_data(only_ath)$TourID)

ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col1), cex = 3) +
  scale_color_viridis_d() +
  theme_bw()
  
#################################
# Set response variable as MDS1
#################################
my.response1 <- my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))] 
#my.response1 <- my.response[match(data_frame_predictors$Sequence_ID, rownames(my.response))] 
#data_frame_predictors <- sample_data(GP_at15_all) %>% select (-c(PlantID, Subject))
my.total.matrix <- cbind(data_frame_predictors, "otu" = my.response1)
my.plantID <- my.total.matrix %>% select(c(PDSI, Tour_ID, Plant_ID))
rownames(my.plantID) <-my.total.matrix$Plant_ID
my.total.matrix <- filter(my.total.matrix, is.na(otu) == FALSE) %>% select (-c(Plant_ID, Sequence_ID))

# calculate correlation matrix. Only possible with numeric predictors
nums <- unlist(lapply(data_frame_predictors, is.numeric))  
correlationMatrix <- cor(data_frame_predictors[,nums], use = "pairwise.complete.obs")
heatmap.2(correlationMatrix, scale = "none", density.info="none", trace="none")

# summarize the correlation matrix
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)

#################################
# Step 2: Preprocess predictor data into dummy variables
#################################

#need to get the contrast
#factor.vars <- my.total.matrix %>% select(c(TourID))
#contrasts <- lapply(my.total.matrix[, sapply(my.total.matrix, is.factor)], contrasts, contrasts = FALSE)
#dummy.model <- model.matrix(otu ~ ., data = my.total.matrix, contrasts.arg = contrasts)

# Tutorial on RFE in caret: http://topepo.github.io/caret/recursive-feature-elimination.html

#center and scale the predictors
my.total.matrix$Lat = as.numeric(as.character(my.total.matrix$Lat))
my.total.matrix$Long = as.numeric(as.character(my.total.matrix$Long))
my.total.matrix[which(is.na(my.total.matrix$Land_cov)),]$Land_cov = "5000"

normalization <- preProcess(my.total.matrix %>% select(-otu),
                            method = c("knnImpute", "center", "scale"),
                            na.remove = TRUE)

#normalization <- preProcess(my.total.matrix, method = c("center", "scale"), na.action = na.omit)
subsets <- c(1:20,25, 30, 33)

x <- predict(normalization, my.total.matrix %>% select(-otu))

x <- as.data.frame(x)

y <- my.total.matrix$otu

#rm NA predictors
x <- x[,colSums(is.na(x))==0,] %>% select(-c(Site_ID)) %>% select(-c(Site_name))


#################################
# Step 3: Choose features via recursive feature elimination
#################################

set.seed(16)

# Create a control object to set the details of feature selection that will occur next
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   saveDetails = TRUE,
                   verbose = TRUE)
#I was getting a weird error with lmfuncs. Not clear if the problem was me or a bug. The following error was also found in this blog:https://community.rstudio.com/t/rfe-error-logistic-classification-undefined-columns-selected-and-match-requires-vector-arguments/23988

rfProfile <- rfe(x = x, y = y,
                 sizes = subsets,
                 rfeControl = ctrl, #rfeControl(functions = caretFuncs))
                 na.action = na.omit)

my.predictors = predictors(rfProfile)

trellis.par.set(caretTheme())
plot(rfProfile, type = c("g", "o"))


#################################
# Step 4: Esimate importance of features via random forest
#################################
# Subset predictors to those chosen in Step 1
set.seed(116)
metric <- "RMSE"
fin_predictors <- x %>% select(my.predictors)
fin_predictors <- cbind(fin_predictors, response = y)

control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats = 3,
                        verboseIter = TRUE,
                        savePredictions = TRUE)

rf.output <- caret::train(response~., 
                   data = fin_predictors, 
                   method="rf", 
                
                   importance = "permutation",
                   metric=metric, 
                   trControl = control,
                   verbose = TRUE)
                   

importance <- varImp(rf.output, scale=FALSE)
import_MDS1 <- plot(importance)

# combine all data sets

col1 = as.factor(sample_data(only_ath)$Clim)
col2 = my.plantID[rownames(my.responseorig)]

p1 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col1), cex = 3) +
  scale_color_viridis_d() +
  theme_bw()

p2 <- ggplot(data = my.responseorig, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col2), cex = 3) +
  scale_colour_gradient2() +
  theme_bw()


plot_grid(import_MDS1, p2)

####
#Now do the same thing as above but with MDS2




#Notes: rfe is a simple backwards selection (recursive feature elimination algorith)

# problems with feature selection: https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-validation