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


run_rf <-function(fin_predictors,  metric = "RMSE"){
  set.seed(116)
  control <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats = 3,
                          verboseIter = TRUE,
                          savePredictions = TRUE)
  
  rf.output <- caret::train(response~., 
                            data = fin_predictors, 
                            method="rf", 
                            #importance = "permutation",
                            metric=metric, 
                            trControl = control,
                            verbose = TRUE)
  return(rf.output)
  
}

# Create a control object to set the details of feature selection that will occur next
feature_elim <- function(x, y, subsets = c(1:20,25, 30, 33)){
  set.seed(16)
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
  return(rfProfile)
}


#################################
# Step 1: Read in metadata and OTU data and filter to only plant
#################################
#metadata was edited and such in the prep_metadata.R script in the climate folder
#load("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metadata.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")

plant_val = which(OTU_clim$clim_data$Host_Species=="Ath")

cap_val = which(OTU_clim$clim_data$Host_Species!="Soil")

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

# plant_clim <- list(otu_table = OTU_clim$otu_table[plant_val,], 
#                    clim_data = OTU_clim$clim_data[plant_val,], 
#                    tax_table = OTU_clim$tax_table, 
#                    phy_tree = OTU_clim$phy_tree, 
#                    refseq = OTU_clim$refseq)

cap_clim <- list(otu_table = OTU_clim$otu_table[cap_val,], 
                 clim_data = OTU_clim$clim_data[cap_val,], 
                 tax_table = OTU_clim$tax_table, 
                 phy_tree = OTU_clim$phy_tree, 
                 refseq = OTU_clim$refseq)

#################################
# Step 2: Read in response variable 
#################################
#load(file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")

otu_name = "seq_1"
#my.response <- otu_table(GP_at15_all)[,otu_name]
#only_ath <- subset_samples(GP_at15_all, Subject %in% all_metadata$Sequence_ID)

  
#################################
# Match response variable and metadata then filter matrix and convert to factors and numerics
#################################

#MDS <- sqrt(sqrt(plant_clim$otu_table)) %>% dist() %>% cmdscale(eig = TRUE)

#my.responseorig <- MDS$points

my.var <- as.factor(plant_clim$clim_data$cluster) #my.responseorig[,1]

names(my.var) <- plant_clim$clim_data$Sequence_ID #rownames(my.responseorig)

data_frame_predictors <- plant_clim$clim_data #%>% select (-c(PlantID, Subject))

my.response1 <- my.var[match(data_frame_predictors$Plant_ID, names(my.var))] 

data_frame_predictors <- data_frame_predictors %>% select(-c(Plant_ID, Sequence_ID, cluster, Date, Tour_ID, Host_Species))

my.total.matrix <- cbind(data_frame_predictors, "otu" = my.response1)

my.total.matrix <- filter(my.total.matrix, is.na(otu) == FALSE) #%>% select (-c(Plant_ID, Sequence_ID))

my.total.matrix$Lat = as.numeric(as.character(my.total.matrix$Lat))

my.total.matrix$Long = as.numeric(as.character(my.total.matrix$Long))

my.total.matrix[which(is.na(my.total.matrix$Land_cov)),]$Land_cov = "5000"


#################################
# Step 3: Preprocess predictors (center and scale)
#################################
# Tutorial on RFE in caret: http://topepo.github.io/caret/recursive-feature-elimination.html

#center and scale the predictors
normalization <- preProcess(my.total.matrix %>% select(-otu),
                            method = c("knnImpute", "center", "scale"),
                            na.remove = TRUE)

x <- predict(normalization, my.total.matrix %>% select(-otu))
x <- as.data.frame(x)
x <- x[,colSums(is.na(x))==0,] %>% select(-c(Site_ID)) %>% select(-c(Site_name))

y <- as.factor(my.total.matrix$otu)

#################################
# Step 4: Choose features via recursive feature elimination
#################################

rfProfile <- feature_elim(x, y, subsets = c(1:20,25, 30, 33))

my.predictors = predictors(rfProfile)

#trellis.par.set(caretTheme())
#plot(rfProfile, type = c("g", "o"))

#################################
# Step 5: Esimate importance of features via random forest
#################################
# Subset predictors to those chosen in Step 1

fin_predictors <- x %>% select(my.predictors)

fin_predictors <- cbind(fin_predictors, response = y)

rf.output <-run_rf(fin_predictors, metric = "RMSE")

importance <- varImp(rf.output, scale=TRUE)

import_MDS1 <- plot(importance)

#################################
# Step 6: Plotting
#################################

# First the plot for the MDS on the total dataset
MDS.all <- (sqrt(OTU_clim$otu_table)) %>% dist() %>% cmdscale(eig = TRUE)
exp1 <-  ((MDS.all$eig) / sum(MDS.all$eig))[1]*100
exp2 <-  ((MDS.all$eig) / sum(MDS.all$eig))[2]*100

colnames(MDS.all$points) = c("MDS1", "MDS2")
col3 = OTU_clim$clim_data$Host_Species

all_data_MDS <- ggplot(data = data.frame(MDS.all$points), aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = col3), cex = 1, alpha = 0.8) +
  scale_color_viridis_d(labels = c(expression(italic("A. thaliana")), expression(italic("Capsella bursa-pastoris")), "Soil")) +
  xlab(paste(paste("MDS1 (", round(exp1), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp2), sep=""),"%)",sep="")) +
  theme_bw() +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.9),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )
        

MDS.cap <- (sqrt(cap_clim$otu_table)) %>% dist() %>% cmdscale(eig = TRUE)
col2 = cap_clim$clim_data$Host_Species
exp3 <-  ((MDS.cap$eig) / sum(MDS.cap$eig))[1]*100
exp4 <-  ((MDS.cap$eig) / sum(MDS.cap$eig))[2]*100
colnames(MDS.cap$points) = c("MDS1", "MDS2")

thaliana_cap_MDS <- ggplot(data = data.frame(MDS.cap$points), aes(x=MDS1, y=MDS2)) + 
  geom_point(cex = 1, alpha = 0.8, aes(col = col2)) +
  scale_color_manual(values = viridis_pal()(3)[c(1,2)], labels = c(expression(italic("A. thaliana")), expression(italic("Capsella bursa-pastoris")))) +
  xlab(paste(paste("MDS1 (", round(exp3), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp4), sep=""),"%)",sep="")) +
  theme_bw() +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.9),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )

legend <- get_legend(all_data_MDS) +
  theme(legend.position = "bottom")

all_capsella <- plot_grid(all_data_MDS + theme(legend.position = "none"), thaliana_cap_MDS + theme(legend.position = "none"))


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/MDS_cap_soil.pdf", useDingbats = FALSE, width = 3.5, height = 2)
plot_grid(all_capsella, legend, ncol = 1, rel_heights = c(1,0.1))
dev.off()

# Now Let's plot MDS with our random forest


my.plantID <- my.total.matrix %>% select(c(PDSI, Tour_ID, Plant_ID))

rownames(my.plantID) <-my.total.matrix$Plant_ID

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





# # calculate correlation matrix. Only possible with numeric predictors
# nums <- unlist(lapply(data_frame_predictors, is.numeric))  
# 
# correlationMatrix <- cor(data_frame_predictors[,nums], use = "pairwise.complete.obs")
# 
# heatmap.2(correlationMatrix, scale = "none", density.info="none", trace="none")
# 
# # summarize the correlation matrix
# # find attributes that are highly corrected (ideally >0.75)
# highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = 0.75)

#Notes: rfe is a simple backwards selection (recursive feature elimination algorith)

# problems with feature selection: https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-evalidation
# #need to get the contrast
#factor.vars <- my.total.matrix %>% select(c(TourID))
#contrasts <- lapply(my.total.matrix[, sapply(my.total.matrix, is.factor)], contrasts, contrasts = FALSE)
#dummy.model <- model.matrix(otu ~ ., data = my.total.matrix, contrasts.arg = contrasts)
