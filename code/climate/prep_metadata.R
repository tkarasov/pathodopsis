#!/usr/bin/env Rscript

# This is the key script for merging the metadata and the phyloseq object


library(dplyr)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(phyloseq)

#the goal of this script is to prep the metadata for random forest

metadata = read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_1_2020_reads.csv", header=T, fill =TRUE) %>% select(-c(X, X.1, X.2))

#choose only A. thaliana
all_metadata = metadata #filter(metadata, Host_Species == "Ath")

#rename Host_species column and remove used_for
all_metadata$Host_Species = as.character(all_metadata$Host_Species)
all_metadata[all_metadata$Used_for == "SOIL",]$Host_Species = "Soil"
all_metadata$Host_Species = as.factor(all_metadata$Host_Species)
all_metadata = all_metadata %>% select(-c(Used_for))

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
facs <- unlist(lapply(data_frame_predictors, is.factor)) 

save(all_metadata, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metadata.rds")


#################################
# Step 1: Read in metadata
#################################
#metadata was edited and such in the prep_metadata.R script in the climate folder
load("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metadata.rds")

#################################
# Step 2: Read in response variables
#################################
load(file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")


#################################
# Step 3: Reorder climate variables according to OTU table
#################################
phy_seq_reorder <- subset_samples(GP_at15_all, Subject %in% all_metadata$Sequence_ID)
all_metadata_reorder <- all_metadata[match(sample_data(phy_seq_reorder)$Subject, all_metadata$Sequence_ID),] 

# otu_name = "seq_1"
# #my.response <- otu_table(GP_at15_all)[,otu_name]
# only_ath <- subset_samples(GP_at15_all, Subject %in% all_metadata$Sequence_ID)


#################################
# Create S3 object
#################################

OTU_clim <- list(otu_table = otu_table(phy_seq_reorder), 
                 clim_data = all_metadata_reorder, 
                 tax_table = tax_table(phy_seq_reorder), 
                 phy_tree = phy_tree(phy_seq_reorder), 
                 refseq = refseq(phy_seq_reorder))

save(OTU_clim, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
