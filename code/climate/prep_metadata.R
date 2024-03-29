#!/usr/bin/env Rscript

# This is the key script for merging the metadata and the phyloseq object

library(devtools)
library(dplyr)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(phyloseq)
library(intrval)

#the goal of this script is to prep the metadata for random forest

metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_2_2020_reads.tsv", header=T, fill =TRUE, sep = "\t") #%>% select(-c(X, X.1, X.2))

#filter out the weird rows that don't start with P
metadata2 = metadata[which(grepl("^P",metadata$Sequence_ID)),]


#choose only A. thaliana
all_metadata = metadata2 #filter(metadata, Host_Species == "Ath")

#rename Host_species column and remove used_for
all_metadata$Host_Species = as.character(all_metadata$Host_Species)
all_metadata[which(all_metadata$Used_for == "SOIL"),]$Host_Species = "Soil"
all_metadata$Host_Species = as.factor(all_metadata$Host_Species)
all_metadata = all_metadata %>% select(-c(Used_for))

#convert factors to factors
factor_var = c("Albugo", "Necrosis", "Land_cov")
all_metadata[,factor_var] <- apply(all_metadata[,factor_var], 2, as.factor)

#remove weird columns
all_metadata <- all_metadata %>% select(-c("Time", "Humidity_ground.1", "Sun.1", "Slope.1", "Site_Name"))

#combine columns over averages
many <- c("tmax", "tmin", "vap", "ppt", "srad", "soil", "ws", "aet", "def", "PDSI", "vpd", "pet")

for(val in many){
  rel = colnames(all_metadata)[startsWith(colnames(all_metadata), val)]
  rel_mat <-data.frame(lapply(all_metadata[,rel], function(x) as.numeric(as.character(x))))
  mean_val = rowMeans(rel_mat, na.rm =T)
  all_metadata <- all_metadata %>% select(-c(rel))
  all_metadata[,val] = mean_val
}

#remove samples for which missing more than 20 predictors
all_metadata<- all_metadata[-which(apply(all_metadata, 1, function(x) sum(is.na(x)))>20),]
data_frame_predictors <- all_metadata

#too many factors
facs <- unlist(lapply(data_frame_predictors, is.factor)) 

#convert factors to factors and other to other
factor_table <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Meta_table_from_rebecca_2_2020 - variable_types.tsv", header = T, sep = "\t", fill = T, row.names = 1)
is_factor <- colnames(factor_table)[which(factor_table["Type",]=="Qualitative")]
not_factor <- colnames(factor_table)[which(factor_table["Type",]=="Quantitative")]
all_metadata[,is_factor] <- lapply(all_metadata[,is_factor], function(x) as.factor(x))
all_metadata[,not_factor] <- lapply(all_metadata[,not_factor], function(x) as.numeric(as.character(x)))
all_metadata$Lat <- as.numeric(as.character(all_metadata$Lat))

to_exclude <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/exclusion_notes.tsv.txt", sep="\t", header = TRUE)

remove_all <- to_exclude %>% filter(samples=="all")

# Two samples had repeatedly weird input problems: PA1336.P207 and PA1339.P207 so were excluded and PA1408, PA1410 were amended because of problems
all_metadata <- all_metadata %>% filter(Sequence_ID %ni% c("PA1336.P207.IRE", "PA1339.P207.IRE"))

remove_PC <- to_exclude %>% filter(what_to_exclude == "PC")
all_metadata <- all_metadata %>% filter(Tour_ID %ni% remove_all$TOUR_ID)
all_metadata <- all_metadata %>% filter(Sequence_ID %ni% remove_PC$samples)

write.table(all_metadata, "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_9_2020_reads.tsv", col.names=TRUE,  sep = "\t", quote=FALSE, row.names=FALSE)

save(all_metadata, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metadata.rds")


#################################
# Step 1: Read in metadata
#################################
#metadata was edited and such in the prep_metadata.R script in the climate folder
load("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metadata.rds")

#################################
# Step 2: Read in response variables
#################################
# Load the reduced dataset
load(file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")

# Also load the original GP50 for downstream comparative analysis
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")

#################################
# Step 3: Reorder climate variables according to OTU table
#################################
phy_seq_reorder <- subset_samples(GP_at15_all, samples.out %in% all_metadata$Sequence_ID)
all_metadata_reorder <- all_metadata[match(sample_data(phy_seq_reorder)$samples.out, all_metadata$Sequence_ID),] 

#################################
# Step 4: Remake GP50 with amended metadata and reduced samples
#################################
phy_seq50_reorder <- subset_samples(GP50, samples.out %in% all_metadata$Sequence_ID)
GP50_all_metadata_reorder <- all_metadata[match(sample_data(phy_seq50_reorder)$samples.out, all_metadata$Sequence_ID),]
rownames(GP50_all_metadata_reorder) <- GP50_all_metadata_reorder$Sequence_ID
sample_data(phy_seq50_reorder) <- sample_data(GP50_all_metadata_reorder)
save(phy_seq50_reorder, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/GP50_subset_fin.rds")

#################################
# Create S3 object OTU_clim
#################################

OTU_clim <- list(otu_table = otu_table(phy_seq_reorder), 
                 clim_data = all_metadata_reorder, 
                 tax_table = tax_table(phy_seq_reorder), 
                 phy_tree = phy_tree(phy_seq_reorder), 
                 refseq = refseq(phy_seq_reorder))

save(OTU_clim, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")


#################################
# Create S3 object plant_clim
#################################

plant_val = which(OTU_clim$clim_data$Host_Species=="Ath")

plant_clim <- list(otu_table = OTU_clim$otu_table[plant_val,], 
                   clim_data = OTU_clim$clim_data[plant_val,], 
                   tax_table = OTU_clim$tax_table, 
                   phy_tree = OTU_clim$phy_tree, 
                   refseq = OTU_clim$refseq)

save(plant_clim, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

