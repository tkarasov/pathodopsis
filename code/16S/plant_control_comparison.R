library(phyloseq)
library(dplyr)
library(ggplot2)
library(dada2)
library(tidyr)
library(DESeq2)
library(reshape2)
library(cowplot)
library(lme4)
library(lmerTest)

# The goal of this script is to look for plate effect.
####################################
#Load data
####################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
demultiplexing_file=read.table(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/plate_location_map_16S_ITS1.txt",
  sep = "\t", header = TRUE)

# add MDS to plant_clim
my.responseorig <- data.frame(sqrt(plant_clim$otu_table)) %>% dist() %>% cmdscale(eig = F) %>% data.frame()
colnames(my.responseorig) = c("MDS1", "MDS2")
plant_clim$clim_data$MDS1 = my.responseorig$MDS1
plant_clim$clim_data$MDS2 = my.responseorig$MDS2

####################################
# Add sequencing information to plate
####################################
tog <- merge(plant_clim$clim_data, demultiplexing_file, 
             by.x = "Plant_ID",
             by.y = "ample_ID",
             all = TRUE)



####################################
# Let's do ANOVA now to look at the relative contribution of plate vs. Tour
####################################
my.responseorig <- data.frame(sqrt(plant_clim$otu_table)) %>% dist() %>% cmdscale(eig = F) %>% data.frame()

colnames(my.responseorig) = c("MDS1", "MDS2")

anova(MDS1 ~ Tour_ID + Sample_Plate, data = tog)

# Now let's see how plate identity influences the loading on MDS1
model1 <- lmer(MDS1 ~ (1|Tour_ID) + (1|Sample_Plate) + Lat, data = tog)
summary(model1)

# Not significant. Great!

####################################
# How many reads in the controls?
####################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP_nofilter.rds")