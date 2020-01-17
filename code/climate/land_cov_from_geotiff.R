#!/usr/bin/env Rscript

#the goal of this script is to take the land coverage info from a geotiff

library(rgdal)
library(raster)
library(dplyr)

####################################################################################
# Load the geotiff and metadata into memory
####################################################################################

land_cov <- raster( x = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/climate_spatial_modelling/data/Globcover2009_V2.3_Global_/GLOBCOVER_L4_200901_200912_V2.3.tif")
land_cov_df <- as.data.frame(land_cov, xy = TRUE)

relevant_crs = crs(land_cov)

#my metadata
metadata0 = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course3.txt", header=T, sep="\t")
metadata = filter(metadata0, is.na(Lat)==F & is.na(Long)==F)
coordinates(metadata) = ~Long + Lat


####################################################################################
# Extract points
####################################################################################

land_value <- extract(land_cov, metadata)
metadata$Land_cov = land_value



####################################################################################
# Rebuild metadata original and write to file
####################################################################################
rownames(metadata0) = metadata0$Sequence_ID
metadata0$Land_cov[metadata$Sequence_ID] = metadata$Land_cov

write.table(metadata0, "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course3.txt", col.names = T, sep="\t", row.names = T, quote = F)
-