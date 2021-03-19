#!/bin/bash
# this script is the overall analysis script to generate the many analyses
scripts=/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/


# Make the generic OTU table
### 16S
$scripts/amplicon_general/after_dada2_make_otutable_generic.R

# Important output files
# This next file is the samples with at least 1000 reads, limited to those ASVs with representation in at least 5% of samples
# /ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds

### ITS
$scripts/amplicon_general/after_dada2_make_otutable_generic_ITS.R
# /ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_ITS_GP1000_at15.rds

### Merge the OTU file with the metadata files to create OTU_clim
# Climate data and OTU data: /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds
$scripts/climate/prep_metadata.R

# Compare plant species with differential abundance using edgeR
$scripts/climate/hierarchical_multiple_testing.R

# Compare correlations in abundance between soil and hosts for all microbes and for pathogens
$scripts/climate/general_spatial_dist.R

# PCA on the subsetted microbiome and kmeans clustering into phylotypes
$scripts/climate/kmeans_clustering.R

# Kriging on the microbiome with OTU5 and 
$scripts/climate/kriging_OTU5_map.R

# Random forest on the microbiome of plants
$scripts/climate/random_forest_iterative_temp.R

# Random forest on the microbiome of soils
$scripts/climate/random_forest_iterative_soil.R

# Plant genotype comparisons to show differences in Fst of immune genes
$scripts/plant_genotype/IBS.R

# Plant health comparisons
$scripts/plant_genotype/plant_health.R

# Host genotype comparisons

# Immune gene comparisons
