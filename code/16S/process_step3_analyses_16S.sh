#!/bin/bash
# this script is the overall analysis script to generate the many analyses
scripts=/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/

# Make the generic OTU table
# 16S
$scripts/amplicon_general/after_dada2_make_otutable_generic.R

# ITS
$scripts/amplicon_general/after_dada2_make_otutable_generic_ITS.R

# Compare plant species
$scripts/climate/hierarchical_multiple_testing.R

# PCA on the subsetted microbiome

# Kriging on the microbiome
$scripts/climate/kmeans_clustering.R

# Random forest on the microbiome
$scripts/climate/random_forest_iterative.R

# Plant genotype comparisons to show differences in Fst of immune genes
$scripts/plant_genotype/IBS.R

# Plant health comparisons
$scripts/plant_genotype/plant_health.R
