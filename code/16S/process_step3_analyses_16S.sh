#!/bin/bash

scripts=/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/scripts/
# Make the generic OTU table
after_dada2_make_otutable_generic.R

# compare plant species
$scripts/hierarchical_multiple_testing.R

# PCA on the subsetted microbiome



# Kriging on the microbiome
$scripts/climate/kmeans_clustering.R

# Random forest on the microbiome
$scripts/climate/random_forest_iterative.R
