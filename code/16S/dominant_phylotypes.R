#!/usr/bin/env Rscript
library(phyloseq)
library(cowplot)
library(DESeq2)
# Goal of this script is to find the dominant phylotypes and write to file

####################################################################################
# Load the subsamples dataset
####################################################################################

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_ITS_GP1000.rds")
#subset to only athaliana
GP_athal <- subset_samples(GP1000, Species == "Ath")


####################################################################################
# Subset 16S to ASVs present in at least 5% of samples
####################################################################################
otus = as.matrix(otu_table(GP_athal))

is_zero = function(x){if(x==0)return(0)else return(1)}

is_there_otus <- apply(otus, 1:2, FUN = is_zero)

keep_otus <- which(colSums(is_there_otus)/dim(is_there_otus)[1] >= 0.05)

GP_at15 <- prune_taxa(names(keep_otus), GP_athal)

GP_Pseud <- subset_taxa(GP_at15, Genus=="Pseudomonas")

GP_at15_all = prune_taxa(colnames(otu_table(GP_at15)), GP1000)

save(GP_at15_all, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")


####################################################################################
# Subset ITS to ASVs present in at least 5% of samples
####################################################################################
otus = as.matrix(otu_table(GP_ITS_athal))

is_zero = function(x){if(x==0)return(0)else return(1)}

is_there_otus <- apply(otus, 1:2, FUN = is_zero)

keep_otus <- which(colSums(is_there_otus)/dim(is_there_otus)[1] >= 0.05)

GP_at15 <- prune_taxa(names(keep_otus), GP_athal)

GP_Pseud <- subset_taxa(GP_at15, Genus=="Pseudomonas")

GP_at15_all = prune_taxa(colnames(otu_table(GP_at15)), GP1000)

save(GP_at15_all, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_ITS_GP1000_at15.rds")
