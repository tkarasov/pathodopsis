#!/usr/bin/env Rscript

library(phyloseq)
library(dada2)
library(dplyr)
library(tidyverse)
library(genefilter)
library(DECIPHER)
library(microbiome)
library(phangorn)

#This script is supposed to take the OTUs from dada2 and create a tree. The tutorial from which this code is pilfered can be found here: https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html


path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/amp_seq_functions.R", sep=""))
output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"
#output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/"

#######################################################################
# Read in Phyloseq object
#######################################################################
# seqtab.nochim = readRDS(paste(output_direc,"/seqtab_final.rds", sep="/"))
# taxa = readRDS(paste(output_direc,"/tax_final.rds", sep="/"))

# Load the phyloseq object generated in after_dada2_make_otutable
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")
sequences<-taxa_names(GP50)
names(sequences)<-sequences

#######################################################################
# Read in Phyloseq object
#######################################################################
print("Now starting alignment")

#This next step of aligning takes A LOT OF MEMORY. I had to reserve 8CPUs at 64Gb RAM per CPU
alignment <- AlignSeqs(refseq(GP50), anchor=NA, processors = 8, verbose = TRUE) 
  #AlignSeqs(DNAStringSet(sequences), anchor=NA, processors = NULL, verbose = TRUE)
save(alignment, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU_alignment.rds")

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phang.align)

treeNJ <- NJ(dm) # Note, tip order != sequence order

fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#######################################################################
# Save tree to image
#######################################################################
save(figGTR, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU_tree.RData")