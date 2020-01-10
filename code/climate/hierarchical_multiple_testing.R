#!/usr/bin/env Rscript

library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)
library(ade4)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(genefilter)
library(DESeq2)
library(BiocParallel)

#the goal of this script is to determine OTUs differentially abundant between  (1) A. thalaina and Capsella and differentially abundant between (2) Environmental Conditions.



#Much of the analysis comes from here: https://f1000research.com/articles/5-1492/v2

#######################################################################
# Read in Phyloseq object
#######################################################################
# Load the phyloseq object generated in after_dada2_make_otutable
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")

# #######################################################################
# # Prune phyloseq to only those samples with at least 1000 reads, and to those ASVs that are zero add 1
# #######################################################################
GP_fam = tax_glom(GP50, "Family")

#Prune samples with unknown family
keep = tax_table(GP_fam)[tax_table(GP_fam)[,"Family"]!="Unknown_Family",]
GP_fam = prune_taxa(rownames(keep), GP_fam)

#rename the taxa with family
taxa_names(GP_fam) = tax_table(GP_fam)[,"Family"]

#add 1 to every sample to enable geometric mean calculation in the variance stabilization
otu_table(GP_fam) = otu_table(GP_fam) + 1


#######################################################################
# Convert Phyloseq to deseq
#######################################################################
ps_dds <- phyloseq_to_deseq2(GP_fam, ~ PlantID + Sample_type)

#######################################################################
# Perform variance stabilizing transformation
#######################################################################
vsd <- varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "local")


#######################################################################
# Perform variance stabilizing transformation
#######################################################################
ps_dds <- estimateSizeFactors(ps_dds)
ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)

# perform DEA
dea_indp <- DESeq(dds_indp)
dea_paired <- DESeq(ps_dds, parallel=TRUE, BPPARAM=MulticoreParam(8), quiet = FALSE)

# result analysis
res_indp <- results(dea_indp)
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/dea_unpaired.rds")
res <- results(dea_paired)


t(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

# p-adjusted value
padj_indp <- res_indp$padj
padj_paired <- res_paired$padj

# gene list
geneList <- rownames(control)

# significantly expressed genes
sig_indp <- geneList[which(padj_indp < 0.05)]
sig_paired <- geneList[which(padj_paired < 0.05)]


