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
library(edgeR)

#the goal of this script is to determine OTUs differentially abundant between  (1) A. thalaina and Capsella and differentially abundant between (2) Environmental Conditions.

#Much of the analysis comes from here: https://f1000research.com/articles/5-1492/v2

#######################################################################
# Read in Phyloseq object
#######################################################################
# Load the phyloseq object generated in after_dada2_make_otutable
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")

# #######################################################################
# # Prune phyloseq to only those samples with at least 1000 reads, and to those ASVs that are zero add 1
# #######################################################################
# GP_fam = tax_glom(GP50, "Family")
# 
# #Prune samples with unknown family
# keep = tax_table(GP_fam)[tax_table(GP_fam)[,"Family"]!="Unknown_Family",]
# GP_fam = prune_taxa(rownames(keep), GP_fam)
# 
# #rename the taxa with family
# taxa_names(GP_fam) = tax_table(GP_fam)[,"Family"]

#add 1 to every sample to enable geometric mean calculation in the variance stabilization
otu_table(GP_at15_all) = otu_table(GP_at15_all) + 1

# #######################################################################
# For differential analysis use limma-voom
# #######################################################################
# https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html
m = as(otu_table(GP_at15_all), "matrix")
# Add one to protect against overflow, log(0) issues.
m = m + 1
# Define gene annotations (`genes`) as tax_table
taxonomy = tax_table(GP_at15_all, errorIfNULL=FALSE)
if( !is.null(taxonomy) ){
  taxonomy = data.frame(as(taxonomy, "matrix"))
} 
# Now turn into a DGEList
d = DGEList(counts=m, genes=taxonomy, remove.zeros = TRUE)

# Calculate the normalization factors
z = calcNormFactors(d, method="RLE")
# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(z$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}






# #######################################################################
# # Convert Phyloseq to deseq
# #######################################################################
# ps_dds <- phyloseq_to_deseq2(GP_at15_all, ~ PlantID + Sample_type)
# 
# #######################################################################
# # Perform variance stabilizing transformation
# #######################################################################
# vsd <- varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "local")
# 
# 
# #######################################################################
# # Perform variance stabilizing transformation
# #######################################################################
# ps_dds <- estimateSizeFactors(ps_dds)
# ps_dds <- estimateDispersions(ps_dds)
# abund <- getVarianceStabilizedData(ps_dds)
# 
# # perform DEA
# dea_indp <- DESeq(dds_indp)
# dea_paired <- DESeq(ps_dds, parallel=TRUE, BPPARAM=MulticoreParam(8), quiet = FALSE)
# 
# # result analysis
# res_indp <- results(dea_indp)
# load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/dea_unpaired.rds")
# res <- results(dea_paired)
# 
# 
# t(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
# 
# # p-adjusted value
# padj_indp <- res_indp$padj
# padj_paired <- res_paired$padj
# 
# # gene list
# geneList <- rownames(control)

# significantly expressed genes
sig_indp <- geneList[which(padj_indp < 0.05)]
sig_paired <- geneList[which(padj_paired < 0.05)]


















####################################################################################
# Let's look at the  abundance of the microbes across environments
####################################################################################

GP_Pseud <- subset_taxa(GP_at15_all, Genus=="Pseudomonas")

GP_Sphing <- subset_taxa(GP_at15_all, Genus == "Sphingomonas")

GP_Burkholderia <- subset_taxa(GP_at15_all, Family == "Burkholderiaceae")

####################################################################################
# Show phylogeny 
####################################################################################
p.pseud = plot_tree(GP_Pseud, 
                    size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species")+
  scale_color_viridis_d() + ggtitle("Pseudomonas")

#+ coord_polar(theta="y")

p.sping = plot_tree(GP_Sphing, 
                    size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species") +
  scale_color_viridis_d() + ggtitle("Sphingomonas")

p.burk = plot_tree(GP_Burkholderia,
                   size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species")+
  scale_color_viridis_d() + ggtitle("Burkholderiaceae")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/genus_level_phyl_abund.pdf", uuseDingbats = FALSE, fonts = "ArialMT", width = 12, height =12)
plot_grid(p.pseud, p.sping, p.burk, nrow = 1)
dev.off()


####################################################################################
# Compare abundance between environment
####################################################################################

