library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)
library(ade4)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
#library(genefilter)
library(limma)
library(edgeR)

# The goal of this script is to associate microbial similarity with distance with genetic similarity
# basic info on handling vcf file and calculating ibs https://bioconductor.statistik.tu-dortmund.de/packages/3.2/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html

setwd("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/")
hue1_25 = c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744")
vcf = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolsGVCF.filtered_snps_final.PASS.bi.vcf"
genot_dir = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/"
output_direc = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"



#######################################################################
# Calculate kinship
#######################################################################
ibs = make_kinship(vcf = vcf)
kinship = ibs$ibs