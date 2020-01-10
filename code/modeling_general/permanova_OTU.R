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

# The goal of this script is to run a permanova on the OTUs in a dataframe with the effect of climate
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

#######################################################################
# Read in Phyloseq object
#######################################################################
# Load the phyloseq object generated in after_dada2_make_otutable
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")

# #######################################################################
# # Read in phyloseq object
# #######################################################################
# otu.table <-paste(output_direc, "otu_table.csv", sep="")
# taxonomy.table <-paste(output_direc, "taxonomy_table.csv", sep="")
# metadata.table <- paste(output_direc, "metadata_table.csv", sep="")
# my_phyloseq <- read_phyloseq(otu.file = otu.table, metadata.file = metadata.table, taxonomy.file = taxonomy.table, type = "simple")
# 
# #######################################################################
# # Prune phyloseq to only those samples with at least 1000 reads, and to those ASVs that 
# #######################################################################
# GP = prune_samples(sample_sums(my_phyloseq)>=1000, my_phyloseq)
# GP = prune_samples(sample_data(GP)$TourID!="NA", GP)
# flist    <- filterfun(kOverA(1, 50))
# GP50 = filter_taxa(GP, flist, TRUE )
# 

#######################################################################
# Rarefy data then do PCoA
#######################################################################
my.rarefy <- rarefy_even_depth(GP50, rngseed = 416)
my.RA <- transform(my.rarefy, "compositional")
#my.sqrt <- transform(my.rarefy, )

#bc.rare <- vegdist(as(otu_table(my.rarefy),"matrix"), method = "bray")
#bc.pcoa <- ape::pcoa(bc.rare, correction = "lingoes")

# hm =  ordinate(my.RA, method = "MDS", distance = "bray")
#The ape package PCoA seems to have a problem with my dataframes: https://github.com/joey711/phyloseq/issues/594 and outputting the eigenvectors from PCoA in non-decreasing order. I haven't been able to get MDS (PCoA) to work as of 1/2020 with the full dataset but when I pruned it worked, and am currently giving up

species = plot_ordination(my.rarefy, ordinate(my.rarefy, method = "MDS", distance = "bray"), type = "samples", color = "Species") + 
  geom_point(size = 5) +
scale_color_manual(values = hue1_25) +
  theme_bw()

clim = plot_ordination(my.rarefy, ordinate(my.rarefy, method = "MDS", distance = "bray"), type = "samples", color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) +
  theme_bw()

samp = plot_ordination(my.rarefy, ordinate(my.rarefy, method = "MDS", distance = "bray"), type = "samples", color = "Sample_type") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/PCoA_ordination_microbiome.pdf",family = "ArialMT", useDingbats = F)
plot_grid(species + theme(axis.text.x = element_blank()), clim, ncol = 1, align = "hv")
dev.off()

save.image(file = "my_work_space.RData")
load("my_work_space.RData")

#######################################################################
# permanova with climate
#######################################################################
metadata <- as(sample_data(GP), "data.frame")
pslog <- transform_sample_counts(my.rarefy, function(x) log(1 + x))
#out.bc.log <- ordinate(pslog, method = "NMDS", distance = "bray")

#bray curtis distance matrix with permanova
rare.bray <- phyloseq::distance(my.rarefy, method = "bray")

rare.jsd <- phyloseq::distance(my.rarefy, method = "jsd")

bray.perm <- adonis(rare.bray ~  Sample_type + Clim,
       data = metadata)

# jsd.perm <- adonis(rare.jsd ~  Sample_type + Clim,
#                    data = metadata)

#######################################################################
# Subset to only metagenomes and perform different ordinations
#######################################################################

meta_only = subset_samples(my.rarefy, Sample_type == "M") 
keep_samples = c(sample_data(meta_only)[,"Subject"])$Subject
metadata_only = metadata[which(metadata$Subject %in% keep_samples),]


rare.bray.M <- phyloseq::distance(meta_only, method = "bray")

# hm.NMDS = ordinate(meta_only, method = "NMDS", distance = "bray", trymax = 50)

hm.MDS = ordinate(meta_only, method = "MDS", distance = "bray")

# hm.DCA = ordinate(meta_only, method = "DCA", distance = "bray")

# clim.meta.NMDS = plot_ordination(meta_only, hm.NMDS, color = "Clim") + 
#   geom_point(size = 5) +
#   scale_color_manual(values = hue1_25) + 
#   theme_bw()

clim.meta.MDS = plot_ordination(meta_only, hm.MDS, color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

species.meta.MDS = plot_ordination(meta_only, hm.MDS, color = "Species") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/plant_PCoA_ordination_microbiome.pdf",family = "ArialMT", useDingbats = F)
plot_grid(species.meta.MDS + theme(axis.text.x = element_blank()), clim.meta.MDS, ncol = 1, align = "hv")
dev.off()

# clim.meta.DCA = plot_ordination(meta_only, hm.DCA, color = "Clim") + 
#   geom_point(size = 5) +
#   scale_color_manual(values = hue1_25) + 
#   theme_bw()

bray.perm.meta <- adonis(rare.bray.M ~  Clim + Species + TourID,
                    data = metadata_only)

#######################################################################
# Subset to only A. thaliana
#######################################################################
athal_rarefy = subset_samples(my.rarefy, Species == "Ath") 
athal_meta = metadata_only[which(metadata_only$Species == "Ath"),]

rare.bray.athal <- phyloseq::distance(athal_rarefy, method = "bray")
athal.MDS = ordinate(athal_rarefy, method = "MDS", distance = "bray")

clim.athal.MDS = plot_ordination(athal_rarefy, athal.MDS, color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

bray.perm.athal <- adonis2(rare.bray.athal ~  as.factor(Clim) + TourID,
                         data = athal_meta)

# #######################################################################
# # Rank threshold ordination
# #######################################################################
# # https://f1000research.com/articles/5-1492/v2#f28
# # Using PCoA on non-euclidean distance is problematic
# abund <- otu_table(pslog)
# abund_ranks <- t(apply(abund, 1, frank))
# 
# abund_ranks[abund_ranks < 1] <- 1
# abund_ranks <- t(abund_ranks) 
# 
# ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 4) #keep 4 axes. Rows are individuals and columns are the different variables 
# row_scores <- data.frame(li = ranks_pca$li,
#                          SampleID = rownames(abund_ranks))
# col_scores <- data.frame(co = ranks_pca$co,
#                          seq = colnames(abund_ranks))
# 
# tax <- tax_table(pslog)@.Data %>%
#   data.frame(stringsAsFactors = FALSE)
# tax$seq <- rownames(tax)


# #######################################################################
# # PCA on deseq transformed data. In reality, it's not OK to use PCoA and I also was not using normalized data. Here I use limma to normalize and do PCA on the normalized data
# #######################################################################
# Susan Holmes suggests just doing PCA on the deseq transformed data https://github.com/joey711/phyloseq/issues/492
#add 1 to every sample to enable geometric mean calculation in the variance stabilization
GP_trans = GP50
otu_table(GP_trans) = otu_table(GP_trans) + 1

m = t(as(otu_table(GP_trans), "matrix"))
taxonomy = data.frame(as(tax_table(GP_trans), "matrix"))
samp_data = (as(sample_data(GP_trans), "matrix"))

#GP50_deseq <- phyloseq_to_deseq2(GP_trans, ~ PlantID + Sample_type)
#trial.DEseq = estimateSizeFactors(GP50_deseq)
#trial.DEseq = estimateDispersions(trial.DEseq,quiet = FALSE)
#trial.vst = getVarianceStabilizedData(trial.DEseq)

# Deseq is too slow. Try with limma/voom.  Tutorial here :https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# Now turn into a DGEList
d = DGEList(counts = m, genes = (taxonomy), samples = samp_data, remove.zeros = TRUE)
#d = DGEList(counts = m)
# Calculate the normalization factors
z = calcNormFactors(d, method="RLE")

# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(z$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}

# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + PlantID + Sample_type, data = data.frame(as(sample_data(GP_trans), "matrix")))

# Voom transforms count data to log2 counts per million, and computes observation-level weights. The data can then be handled for ordiation
v <- voom(z, mm, plot=TRUE)
#fit <- lmFit(v, mm)
#fit <- eBayes(fit, robust=TRUE)

GP_voom = phyloseq(otu_table(t(v$E), taxa_are_rows=FALSE), 
                   sample_data(GP_trans), 
                   tax_table(GP_trans))

save(GP_voom, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP500_voom.rds")

#otu_table(GP_voom) = t(v$E)
pca.all <- prcomp(t(na.omit(v$E)))


