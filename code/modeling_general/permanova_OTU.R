library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
setwd("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/")
hue1_25 = c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744")

# The goal of this script is to run a permanova on the OTUs in a dataframe with the effect of climate

vcf = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolsGVCF.filtered_snps_final.PASS.bi.vcf"
genot_dir = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/"
output_direc = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"
ter
# basic info on handling vcf file and calculating ibs https://bioconductor.statistik.tu-dortmund.de/packages/3.2/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html

#######################################################################
# Calculate kinship
#######################################################################
ibs = make_kinship(vcf = vcf)
snpset = make_snpset(vcf = vcf, ld.threshold = 0.2)

kinship = ibs$ibs

#######################################################################
# Read in phyloseq object
#######################################################################
otu.table <-paste(output_direc, "otu_table.csv", sep="")
taxonomy.table <-paste(output_direc, "taxonomy_table.csv", sep="")
metadata.table <- paste(output_direc, "metadata_table.csv", sep="")
my_phyloseq <- read_phyloseq(otu.file = otu.table, metadata.file = metadata.table, taxonomy.file = taxonomy.table, type = "simple")

#######################################################################
# Rarefy data for some diversity assessments
#######################################################################
my.rarefy <- rarefy_even_depth(my_phyloseq, rngseed = 416)
bc.rare <- vegdist(as(otu_table(my.rarefy),"matrix"), method = "bray")
bc.pcoa <- ape::pcoa(bc.rare, correction = "lingoes")

hm =  ordinate(my.rarefy, method = "MDS", distance = "bray")
#The ape package PCoA seems to have a problem with my dataframes: https://github.com/joey711/phyloseq/issues/594 and outputting the eigenvectors from PCoA in non-decreasing order.

species = plot_ordination(my.rarefy, ordinate(my.rarefy, method = "NMDS", distance = "bray"), type = "samples", color = "Species") + 
  geom_point(size = 5) +
scale_color_manual(values = hue1_25)

clim = plot_ordination_tlk(my.rarefy, ordinate(my.rarefy, method = "DCA", distance = "bray"), "MDS", color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25)

samp = plot_ordination_tlk(my.rarefy, ordinate(my.rarefy, method = "DCA", distance = "bray"), "MDS", color = "Sample_type") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/ordination_microbiome.pdf",family = "ArialMT", useDingbats = F)
plot_grid(species, clim, samp, ncol = 1, align = "hv")
dev.off()


save.image(file = "my_work_space.RData")
load("my_work_space.RData")
#######################################################################
# permanova with climate
#######################################################################
metadata <- as(sample_data(my_phyloseq), "data.frame")
 

#bray curtis distance matrix with permanova
rare.bray <- phyloseq::distance(my.rarefy, method = "bray")
rare.jsd <- phyloseq::distance(my.rarefy, method = "jsd")

bray.perm <- adonis(rare.bray ~  Sample_type + Clim,
       data = metadata)

jsd.perm <- adonis(rare.jsd ~  Sample_type + Clim,
                    data = metadata)

#######################################################################
# Subset to only metagenomes and perform different ordinations
#######################################################################

meta_only = subset_samples(my.rarefy, Sample_type == "M") 
keep_samples = c(sample_data(meta_only)[,"Subject"])$Subject
metadata_only = metadata[which(metadata$Subject %in% keep_samples),]

rare.bray.M <- phyloseq::distance(meta_only, method = "bray")

hm.NMDS = ordinate(meta_only, method = "NMDS", distance = "bray", trymax = 50)

hm.MDS = ordinate(meta_only, method = "MDS", distance = "bray")

hm.DCA = ordinate(meta_only, method = "DCA", distance = "bray")

clim.meta.NMDS = plot_ordination(meta_only, hm.NMDS, color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

clim.meta.MDS = plot_ordination(meta_only, hm.MDS, color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

clim.meta.DCA = plot_ordination(meta_only, hm.DCA, color = "Clim") + 
  geom_point(size = 5) +
  scale_color_manual(values = hue1_25) + 
  theme_bw()

bray.perm <- adonis(rare.bray.M ~  Clim + Species + TourID,
                    data = metadata_only)


#######################################################################
# Rank threshold ordination
#######################################################################



