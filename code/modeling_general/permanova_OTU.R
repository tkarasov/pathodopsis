library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)

# The goal of this function is to run a permanova on the OTUs in a dataframe with the effect of climate
vcf = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolsGVCF.filtered_snps_final.PASS.bi.vcf"
genot_dir = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/"
output_direc = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/otu_table.csv"



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
my.rarefy <- rarefy_even_depth(my_phyloseq)
bc.rare <- vegdist(as(otu_table(my.rarefy),"matrix"), method = "bray")


temp = plot_ordination_tlk(my.rarefy, ordinate(my.rarefy, method = "PCoA", distance = "bray"), "MDS", color = "Species") + 
  geom_point(size = 5)

+ scale_color_manual(values = hue1_25)


#######################################################################
# permanova with climate
#######################################################################
metadata <- as(sample_data(my_phyloseq), "data.frame")
adonis()

#bray curtis distance matrix with permanova
rare.bray <- phyloseq::distance(my.rarefy, method = "bray")
adonis(rare.bray ~  Sample_type + Clim,
       data = metadata)
