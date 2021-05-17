#!/usr/bin/env Rscript
library(phyloseq)
library(vegan)

####################
# Load phyloseq object & pathogen identities
####################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds')
my_bac <- read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_bacteria.csv",
  header = TRUE, row.names = 1)
my_oom <- read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_oomycete.csv",
  header = TRUE, row.names = 1)
my_fungi = read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_fungi.csv",
  header = TRUE, row.names = 1)

my_phylo <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/seqtab_final.rds")
my_phylo_tax <- readRDS('/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/tax_final.rds')

####################
# Build Lat/Long Dataframe
####################
clim <- plant_clim$clim_data
rownames(clim) <- clim$Sequence_ID
plant_phylo <- phyloseq(sample_data(clim), refseq = plant_clim$refseq, tax_table = plant_clim$tax_table, phy_tree = plant_clim$phy_tree, otu_table = plant_clim$otu_table )
sample_data(GP1000)$Lat <- as.numeric(as.character(sample_data(GP1000)$Lat))
sample_data(GP1000)$Long <- as.numeric(as.character(sample_data(GP1000)$Long))
#GP1000_dim <- subset_samples(GP1000, sample_data(GP1000)$samples.out %in% colnames(my_bac))
GP1000_dim <- subset_samples(plant_phylo, sample_data(plant_phylo)$Sequence_ID %in% colnames(my_bac))
reorder <- as.character(sample_data(GP1000_dim)$Sequence_ID)
my_bac <- my_bac[, reorder]
my_oom <- my_oom[, reorder]
my_fungi <- my_fungi[, reorder]
all_microb <- rbind(my_bac, my_oom, my_fungi)

load_tot <- data.frame(plant_ID = colnames(my_bac), 
                       oom_load = colSums(my_oom), 
                       fung_load = colSums(my_fungi), 
                       bac_load = colSums(my_bac, na.rm = TRUE),
                       hpa_load = t(my_oom["Peronosporaceae",]),
                       albug_load = t(my_oom["Albuginaceae",]))

sample_data(GP1000_dim)$oom_load = load_tot$oom_load
sample_data(GP1000_dim)$fung_load = load_tot$fung_load
sample_data(GP1000_dim)$bac_load = load_tot$bac_load
sample_data(GP1000_dim)$hpa_load = load_tot$Peronosporaceae
sample_data(GP1000_dim)$albug_load = load_tot$Albuginaceae
sample_data(GP1000_dim)$albug_load = load_tot$Albuginaceae
sample_data(GP1000_dim)$otu5 <- otu_table(GP1000_dim)[,"seq_10"]/1000

####################
# permanova on major variables
####################
# recalculate date:
new_date <- c(sample_data(GP1000_dim)$Date) - 20180200
i = 0
for(rec in new_date){
  i=i+1
  hm = rec
  if(hm<100) rec1= rec
  if(hm>=100 & hm < 200) rec1= 28 +  rec - 100
  if(hm>=200 & hm < 300) rec1= 28 + 31 +  rec - 200
  if(hm>=300 & hm < 400) rec1= 28 + 31 + 30 + rec - 300
  new_date[i]=rec1
}
sample_data(GP1000_dim)$Date <- new_date

c("PDSI", "srad", "vpd", "Elevation", "ws", "pet", "Lat", "Date" )
metadata <- as(sample_data(GP1000_dim), "data.frame")

is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}


metadata[is.nan.data.frame(metadata)] <- 0
GP_poo <- GP1000_dim
otu_table(GP1_poo) <- otu_table(GP1000_dim)/1000[,c(1:50)]
bray <- phyloseq::distance(GP1000_dim, method="bray")

perm <- adonis2(bray ~ Lat +  Date, data = metadata, na.action = na.omit)

Date + PDSI + srad + vpd + Elevation + ws + pet + 



