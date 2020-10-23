library(phyloseq)
library(dada2)
library(dplyr)
library(tidyverse)
library(genefilter)
library("RColorBrewer")
library(microbiome)

hue1_25 = c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744")
#hue1_25 = c("#bf82cc","#5bc14d","#bf4db3","#9fba36","#7861ce","#4d8c28","#d83e76","#44c181","#d0452f","#4aadd6","#d6812c","#667fc7","#cbaa3b","#9c4769","#7dba6f","#dd809f","#3e8148","#c25d4d","#59c5b0","#de986d","#2f8a72","#91692e","#afb16c","#5f6c2b","#84892d")

#######################################################################
# Merge metadata
#######################################################################

metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_7_2020_reads.tsv",
                      header=T, sep="\t", fill =TRUE)

#######################################################################
# Reading in data and setting up relevant formats
#######################################################################

path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/amp_seq_functions.R", 
             sep=""))
output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"

seqtab.nochim = readRDS(paste(output_direc,"/seqtab_final.rds", sep="/"))

taxa = readRDS(paste(output_direc,"/tax_final.rds", sep="/"))

taxa = addSpecies(taxa, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/silva_species_assignment_v132.fa.gz")

#######################################################################
# Data is read in, now we want to subset it properly to rename the controls 
#######################################################################

control <-"Control|CONTROL|control"
controls <- grep(control, rownames(seqtab.nochim),ignore.case=TRUE,value=TRUE)
control_index <- which(rownames(seqtab.nochim) %in% controls)
empty <- grep("empty", rownames(seqtab.nochim),ignore.case=TRUE,value=TRUE)
empty_index <- which(rownames(seqtab.nochim) %in% empty)
samples.out <- sapply(strsplit(rownames(seqtab.nochim),"_plate"), "[",1)

# keep the control names in the name 
samples.out[control_index]= paste(controls, c(1:length(controls)), sep = "_")
samples.out[empty_index] = paste(empty, c(1:length(empty)), sep = "_")
rownames(seqtab.nochim) = samples.out
metadata_keep = metadata[metadata$Plant_ID%in%samples.out,]
meta_unique = metadata_keep %>% distinct()
metadata_organized =
	merge(data.frame(samples.out), 
		meta_unique, 
		by.x="samples.out", 
		by.y="Sequence_ID", all.x=TRUE)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

#remove the three duplicated samples

metadata_organized<-metadata_organized[which(!duplicated(metadata_organized$samples.out)),]

samdf <- data.frame(metadata_organized)

rownames(samdf) <- samdf$samples.out
  #samdf$Subject
sample_names(seqtab.nochim)=samples.out


#Now let's make the otu_table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

rm(seqtab.nochim,taxa)


#######################################################################
# Initial Filtering of data to remove samples with fewer than 1000 reads
#######################################################################
#Remove samples with fewer than 1000 reads
GP0 = prune_samples(sample_sums(ps)>=1000, ps)
GP = prune_samples(is.na(sample_data(GP0)$Tour_ID)==FALSE, GP0)
rm(ps)

#######################################################################
# Now for taxa names
#######################################################################
dna <- Biostrings::DNAStringSet(taxa_names(GP))
names(dna) <- taxa_names(GP)
GP <- merge_phyloseq(GP, dna)
taxa_names(GP) <- paste0("seq_", seq(ntaxa(GP)))
rm(dna)

#FYI: you can access the sequence by refeseq(GP)

#######################################################################
# Now filter mitochondria and taxa that are never observed with 50 reads in a sample
#######################################################################
mito = colnames(otu_table(GP))[which(tax_table(GP)[,5] != "Mitochondria")]
GP = prune_taxa(mito, GP)
flist    <- filterfun(kOverA(1, 50))
GP50 = filter_taxa(GP, flist, TRUE )

#######################################################################
# Build the tree
#######################################################################
# the first time this was run I did an alignment but now can just load which outputs the OTU_tree.RData file
# source("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/amplicon_general/after_dada2_make_tree.R)

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU_tree.RData")

GP50 = merge_phyloseq(GP50, fitGTR$tree)
GP_fam=tax_glom(GP50, "Family")

set.seed(4)
GP1000 = rarefy_even_depth(GP50, sample.size = 1000)

save(GP1000, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds")
save(GP50, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")
save(GP, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP.rds")
rm(GP)

# #######################################################################
# # Get dominant phylotypes
# #######################################################################
source("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/dominant_phylotypes.R")

# #######################################################################
# # Write ASVs to file and OTU table
# #######################################################################

# Extract tables
extracted_GP50 = as(otu_table(GP50), "matrix")

#######################################################################
# Merge ASV table on genus and family 
#######################################################################
#Merge ASV table on genus and family 
#Rename ASV for family and rename duplicates (but why are there duplicates?)
dups = which(duplicated(as.character(tax_table(GP_fam)[,5])))
the_name = paste(as.character(tax_table(GP_fam)[,5])[dups], as.character(dups), sep="_")
tax_table(GP_fam)[,5][dups] = the_name
taxa_names(GP_fam) = as.character(tax_table(GP_fam)[,5])

extracted_GP50_family = as(otu_table(GP_fam), "matrix")
write.table(extracted_GP50_family, paste(output_direc, "/sixteenS_family5.csv", sep=""), sep=",", col.names = T, row.names=T, quote = F)

