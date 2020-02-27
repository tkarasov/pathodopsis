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
clim_gen = read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/terraclim_2_27_2020.csv", sep = " ", header = T)

#my metadata
metadata_old = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course3.txt", header=T, sep="\t") %>% 
  select(c("Sequence_ID", "Tour_ID", "Site_ID", "Plant_ID", "Albugo", "Necrosis", "Used_for", "Host_Species", "ClimateZ", "Land_cov"))

metadata_rebecca = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Meta_table_from_rebecca_2_2020 - Data.tsv", sep ="\t", header = T)

# #merge all metadata
all_metadata = metadata_old %>% full_join(clim_gen, by = c("Site_ID", "Tour_ID")) %>% full_join(metadata_rebecca)
write.table(all_metadata, "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_2_2020_reads.tsv", col.names = TRUE, sep = '\t', quote = FALSE, row.names = F)
#Make variable that is only first letter of Koppen-Geiger
#all_metadata$First_Kop = substr(all_metadata$ClimateZ, 1,1)

metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_2_2020_reads.tsv", header=T, sep="\t", fill =TRUE)

#######################################################################
# Reading in data and setting up relevant formats
#######################################################################

path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/amp_seq_functions.R", sep=""))
output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"
#output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/"

seqtab.nochim = readRDS(paste(output_direc,"/seqtab_final.rds", sep="/"))
taxa = readRDS(paste(output_direc,"/tax_final.rds", sep="/"))
#metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course3.txt", header=T, sep="\t")

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU_tree.RData")

#Data is read in, now we want to subset it properly to rename the controls (1:26) and give metadata
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
metadata_keep=metadata[metadata$Plant_ID%in%samples.out,]
meta_unique = metadata_keep %>% distinct()
metadata_organized=merge(data.frame(samples.out), meta_unique, by.x="samples.out", by.y="Sequence_ID", all.x=TRUE)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

#remove the three duplicated samples

metadata_organized<-metadata_organized[which(!duplicated(metadata_organized$samples.out)),]


# samdf <- data.frame(Subject=metadata_organized$samples.out, 
#                     Latitude=metadata_organized$Lat, 
#                     Longitude=metadata_organized$Long, 
#                     Altitude=metadata_organized$Altitude, 
#                     hpa=metadata_organized$HpA_plant, 
#                     TourID=metadata_organized$Tour_ID, 
#                     Clim=metadata_organized$ClimateZ,
#                     PlantID = metadata_organized$Plant_ID,
#                     Sample_type = metadata_organized$Used_for,
#                     Species = metadata_organized$Host_Species)

samdf <- data.frame(Subject=metadata_organized$samples.out, 
                    Latitude=metadata_organized$Lat, 
                    Longitude=metadata_organized$Long, 
                    Elevation=metadata_organized$Elevation, 
                    #hpa=metadata_organized$HpA_plant, 
                    TourID=metadata_organized$Tour_ID, 
                    Clim=metadata_organized$ClimateZ,
                    PlantID = metadata_organized$Plant_ID,
                    SiteID = metadata_organized$Site_ID,
                    Sample_type = metadata_organized$Used_for,
                    Species = metadata_organized$Host_Species,
                    Soil_temp = metadata_organized$Soil_temp,
                    ClimateZ = metadata_organized$ClimateZ,
                    Land_cov = metadata_organized$Land_cov,
                    tmax = metadata_organized$tmax.6,
                    tmin = metadata_organized$tmax.6,
                    rosette_diam = metadata_organized$R_diameter,
                    herbivory = metadata_organized$Herbivory)

rownames(samdf) <- samdf$Subject
sample_names(seqtab.nochim)=samples.out


#Now let's may the otu_table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

rm(seqtab.nochim,taxa)


#######################################################################
# Filtering data
#######################################################################
#Remove samples with fewwer than 1000 reads
GP0 = prune_samples(sample_sums(ps)>=1000, ps)
GP = prune_samples(is.na(sample_data(GP0)$TourID)==FALSE, GP0)
rm(ps)

# #rename species for soil to "Soil". This doesn't work. Need to add factor level but don't see how to do this in phyloseq object.
# is_soil = which(sample_data(GP)[,"Sample_type"] == "SOIL")
# sample_data(GP)[is_soil,"Species"] = as.character(sample_data(GP)[is_soil,"Species"])
# sample_data(GP)[is_soil,"Species"] <- as.factor("SOIL")


#Now for taxa names
dna <- Biostrings::DNAStringSet(taxa_names(GP))
names(dna) <- taxa_names(GP)
GP <- merge_phyloseq(GP, dna)
taxa_names(GP) <- paste0("seq_", seq(ntaxa(GP)))
rm(dna)

#FYI: you can access the sequence by refeseq(GP)

#Remove mitochondria and ASVs that have fewer than 50 reads in any sample
mito = colnames(otu_table(GP))[which(tax_table(GP)[,5] != "Mitochondria")]
GP = prune_taxa(mito, GP)

flist    <- filterfun(kOverA(1, 50))
GP50 = filter_taxa(GP, flist, TRUE )
GP50 = merge_phyloseq(GP50, fitGTR$tree)
GP_fam=tax_glom(GP50, "Family")

set.seed(4)
GP1000 = rarefy_even_depth(GP50, sample.size = 1000)
#qGPr  = transform_sample_counts(GP, function(otu) otu/sum(otu))


save(GP1000, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds")
save(GP50, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")
save(GP, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP.rds")
rm(GP)
#######################################################################
# Write ASVs to file and OTU table
#######################################################################
#write fasta file of taxa names
taxa_seq = taxa_names(GP)

cat("", file=paste(output_direc, "/all_ASVs.fasta", sep=""), sep='', append=FALSE)

for(i in 1:length(taxa_seq)){
  cat(paste(">seq_",i, sep=""), file=paste(output_direc, "/all_ASVs.fasta", sep=""), sep='\n', append=TRUE)
  cat(taxa_seq[i],file=paste(output_direc, "/all_ASVs.fasta", sep=""), sep='\n', append=TRUE)
}

write_phyloseq_tlk(otu_table(GP), path = output_direc, type = "OTU")
write_phyloseq_tlk(tax_table(GP), path = output_direc, type = "TAXONOMY")
write_phyloseq_tlk(sample_data(GP), path = output_direc, type = "METADATA")
  
# Extract tables
extracted_GP50 = as(otu_table(GP50), "matrix")


#Merge ASV table on genus and family 
#Rename ASV for family and rename duplicates (but why are there duplicates?)
dups = which(duplicated(as.character(tax_table(GP_fam)[,5])))
the_name = paste(as.character(tax_table(GP_fam)[,5])[dups], as.character(dups), sep="_")
tax_table(GP_fam)[,5][dups] = the_name
taxa_names(GP_fam) = as.character(tax_table(GP_fam)[,5])

extracted_GP50_family = as(otu_table(GP_fam), "matrix")
write.table(extracted_GP50_family, paste(output_direc, "/sixteenS_family5.csv", sep=""), sep=",", col.names = T, row.names=T, quote = F)




# 
# #######################################################################
# # Basic plot of abundance
# #######################################################################
# #basic plot
# top20 <- names(sort(taxa_sums(GP50), decreasing=TRUE))[1:20]
# ps.top20 <- transform_sample_counts(GP50, function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20, ps.top20)
# 
# plot_bar(ps.top20, x="Subject", fill="Family") + facet_wrap(~Clim, scales="free_x") + theme(axis.line.x = element_blank())
# 
# GP50.prop <- transform_sample_counts(GP50, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(GP50.prop, method="NMDS", distance="bray")
# 
# plot_ordination(GP50.prop, ord.nmds.bray, color="Clim", title="Bray NMDS") +
#   theme_light() +
#   scale_fill_manual(values = hue1_25)
# 
# plot_ordination(GP50.prop, ord.nmds.bray, color="Species", title="Bray NMDS") +
#   theme_light() +
#   scale_fill_manual(values = hue1_25)
# 
# nmds = plot_ordination(GP50.prop, ord.nmds.bray, justDF = T)
# 
# summary(aov(NMDS1~Clim,data=nmds))
# 
# #The richness plotting shouldn't be done with filtered data
# plot_richness(GP, x="Subject", measures=c("Shannon", "Simpson"), color="Sample_type") +  
#   scale_fill_manual(values = hue1_25)
# 
# 


