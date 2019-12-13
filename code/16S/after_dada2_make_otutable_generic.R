library(phyloseq)
library(dada2)
library(dplyr)
library(tidyverse)
library(genefilter)
library("RColorBrewer")

hue1_25 = c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744")
#hue1_25 = c("#bf82cc","#5bc14d","#bf4db3","#9fba36","#7861ce","#4d8c28","#d83e76","#44c181","#d0452f","#4aadd6","#d6812c","#667fc7","#cbaa3b","#9c4769","#7dba6f","#dd809f","#3e8148","#c25d4d","#59c5b0","#de986d","#2f8a72","#91692e","#afb16c","#5f6c2b","#84892d")

#######################################################################
# Reading in data and setting up relevant formats
#######################################################################

path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/amp_seq_functions.R", sep=""))
#output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"
output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/"

seqtab.nochim = readRDS(paste(output_direc,"/seqtab_final.rds", sep="/"))
taxa=readRDS(paste(output_direc,"/tax_final.rds", sep="/"))
metadata=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course2.txt", header=T, sep="\t")


#Data is read in, now we want to subset it properly to rename the controls (1:26) and give metadata
controls <- grep("control", rownames(seqtab.nochim),ignore.case=TRUE,value=TRUE)
control_index <- which(rownames(seqtab.nochim) %in% controls)
samples.out <- sapply(strsplit(rownames(seqtab.nochim),"_"), "[",1)

# keep the control names in the name 
samples.out[control_index]= controls
rownames(seqtab.nochim) = samples.out
metadata_keep=metadata[metadata$Plant_ID%in%samples.out,]
meta_unique = metadata_keep %>% distinct()
metadata_organized=merge(data.frame(samples.out), meta_unique, by.x="samples.out", by.y="Plant_ID", all.x=TRUE)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
samdf <- data.frame(Subject=metadata_organized$samples.out, 
                    Latitude=metadata_organized$Latitude, 
                    Longitude=metadata_organized$Longitude, 
                    Altitude=metadata_organized$Altitude, 
                    hpa=metadata_organized$HpA_plant, 
                    TourID=metadata_organized$Tour_ID, 
                    Clim=metadata_organized$ClimateZ)
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
GP = prune_samples(sample_sums(ps)>=1000, ps)
GP = prune_samples(sample_data(GP)$TourID!="NA", GP)
rm(ps)

#Now for taxa names
dna <- Biostrings::DNAStringSet(taxa_names(GP))
names(dna) <- taxa_names(GP)
GP <- merge_phyloseq(GP, dna)
taxa_names(GP) <- paste0("seq_", seq(ntaxa(GP)))
rm(dna)

#FYI: you can access the sequence by refeseq(GP)

#Remove mitochondria and samples that have fewer than 50 reads in any sample
mito = colnames(otu_table(GP))[which(tax_table(GP)[,5] != "Mitochondria")]
GP = prune_taxa(mito, GP)
flist    <- filterfun(kOverA(1, 50))
GP50 = filter_taxa(GP, flist, TRUE )
GP_fam=tax_glom(GP50, "Family")
#qGPr  = transform_sample_counts(GP, function(otu) otu/sum(otu))
rm(GP)


#######################################################################
# Write OTU table
#######################################################################
#write fasta file of taxa names
taxa_seq = taxa_names(GP)

cat("", file=paste(output_direc, "/16S/16S_all/all_runs/demult_python/all_ASVs.fasta", sep=""), sep='', append=FALSE)

for(i in 1:length(taxa_seq)){
  cat(paste(">seq_",i, sep=""), file="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/all_ASVs.fasta", sep='\n', append=TRUE)
  cat(taxa_seq[i],file="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/all_ASVs.fasta", sep='\n', append=TRUE)
}
  

#######################################################################
# Basic plot of abundance
#######################################################################
#basic plot
top20 <- names(sort(taxa_sums(GP50), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(GP50, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

plot_bar(ps.top20, x="Subject", fill="Family") + facet_wrap(~Clim, scales="free_x") + theme(axis.line.x = element_blank())

GP50.prop <- transform_sample_counts(GP50, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(GP50.prop, method="NMDS", distance="bray")

plot_ordination(GP50.prop, ord.nmds.bray, color="Clim", title="Bray NMDS") +
  theme_light() +
  scale_fill_manual(values = hue1_25)

nmds = plot_ordination(GP50.prop, ord.nmds.bray, justDF = T)

summary(aov(NMDS1~Clim,data=nmds))

#The richness plotting shouldn't be done with filtered data
plot_richness(GP50.prop, x="Subject", measures=c("Shannon", "Simpson"), color="Clim") +  
  scale_fill_manual(values = hue1_25)


# Extract tables
extracted_GP50 = as(otu_table(GP50), "matrix")


#Merge ASV table on genus and family 
#Rename ASV for family and rename duplicates (but why are there duplicates?)
dups = which(duplicated(as.character(tax_table(GP_fam)[,5])))
the_name = paste(as.character(tax_table(GP_fam)[,5])[dups], as.character(dups), sep="_")
tax_table(GP_fam)[,5][dups] = the_name
taxa_names(GP_fam) = as.character(tax_table(GP_fam)[,5])

extracted_GP50_family = as(otu_table(GP_fam), "matrix")
write.table(extracted_GP50_family, paste(output_direc, "/16S/16S_all/all_runs/demult_python/sixteenS_family5.csv", sep=""), sep=",", col.names = T, row.names=T, quote = F)

