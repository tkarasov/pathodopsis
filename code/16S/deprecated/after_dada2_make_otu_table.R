#library(phyloseq)
library(dada2)
library(dplyr)
library(tidyverse)
library(fossil)

#library(msa)
#library(DECIPHER)

library(genefilter)
library(phangorn)
library("RColorBrewer")
library(gplots)
library(sjstats)
library(nlme)

#path="/Users/tkarasov/work_main"

#choose whether working on home computer or at work
path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/scripts/16S/amp_seq_functions.R", sep=""))

#https://f1000research.com/articles/5-1492/v1 This is the dada2 file

output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all"

seqtab.nochim = readRDS(paste(output_direc,"seqtab_final.rds", sep="/"))

taxa=readRDS(paste(output_direc,"tax_final.rds", sep="/"))

#metadata=read.table(paste(path,"/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/v1_22_5_merged.txt", sep=""), header=T, sep=",")
#koppen_geiger=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course2.txt", header=T, sep="\t")
metadata=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course2.txt", header=T, sep="\t")
#metadata=merge(metadata, koppen_geiger, by=c("Site.ID", "Latitude", "Longitude"))

samples.out <- rownames(seqtab.nochim)

metadata_keep=metadata[metadata$Plant_ID%in%samples.out,]

meta_unique = metadata_keep %>% distinct()

metadata_organized=merge(data.frame(samples.out), meta_unique, by.x="samples.out", by.y="Plant_ID", all.x=TRUE)

subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

samdf <- data.frame(Subject=metadata_organized$samples.out, Latitude=metadata_organized$Latitude, Longitude=metadata_organized$Longitude, Altitude=metadata_organized$Altitude, hpa=metadata_organized$HpA_plant, TourID=metadata_organized$Tour_ID, Clim=metadata_organized$ClimateZ)

rownames(samdf) <- samdf$Subject

sample_names(seqtab.nochim)=samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#remove samples with fewwer than 1000 reads
GP = prune_samples(sample_sums(ps)>=1000, ps)

#now for taxa names
dna <- Biostrings::DNAStringSet(taxa_names(GP))
names(dna) <- taxa_names(GP)
GP <- merge_phyloseq(GP, dna)
taxa_names(GP) <- paste0("ASV", seq(ntaxa(GP)))

#you can access the sequence by refeseq(GP)

#remove mitochondria
mito = colnames(otu_table(GP))[which(tax_table(GP)[,5] != "Mitochondria")]
GP = prune_taxa(mito, GP)
# remove samples that have fewer than 50 reads in any sample
flist    <- filterfun(kOverA(1, 50))
GP50 = filter_taxa(GP, flist, TRUE )

qGPr  = transform_sample_counts(GP, function(otu) otu/sum(otu))


#basic plot
top20 <- names(sort(taxa_sums(GP50), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(GP50, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Subject", fill="Family") + facet_wrap(~Clim, scales="free_x")

GP50.prop <- transform_sample_counts(GP50, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(GP50.prop, method="NMDS", distance="bray")
plot_ordination(GP50.prop, ord.nmds.bray, color="TourID", title="Bray NMDS")

plot_richness(GP50.prop, x="Subject", measures=c("Shannon", "Simpson"), color="Clim")


# Extract tables
extracted_GP50 = as(otu_table(GP50), "matrix")

#merge on genus and family 
GP_genus=tax_glom(GP50, "Genus")
GP_fam=tax_glom(GP50, "Family")

##rename ASV for family


extracted_GP50_family = as(otu_table(GP_fam), "matrix")


#make phylogenetic tree. DECIPHER Alignseqs scales linearly rather than exponentially. Cannot be run on my computer. 
#seqs <- getSequences(seqtab.nochim)
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#mult <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=TRUE, processors=16)
#align seqs didn't maintain names

#now write decipher alignment to file 
#https://github.com/benjjneb/dada2/issues/204
#writeXStringSet(mult, file= paste(output_direc,"/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/16S_otus_aligned.fasta", sep=""))

#fasttree
#system("export OMP_NUM_THREADS=16")
#system("/usr/bin/fasttreeMP -fastest -noml -gtr -nt /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/16S_otus_aligned.fasta > /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/16S_otus_aligned.FastTree.tre")

# Alternative for trees...as of right now it seems like we need to use alignseqs to align the OTUs but then put throught raxml
#alignment.rax.gtr <- raxml(DNAbin(mult), m="GTRGAMMAIX", # model f="a", # best tree and bootstrap p=1234, # random number seed x=2345, # random seed for rapid bootstrapping N=100, # number of bootstrap replicates file="alignment", # name of output files exec="raxmlHPC", # name of executablethreads=16S)

#plants not in the metadata  "PA0824" "PA0825" "PA1336" "PA1339" "PA1749" "PA1753" "PA1756" "PA1757" "PA1759" "PA1761" "PC0029"
#st <- seqtab.nochim[rowSums(seqtab.nochim) >= 500,]
#seqtab.nochim <- st
