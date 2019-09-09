library(dada2)
library(intrval)
#tutorial largely taken from here

path="/Users/tkarasov/work_main"
#path="/ebio"
proc_path=(paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/", sep=""))
setwd(proc_path)
filtered_path=(paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/filtered_reads/", sep=""))

###REMOVED PA0317 only becayse reverse read was not built

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(proc_path, pattern="16S_R1", recursive= FALSE))
fnRs <- sort(list.files(proc_path, pattern="16S_R2", recursive = FALSE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Now Filter and trim
filtFs <- file.path(filtered_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#this next step seemed to have a memory issue so i turned multithreading off
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,210),
                    maxN=0,  truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=3) # On Windows set multithread=FALSE

write.table(out, paste(path,"after_filtering_16S.csv", sep=""), sep=",", quote=F)
#maxEE=c(5,5) this was removing so many of the reads. It will slow down the filtering but...
head(out)
#we are losing about 20-50% of our reads

#bad files 
#PA0754
#PA1923
#PC0153


#Learn the Error Rates
#The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).
#The following runs in about 3 minutes on a 2013 Macbook Pro:
  
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


#233 which(filtFs=="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/filtered_reads//filtered/PA0754_F_filt.fastq.gz")
#573 which(filtFs=="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/filtered_reads//filtered/PA1923_F_filt.fastq.gz")
#598 which(filtFs=="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/filtered_reads//filtered/PC0153_F_filt.fastq.gz")

filtFs_new=filtFs[-c(573,233,598)][c(1:631)]
filtRs_new=filtRs[-c(573,233,598)][c(1:631)]
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(filtFs_new))
names(mergers)= sapply(strsplit(basename(filtFs_new), "_"), `[`, 1)
#names(mergers) <- sample.names
for(sam in names(mergers)) {
  cat("Processing:", sam, "\n")
  full_name=paste(paste(filtered_path,"filtered/", sep=""), sam, sep="")
  derepF <- derepFastq(paste(full_name,"_F_filt.fastq.gz", sep=""))
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(paste(full_name,"_R_filt.fastq.gz", sep=""))
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
 rm(derepF); rm(derepR)


#now make the sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#now remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



#build phylogeny of sequences
#make phylogenetic tree
#seqs <- getSequences(seqtab.nochim)
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#mult <- msa(seqs, method="ClustalW", type="dna", order="input")

#library("phangorn")
#phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
#dm <- dist.ml(phang.align)
#treeNJ <- NJ(dm) # Note, tip order != sequence order
#fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
#detach("package:phangorn", unload=TRUE)


saveRDS(seqtab.nochim, paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/16S_seqtab.rds", sep=""))
library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
#st1 <- readRDS("path/to/run1/output/seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
#seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/silva_nr_v132_train_set.fa.gz", sep=""), multithread=8, verbose=TRUE)
spec <- assignSpecies(seqtab.nochim, paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/silva_species_assignment_v132.fa.gz",sep=""), verbose=TRUE)
# Write to disk
saveRDS(tax, paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/16S_tax_final.rds", sep="")) # CHANGE ME ...
# saveRDS(spec, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/16S_spec_final.rds")
#write out fasta
uniquesToFasta(seqtab.nochim,paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/16S_tax_final_uniques.fasta", sep=""))

