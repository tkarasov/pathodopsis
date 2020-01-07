library(dada2)
library(intrval)
#tutorial largely taken from here
#https://benjjneb.github.io/dada2/ITS_workflow.html
#Note that the ITS protocol has some specifics because must be sure to remove primers.
path=("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/")
setwd(path)
filtered_path=("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/filtered_reads")

###REMOVED PA0317 only becayse reverse read was not built

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="ITS_R1", recursive= FALSE))
fnRs <- sort(list.files(path, pattern="ITS_R2", recursive = FALSE))

FWD <- "ACCTGCGGARGGATCA"  ## CHANGE ME to your forward primer sequence
REV <- "GAGATCCRTTGYTRAAAGTT"  ## CHANGE ME...


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients







# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Now Filter and trim
filtFs <- file.path(filtered_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filt_names=sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#missing
hm=gsub(".gz", "", gsub( "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/filtered_reads//filtered/", "", filtRs))

#this next step seemed to have a memory issue so i turned multithreading off
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    maxN=0, maxEE=c(2,2), truncQ=2, minLen=100, rm.phix=TRUE,
                    compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
write.table(out, paste(path,"after_filtering_ITS.csv", sep=""), sep=",", quote=F)
#we are losing about 20-50% of our reads

#bad files -
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

saveRDS(seqtab.nochim, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/ITS_seqtab.rds")
library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/ITS_seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
#seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(st1, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/silva_nr_v132_train_set.fa.gz", multithread=4)
# Write to disk
# saveRDS(tax, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/ITS_tax_final.rds") # CHANGE ME ...
