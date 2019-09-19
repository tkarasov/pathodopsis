library(phyloseq)
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

path="/Users/tkarasov/work_main"

#choose whether working on home computer or at work
path="/ebio"
source(paste(path, "/abt6_projects9/pathodopsis_microbiomes/scripts/16S/amp_seq_functions.R", sep=""))

#https://f1000research.com/articles/5-1492/v1 This is the dada2 file

output_direc="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all"
seqtab.nochim = readRDS(paste(output_direc,"seqtab_final.rds", sep="/"))
taxa=readRDS(paste(output_direc,"tax_final.rds", sep="/"))
metadata=read.table(paste(path,"/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/v1_22_5_merged.txt", sep=""), header=T, sep=",")

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
samples.out <- rownames(seqtab.nochim)
meta_unique = metadata_keep %>% distinct()

metadata_organized=merge(data.frame(samples.out), meta_unique, by.x="samples.out", by.y="Plant_ID", all.x=TRUE)

subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
samdf <- data.frame(Subject=metadata_organized$samples.out, Latitude=metadata_organized$Latitude, Longitude=metadata_organized$Longitude, Altitude=metadata_organized$Altitude.x, temp=metadata_organized$Air.temp, humid=metadata_organized$Air.humidity, hpa=metadata_organized$HpA_plant, TourID=metadata_organized$Tour.ID)
rownames(samdf) <- samdf$Subject
sample_names(seqtab.nochim)=samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#remove samples with fewwer than 1000 reads
GP = prune_samples(sample_sums(ps)>=500, ps)
GPr  = transform_sample_counts(GP, function(otu) otu/sum(otu))

#Filter for only those OTUs that make up at least .1
GP_0.01 = filter_taxa(GPr, kOverA(4, 0.01), TRUE)

#Plot the ranks of the most abundant phyla
pdf(paste(path,"/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/raks_genera.pdf", sep=""))
plot_taxa_summary(GP_0.01, "Genus")
dev.off()

summarize_taxa(GP_0.01, "Genus")
hm=as.matrix(data.frame(otu_table(GP_0.01)))
#remove second column because it maps to B.olerea (the capsella sequences)
hm=hm[,-c(2)]
#choose only those samples that are included at present
meta_subset=metadata_organized[which(metadata_organized$samples.out%in%rownames(hm)),]

#we want to sort hm according first to trip then to latitude
meta_sorted <- arrange(meta_subset, Tour.ID)
hm_sorted = hm[match(as.character(meta_sorted$samples.out),rownames(hm)),]
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
x=as.factor(meta_sorted$Tour.ID)
levels(x)=1:length(levels(x))
x=as.numeric(x)

pdf(paste(path,"/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/otu_heatmap.pdf", sep=""))
heatmap.2(log10(hm_sorted+0.00000001),  trace='none', col=brewer.pal(11,"YlOrRd"), Colv = FALSE, Rowv = FALSE, RowSideColors=col_vector[x])
dev.off()
#

#test for effect of latitude and longitude for every OTU
spatial_sig<-function(hm_sorted, meta_sorted){
  otu_vector=numeric(length=dim(hm_sorted)[2])
  for(i in 1:dim(hm_sorted)[2]){
    temp=cbind(hm_sorted[,i], meta_sorted)
    colnames(temp)[1]="otu"
    fit = lm(otu~Latitude, data=temp)
    sig= summary(fit)$coefficients[,4][2]
    #unrestricted_fit=lme(otu~(1|Pop_size)+Latitude, data=temp)
    #restricted_fit=lmer(otu~(1|Pop_size), data=temp)
    #AIC_temp=AIC(restricted_fit)-AIC(unrestricted_fit)
    #print(y)
    otu_vector[i]=sig
  }
  
}
  
  spatial_coeff<-function(hm_sorted, meta_sorted){
    otu_vector=numeric(length=dim(hm_sorted)[2])
    for(i in 1:dim(hm_sorted)[2]){
      temp=cbind(hm_sorted[,i], meta_sorted)
      colnames(temp)[1]="otu"
      fit = lm(otu~Latitude, data=temp)
      sig= summary(fit)$coefficients[,1][2]
      #unrestricted_fit=lme(otu~(1|Pop_size)+Latitude, data=temp)
      #restricted_fit=lmer(otu~(1|Pop_size), data=temp)
      #AIC_temp=AIC(restricted_fit)-AIC(unrestricted_fit)
      #print(y)
      otu_vector[i]=sig
    }
  
  
  #FDR correct
  otu_fin=p.adjust(otu_vector, "BH")
  return(otu_vector)
}

otu_sig= spatial_sig(hm_sorted, meta_sorted)
otu_sig_map=rbind(hm_sorted, otu_sig)
tax_as=tax_table(GP_0.01)[,5][-c(2)]
otu_sig_map=rbind(otu_sig_map,t(tax_as))
rownames(otu_sig_map)[c((dim(otu_sig_map)[1]-1):(dim(otu_sig_map)[1]))]=c("FDR_padjust","tax_class")

#now put taxonomy onto the significance map
FDR=which(rownames(otu_sig_map)=="FDR_padjust")
is_sig=otu_sig_map[,which(otu_sig_map[FDR,]<=0.05)]
#t_is_sig=t(is_sig)[,-c(FDR, FDR+1)]
t_is_sig=data.frame(t(is_sig))
#t_is_sig=apply(t_is_sig, 2, as.numeric)
#t_ist_sig=subset(t_is_sig, select= !c("FDR_padjust","tax_class"))
lat_long=metadata[match(rownames(is_sig), metadata$Plant_ID),][,c("Plant_ID", "Latitude", "Longitude")]
is_sig_lat_long <- cbind(lat_long, is_sig)

temp=is_sig_lat_long[,c(1:4)]
#29/224 OTUs show significant geographic patterning
#Relationship between geographic distance and microbiome distance
#geographic distance
lat_long=metadata_organized[,c( "Longitude", "Latitude")]
rownames(lat_long)=metadata_organized$samples.out
dist_mat=as.matrix(earth.dist(lat_long))
diag(dist_mat)=NA
rownames(dist_mat)=rownames(lat_long)
colnames(dist_mat)=rownames(lat_long)


###Look at principal coordinates of data and perform regression with environmental variables
#Transform data to proportions as appropriate for Bray-Curtis distances
#Hellinger transform
GP_0.01_H <-GP_0.01
otu_table(GP_0.01_H)<-otu_table(sqrt(otu_table(GP_0.01)))
ord.pcoa.bray <- ordinate(GP_0.01_H, method="PCoA", distance="bray")
pdf(paste(path,"/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/bray_16S_pcoa.pdf", sep=""))
plot_ordination(GP_0.01_H, ord.pcoa.bray, color="Latitude", title="Bray PCoA") +theme_bw()
dev.off()
#plot_ordination_utils(GP_0.01, ord.nmds.bray, color="Latitude", title="Bray NMDS") +theme_bw()
axis1=ord.pcoa.bray$vectors[,1]
axis2=ord.pcoa.bray$vectors[,2]
meta_axis=cbind(sample_data(GP_0.01_H), axis1)
meta_axis=cbind(meta_axis, axis2)
summary(lm(axis1~Latitude+Longitude+temp, data=meta_axis))

r2(lmer(axis1~(1|TourID)+Latitude+Longitude, data=meta_axis))
r2(lmer(axis2~(1|TourID)+Latitude+Longitude, data=meta_axis))

m1 <- gls(axis1 ~ Latitude, data=meta_axis, na.action=na.omit)
re_var(m, adjusted = TRUE)
m2 <- gls(axis1 ~ Latitude+Longitude, data=meta_axis, na.action=na.omit)
re_var(m, badjusted = TRUE)
m3 <- lme(axis1 ~ 1, random = ~ 1|TourID, data=meta_axis, na.action=na.omit)
re_var(m, adjusted = TRUE)
my_aov <- anova(update(m1, . ~ ., method = "ML"),update(m2, . ~ ., method = "ML"), update(m3, . ~ ., method = "ML") )











#
bac=subset_taxa(GP, Kingdom=="Bacteria")
#bac_melt=psmelt(bac)



#Relationship between geographic distance and microbiome distance
#geographic distance
lat_long=metadata_organized[,c( "Longitude", "Latitude")]
rownames(lat_long)=metadata_organized$samples.out
dist_mat=as.matrix(earth.dist(lat_long))
diag(dist_mat)=NA
rownames(dist_mat)=rownames(lat_long)
colnames(dist_mat)=rownames(lat_long)

#bray curtis distance
bc=as.matrix(phyloseq::distance(otu_table(GP_0.01), method="bray"))

#Relationship between bray-curtis and geographic
bc_dist=as.data.frame(cbind(c(bc), c(dist_mat)))
colnames(bc_dist)=c("bray_curtis", "distance_km")
bc_dist=bc_dist[is.na(bc_dist[,1])==F &is.na(bc_dist[,2])==F, ]
smoothed = loess(c(bc_dist[,1])~c(bc_dist[,2]), span=.5)
smoothed10=predict(smoothed)
plot(x=c(bc_dist[,2]), y=c(bc_dist[,1]),  main="Loess Smoothing and Prediction", xlab="Distance (km)", ylab="Bray-Curtis Distance", pch=20)
lines(smoothed10, x=c(bc_dist[,2]), col="red")

locpoly(x=c(dist_mat), y=c(bc))


#within 100km
bc_10=bc_dist[bc_dist$distance_km<=10,]


















###############Probably don't use

#msa(all_seq, method="ClustalW", type="dna", order="input")
library("phangorn")
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab.nochim))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)






###################
dist_methods <- unlist(distanceMethodList)
dist_methods[(1:3)]
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# This is the user-defined method:
dist_methods["designdist"]
print(dist_methods)
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(GP_0.01_H, method=i)
  # Calculate ordination
  iMDS  <- ordinate(GPr, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(GPr, iMDS, color="Latitude", shape="hpa")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}


#GPr of Pseudomonas
GP.Ps=subset_taxa(GP_0.01, Genus=="Pseudomonas")
#GP.Ps.filter = filter_taxa(GP.Ps, function(x) sum(x > 3) > (0.001*length(x)), TRUE)
#filterfun_sample(GP)
#plot_tree(GP.Ps, color = "SampleType", shape = "Genus", label.tips = "Genus", size = "abundance", plot.margin = 0.5, ladderize = TRUE)
#GPr.ps  = transform_sample_counts(GP.Ps.filter, function(otu) otu/sum(otu))
#hm=data.frame(otu_table(GPr))
#GPfr = filter_taxa(GP, function(x) mean(x) > 1e-100, TRUE)
#plot_bar(GPr, fill="Genus")
plot_richness(GP.Ps, x="Latitude", measures=c("Shannon", "Simpson"))
# Transform data to proportions as appropriate for Bray-Curtis distances
ord.nmds.bray <- ordinate(GP.Ps, method="NMDS", distance="bray")
ord.nmds.jaccard <- ordinate(GP.Ps, method="NMDS", distance="jaccard")
plot_ordination(GP.Ps, ord.nmds.bray, color="Latitude", title="Bray NMDS") +theme_bw()



plot_heatmap(GP.Ps, "PCA", "bray",taxa.label="")

########
#filter seqtab.nochim
flist<- filterfun(kOverA(5, 2e-05))
seqtab.nochim_filt = filter_taxa(ps, flist, prune=TRUE)
#write_phyloseq(seqtab.nochim_filt, type = "all", path = getwd())
otus <- otu_table(seqtab.nochim_filt)
write.csv(otus, file='/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/otus_filt.csv', quote = F)
#write out as 
otus=read.csv('/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demult_python/otus_filt.csv', header=T, row.names = 1)
all_seq=seqtab.nochim #getSequences(as.matrix(otus))


#hm_sphingo
sphingo="TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCTCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGAGCTAGCGTTATTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAACTCCAGAATTGCCTTTAAGACTGCATCGCTTGAATCCAGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGACTGGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"
hm_sphingo=hm[,sphingo]
lat_sphing=lat_long[names(hm_sphingo),]
sphin_keep=cbind(lat_sphing, hm_sphingo)
p = ggplot(data = sphin_keep, aes(x = Latitude, y = hm_sphingo*100)) + 
  geom_point(color='gray') +
  geom_smooth(method = "lm", se = FALSE, color=1) + theme_bw() + ylab("Percent  Relative Abundance")

pdf("~/Desktop/sphingo_lat.pdf", height=5, width=7,  useDingbats=FALSE)
lb1 <- paste("R == ", round(-0.1485174,2))
p+annotate("text", label = "Sphingomonas\n ", x = 60, y = 5, size = 4, fontface="italic")+annotate("text", label = "Sphingomonas\n ", x = 60, y = 5, size = 4, fontface="italic") + annotate("text", x=60, y=4.7, label=lb1, parse=TRUE)
dev.off()

p_sp_long=ggplot(data = sphin_keep, aes(x = Longitude, y = hm_sphingo*100)) + 
  geom_point(color='gray') +
  geom_smooth(method = "lm", se = FALSE, color=1) + theme_bw() + ylab("Percent  Relative Abundance")

pdf("~/Desktop/sphingo_long.pdf", height=5, width=7,  useDingbats=FALSE)
cor_long_sp=cor.test(sphin_keep$Longitude, sphin_keep$hm_sphingo)[4]$estimate
lb1 <- paste("R == ", round(cor_long_sp,2))
p_sp_long+annotate("text", label = "Sphingomonas\n ", x = 40, y = 5, size = 4, fontface="italic")+ annotate("text", x=40, y=4.8, label=lb1, parse=TRUE)
dev.off()




#beij
beij="TGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTGTCCGGGACGATAATGACGGTACCGGAAGAATAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCCATTCAAGTCGGGGGTGAAAGCCTGTGGCTCAACCACAGAATTGCCTTCGATACTGTTTGGCTTGAGTTTGGTAGAGGTTGGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCCAACTGGACCAATACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG"
hm_beij=hm[,beij]
lat_beij=lat_long[names(hm_beij),]
beij_keep=cbind(lat_beij, hm_beij)

p2 = ggplot(data = beij_keep, aes(x = Latitude, y = hm_beij*100)) + 
  geom_point(color='gray') +
  geom_smooth(method = "lm", se = FALSE, color=1) + theme_bw() + ylab("Percent Relative Abundance")

pdf("~/Desktop/methylo_lat.pdf",  height=5, width=7,  useDingbats=FALSE)
lb1 <- paste("R == ", round(-0.1791008,2))
p2+annotate("text", label = "Methylobacteria", x = 60, y = 7.6, size = 4, fontface="italic") + annotate("text", x=60, y=7.1, label=lb1, parse=TRUE)
dev.off()

p_meth_long=ggplot(data = meth_keep, aes(x = Longitude, y = hm_sphingo*100)) + 
  geom_point(color='gray') +
  geom_smooth(method = "lm", se = FALSE, color=1) + theme_bw() + ylab("Percent  Relative Abundance")

pdf("~/Desktop/meth_long.pdf", height=5, width=7,  useDingbats=FALSE)
cor_long_met=cor.test(beij_keep$Longitude, beij_keep$hm_beij)[4]$estimate
lb1 <- paste("R == ", round(cor_long_met,2))
p_sp_long+annotate("text", label = "Methylobacteria\n ", x = 40, y = 5, size = 4, fontface="italic")+ annotate("text", x=40, y=4.8, label=lb1, parse=TRUE)
dev.off()



