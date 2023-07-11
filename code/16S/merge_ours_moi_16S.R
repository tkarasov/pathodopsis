library(phyloseq)
library(dplyr)
library(ggplot2)
library(dada2)
library(tidyr)
library(DESeq2)
library(reshape2)
library(cowplot)
library(vegan)


# This script takes the seqtabs from our whole experiment and the Carnegie experiment from Moi's lab in 2023 and looks for overlap  in ASV's
####################################
#Load data
####################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
my_phylo <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/seqtab_final.rds")
my_phylo_tax <- readRDS('/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/tax_final.rds')
moi <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_moi_5_2023/processed_reads/seqtab_final.rds")
moi_tax <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_moi_5_2023/processed_reads/tax_final.rds")

#basic manipulation of data
moi_metadata <- read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_moi_5_2023/sample_info/meta_6_20_2024_karasov_moi_leaf_microbe_ids_mergedd.tsv", header=T, sep = "\t")

#there are several NA rows in the metadata sheet. We need to get rid of these. Let's just delete them
hm <- paste("samp_", moi_metadata$Sample_ID, sep ="")
hm <- paste(hm, moi_metadata$tray, sep="_")
hm <- paste(hm, moi_metadata$tray_position, sep = "_")
hm <- paste(hm, moi_metadata$Actual_barcode, sep = "_")
moi_metadata$samp_rename <- hm 

#find duplicate rows and delete one of them
dup <- which(duplicated(moi_metadata$samp_rename))
#'%!in%' <- function(x,y)!('%in%'(x,y))
moi_metadata <- moi_metadata %>% filter(row_number() %!in% dup)
#hm <- hm % filter( %!in% dup)
rownames(moi_metadata) <- moi_metadata$samp_rename

####################################
# Make phyloseq object of combined datasets
####################################
metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_2_2020_reads.tsv",
                      header=T, sep="\t", fill =TRUE)
# new_metadata <- data.frame(sample = rownames(st.all))
# samps <- as.character(new_metadata$sample)
# names(samps) <- c(samps)
# new_metadata$ours <-startsWith(samps, "PA")
# new_metadata$kemen <- startsWith(samps, "SRR")
keep = data.frame(OTU_clim$refseq)[,1]

####################################
# Subset to 1000 reads and make final table
####################################
st.phylo <- phyloseq(otu_table(moi, taxa_are_rows = FALSE), sample_data(moi_metadata))
st.phyo <- rarefy_even_depth(st.phylo, sample.size = 1000, rngseed = 4)
colnames(plant_clim$otu_table) <- data.frame(OTU_clim$refseq)[,1]
rownames(plant_clim$clim_data) <- plant_clim$clim_data$Sequence_ID
plant_phyl <- phyloseq(otu_table(plant_clim$otu_table, taxa_are_rows = TRUE), sample_data(plant_clim$clim_data))
fin.moi <- phyloseq::prune_taxa(keep, st.phyo)
fin.ours <- prune_taxa(colnames(otu_table(fin.moi)), plant_phyl)

# This is still too many ASVs (more ASVs than samples)

####################################
# Convert phyloseq to vegan
####################################
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


# We cannot reasonably test teh effect of genotype given that there is low replication (109 genotypes and 171 samples from Moi's experiment that pass the filter of enough reads)
# instead let's see which of these ASVs are influenced by the environmental treatment. Then see which of these ASVs shows a latitudinal gradient. 











####################################
# And Let's plot the MDS with the clusters
####################################
rop <-rowSums(otu_table(fin.ours))
val = 1000
MDS_mine <- sqrt(otu_table(fin.ours)/1000) %>% dist() %>% cmdscale() #(eig = TRUE,  k = 3) #)dim(plant_clim$otu_table)[1]-1))
MDS_mine_tot <- sqrt(otu_table(fin.ours)/1000) %>% dist() %>% cmdscale(eig = TRUE,  k = 3)

# Solve for betas in pcoa
x1 <- cbind(1, sqrt(otu_table(fin.ours)/1000)) # add intercept
B <- solve(t(x1) %*% x1) %*% t(x1) %*% MDS_mine # Betas

# Project kemen data into pcoa space 
kemen2 <- cbind(1, (sqrt(otu_table(fin.kemen)/1000)))
new_kemen <- data.frame(kemen2 %*% B)
colnames(new_kemen) <- c("eig1", "eig2")
new_ours <- data.frame(x1 %*% B)
colnames(new_ours) <- c("eig1", "eig2")

# rename points
MDS.points <- data.frame(MDS_mine_tot$points)
colnames(MDS.points)[c(1:2)] = c("MDS1", "MDS2")

#calculate percentage explained
exp3 <-  ((MDS_mine_tot$eig) / sum(MDS_mine_tot$eig))[1]*100
exp4 <-  ((MDS_mine_tot$eig) / sum(MDS_mine_tot$eig))[2]*100

# plot basic kmeans
MDS_plot_kmeans <- 
  ggplot(data = MDS.points, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(col=as.factor(sample_data(fin.ours)$cluster)), cex = 3, alpha = 0.1) +
  scale_colour_manual(name = "Cluster", values = c("#D95F02", "#1B9E77", "#E31A1C", "#1F78B4")) +
  theme_bw() +
  xlab(paste(paste("MDS1 (", round(exp3), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp4), sep=""),"%)",sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey70") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey70") +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kemen_pcoa.pdf",useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3)
proj <- MDS_plot_kmeans + 
  geom_point(data = new_kemen, aes(x=eig1, y=eig2,col = as.factor(sample_data(fin.kemen)$Season)), cex = 2) 
  #geom_point(data=MDS.points[as.character(germans$Sequence_ID),], aes(x=MDS1, y=MDS2))
  #geom_point(data = new_ours, aes(x=eig1, y=eig2) ) 
dev.off()

####################################
# pcoa of everything together
####################################
all_otus = rbind(otu_table(fin.kemen), otu_table(fin.ours))
run_pca <- data.frame(sqrt(all_otus/1000)) %>% dist() %>% cmdscale(eig = TRUE,  k = 3)
flat_pca <-data.frame(sqrt(all_otus/1000) %>% dist() %>% cmdscale())
flat_pca$col <- "Pathodopsis"
flat_pca[1:dim(otu_table(fin.kemen))[1],]$col <- "Kemen"
flat_pca$Season <- "Spring"
flat_pca[1:dim(otu_table(fin.kemen))[1],]$Season <- sample_data(fin.kemen)$Season

# rename points
all.points <- data.frame(run_pca$points)
colnames(all.points)[c(1:2)] = c("MDS1", "MDS2")

#calculate percentage explained
exp3 <-  ((run_pca$eig) / sum(run_pca$eig))[1]*100
exp4 <-  ((run_pca$eig) / sum(run_pca$eig))[2]*100

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kemen_total_divided_pcoa.pdf",useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3)
ggplot(data=flat_pca, aes(x=X1, y=X2, col=col)) +
  geom_point(data=flat_pca, aes(shape=Season), cex = 2) +
  #scale_colour_brewer(name = "Cluster", palette = "Dark2") +
  theme_bw() +
  scale_color_viridis_d() + 
  xlab(paste(paste("MDS1 (", round(exp3), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp4), sep=""),"%)",sep="")) +
  #geom_hline(yintercept=0, linetype="dashed", color = "grey70") +
  #geom_vline(xintercept=0, linetype="dashed", color = "grey70") +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )
dev.off()

####################################
# Placement as a reflection of season
####################################
all_kemen <- cbind(new_kemen, sample_data(fin.kemen))
effect_1 <- lm(eig1 ~ Season + Year + lat_lon, data = all_kemen)
effect_2 <- lm(eig2 ~ Season , data = all_kemen)


####################################
# Differential expression of genes between season and by cluster
####################################
metadata <- metadata[-which(duplicated(metadata$Sequence_ID)),] 
metadata2 <- metadata[1:1074,]#metadata[-which(metadata$Sequence_ID=="NA"),]
rownames(metadata2) <- metadata2$Sequence_ID

subset_kemen <- otu_table(st.phylo)[,which(colnames(otu_table(st.phylo)) %in% keep)]
subset_mine <- my_phylo[,which(colnames(otu_table(st.phylo)) %in% keep)]

#Kemen
phylo_kemen <- phyloseq(otu_table(subset_kemen, taxa_are_rows = FALSE)+1, sample_data(st.phylo))
diagdds.kemen = phyloseq_to_deseq2(phylo_kemen, ~Season)
diagdds.kemen = DESeq(diagdds.kemen, test="Wald", fitType="parametric")

#Mine
plant_clim$clim_data$cluster <- as.factor(plant_clim$clim_data$cluster)
my_total <-phyloseq(otu_table(subset_mine, taxa_are_rows = FALSE) + 1, 
                    tax_table(my_phylo_tax), sample_data(plant_clim$clim_data))
diagdds.ours = phyloseq_to_deseq2(my_total, ~1+(cluster))
diagdds.ours = DESeq(diagdds.ours, test="Wald", fitType="parametric")

# Identify ASVs that differ in abundance between clusters 1 and 2
results.kemen <- results(diagdds.kemen, name="Season_Winter_vs_Spring")
results.ours <- results(diagdds.ours, name="cluster_2_vs_1")

# Identify how many of them change with the season
results_df <- data.frame(gene=results.ours@rownames, 
                         kemen=results.kemen[results.ours@rownames,]$log2FoldChange,
                         ours_1_3=results.ours[results.ours@rownames,]$log2FoldChange,
                         pval_kemen=results.kemen[results.ours@rownames,]$padj,
                         pval_ours=results.ours[results.ours@rownames,]$padj)
results_df$joint_pval <- as.numeric(as.character(apply(results_df, 1, function(x) max(x[4:5]))))
results_df$joint_pval <- results_df$joint_pval<0.05
results_df <- results_df[which(is.na(results_df$gene)==FALSE),]


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/FC_kemen_ours.pdf",useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3)
FC_kemen <- ggplot(data=results_df, aes(x=kemen, y=ours_1_3, col=joint_pval))+
  geom_point() +
  xlab("log2(FC Seasons)") +
  ylab("log2(FC Clusters)") +
  theme_bw() +
  scale_color_viridis_d()+
  geom_hline(yintercept = 0, lty="dashed")+
  geom_vline(xintercept = 0, lty="dashed")+
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.8),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
FC_kemen 
 
dev.off()

kmeans <- MDS_plot_kmeans #readRDS( "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/MDS_plot_kmeans.rds")
kmeans_map <- readRDS( "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kmeans_map.rds")
all_MDS<- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/MDS_all.rds")
thal_cap <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/thal_cap_mds.rds" )

poop <- plot_grid(all_MDS, thal_cap)
p1 <- plot_grid( all_MDS + theme(legend.position = "none"), thal_cap + theme(legend.position = "none"),
                 kmeans, FC_kemen + xlab("log2(FC Seasons") + ylab("log2(FC Clusters"), align = 'hv', ncol = 2)

p2 <- plot_grid(kmeans_map, proj, rel_widt = c(0.7, 0.3), align = 'hv', ncol =2)
fig2 <- plot_grid(p1, p2,
          #proj, FC_kemen, 
          align = 'hv', nrow = 2)
  
# pdf(fig2,"/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/fig2_all.pdf", useDingbats = FALSE, 
#    font = "ArialMT", width = 3.5)


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kemen_comparisons_fig2.pdf", useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3.5/2)
plot_grid(proj, FC_kemen)
dev.off()


####################################
# Seasonal Fluctuations of Kemen microbiome
####################################
kemen_sig <- results_df %>% filter(pval_kemen<0.05) %>% select(gene) 
ot <- data.frame(otu_table(phylo_kemen))
phylo_sig <- apply(ot, 2, function(i)i/sum(i))
phylo_sig <- data.frame(phylo_sig[, which(as.character(kemen_sig[,1])%in%colnames(otu_table(phylo_kemen)))])

#Subset to 20
phylo_sig <- phylo_sig[,1:10]

#phylo_sig <- phylo_sig
phylo_sig$Year <-as.numeric(as.character(sample_data(phylo_kemen)$Year))
phylo_sig$Season <- sample_data(phylo_kemen)$Season



phylo_sig[which(phylo_sig$Season=="Winter"),]$Season<-.12
phylo_sig[which(phylo_sig$Season=="Spring"),]$Season<-.04
phylo_sig$Date <- as.factor(as.character(paste(phylo_sig$Year, phylo_sig$Season, sep="")))
melted_sig <- melt(phylo_sig,id.vars="Date")
melted_sig <- melted_sig%>%filter(variable!="Year")%>%filter(variable!="Season")
melted_sig$value <- as.numeric(as.character(melted_sig$value))

df2 <- aggregate(value ~ Date + variable, melted_sig, mean)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kemen_sig_fluctuations.pdf",useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3)
ggplot(df2, aes(x=Date, y=as.numeric(as.character(value)), group = variable, col = variable))+
  geom_line() +
  scale_color_viridis_d() +
  theme_bw() +
  ylab("Relative Abundance") +
  theme(legend.position="none",
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

####################################
# Count per quadrant
####################################
pos_pos <- dim(results_df %>% filter(kemen >= 0) %>% filter(ours_1_3 >=0))[1]
pos_neg <- dim(results_df %>% filter(kemen >= 0) %>% filter(ours_1_3 <0))[1]
neg_pos <- dim(results_df %>% filter(kemen < 0) %>% filter(ours_1_3 >=0))[1]
neg_neg <- dim(results_df %>% filter(kemen < 0) %>% filter(ours_1_3 < 0))[1]
tot <- pos_pos + pos_neg + neg_pos + neg_neg
sub.freqs <- c(pos_pos, pos_neg, neg_pos, neg_neg)

# Calculate the multinomial probability confidence intervals (alpha = 0.05, two-sided)
multinom_prob <- MultinomCI(sub.freqs, conf.level = 0.99, method = "wald")

