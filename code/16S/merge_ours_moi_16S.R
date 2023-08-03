library(phyloseq)
library(dplyr)
library(ggplot2)
library(dada2)
library(tidyr)
library(DESeq2)
library(reshape2)
library(cowplot)
library(vegan)
library(microbiome)
library(PathoStat)


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
moi_metadata <- read.csv("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_moi_5_2023/sample_info/meta_7_19_2023_karasov_moi_leaf_microbe_ids_merged_drought.tsv", header=T, sep = "\t")

#there are several NA rows in the metadata sheet. We need to get rid of these. Let's just delete them
hm <- paste("samp_", moi_metadata$Sample_ID, sep ="")
hm <- paste(hm, moi_metadata$tray, sep="_")
hm <- paste(hm, moi_metadata$tray_position, sep = "_")
hm <- paste(hm, moi_metadata$Actual_barcode, sep = "_")
moi_metadata$samp_rename <- hm 

#find duplicate rows and delete one of them
dup <- which(duplicated(moi_metadata$samp_rename))
'%!in%' <- function(x,y)!('%in%'(x,y))
moi_metadata <- moi_metadata %>% filter(row_number() %!in% dup)
#hm <- hm % filter( %!in% dup)
rownames(moi_metadata) <- moi_metadata$samp_rename

####################################
# Make phyloseq object of combined datasets
####################################
metadata = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_2_2020_reads.tsv",
                      header=T, sep="\t", fill =TRUE)
#keep the ASVs specified in OTU_clim
keep = data.frame(OTU_clim$refseq)[,1]

#here are the seqid names
keep_seqid = OTU_clim$refseq@ranges@NAMES

####################################
# Subset Moi's data to samples with at least 1000 reads and ASVs found in at least 10% of samples
####################################
st.phylo <- phyloseq(otu_table(moi, taxa_are_rows = FALSE), sample_data(moi_metadata))
st.phyo <-   prune_samples(sample_sums(st.phylo) >= 1000, st.phylo)

# prune also for those that are found in at least 10% of samples
who <- genefilter_sample(st.phyo, filterfun_sample(function(x)x > 5), A=0.1*nsamples(st.phyo))
                                  
who_keep <- names(which(who==TRUE))
st.phylo <- prune_taxa( who_keep, st.phyo)
colnames(plant_clim$otu_table) <- data.frame(OTU_clim$refseq)[,1]
rownames(plant_clim$clim_data) <- plant_clim$clim_data$Sequence_ID
plant_phyl <- phyloseq(otu_table(t(plant_clim$otu_table), taxa_are_rows = TRUE),  sample_data(plant_clim$clim_data))
tax_table(plant_phyl) <- tax_table(matrix(plant_clim$tax_table))
fin.moi <- phyloseq::prune_taxa(keep, st.phylo) #
the_common_set <- colnames(otu_table(fin.moi))
fin.ours <- prune_taxa(the_common_set, plant_phyl) #plant_phyl has 570 taxa and 454 samples
#rownames(otu_table(fin.ours)) <- keep_seqid


# We cannot reasonably test teh effect of genotype given that there is low replication (109 genotypes and 171 samples from Moi's experiment that pass the filter of enough reads)
# instead let's see which of these ASVs are influenced by the environmental treatment. Then see which of these ASVs shows a latitudinal gradient. It's not clear to me what the zone variable is telling us but let's see...

####################################
# Differential abundance  of ASVs between treatments
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#wald-test-individual-steps
####################################
metadata <- metadata[-which(duplicated(metadata$Sequence_ID)),] 
metadata2 <- metadata[1:1074,]#metadata[-which(metadata$Sequence_ID=="NA"),]
rownames(metadata2) <- metadata2$Sequence_ID


# Get subsetted otu table 
subset_moi <- otu_table(fin.moi) #otu_table(st.phylo)[,which(colnames(otu_table(st.phylo)) %in% keep)]
subset_mine <- otu_table(fin.ours)
                                                   
#Moi's phyloseq object
phylo_moi <- phyloseq(otu_table(subset_moi, taxa_are_rows = FALSE)+1, sample_data(st.phylo))

#this model just tests drought treatment
diagdds.moi = phyloseq_to_deseq2(phylo_moi, ~1+(treatment.x))
diagdds.moi = DESeq(diagdds.moi, test="Wald", fitType="parametric")

#this model just tests PRS
diagdds.moi_prs = phyloseq_to_deseq2(phylo_moi, ~1+(PRS))
diagdds.moi_prs = DESeq(diagdds.moi_prs, test="Wald", fitType="parametric")

#this model test treatment and PRS and interaction
diagdds.moi_joint = phyloseq_to_deseq2(phylo_moi, ~1+(PRS) + (treatment.x) + PRS:treatment.x)
#diagdds.moi_joint = DESeq(diagdds.moi_joint, test="LRT", fitType="parametric", reduced=~1)
diagdds.moi_joint = DESeq(diagdds.moi_joint, test="Wald", fitType="parametric")

#My all data phyloseq object
sample_data(fin.ours)$cluster <- as.factor(sample_data(fin.ours)$cluster)
my_total <-phyloseq(otu_table(subset_mine, taxa_are_rows = FALSE) + 1, 
                   sample_data(fin.ours))
#sample_data(my_total)$cluster <- sample(sample_data(my_total)$cluster)
diagdds.ours = phyloseq_to_deseq2(my_total, ~1+(cluster))
diagdds.ours = DESeq(diagdds.ours, test="Wald", fitType="parametric")

# Identify ASVs that differ in abundance between clusters 1 and 2
results.moi <- results(diagdds.moi, name="treatment.x_Watered_vs_Drought")
results.moi_prs <- results(diagdds.moi_prs, name="PRS_lps_vs_hps")
#results.moi_joint <- results(diagdds.moi_joint, name="PRSlps.treatment.xWatered") #20% of the ASVs (17/85) had a different response to drought depending on which genotype they were. 
results.moi_joint <- results(diagdds.moi_joint, contrast = list("PRS_lps_vs_hps", "PRSlps.treatment.xWatered")) # 21/103 had a significant contrast in this. 
results.ours <- results(diagdds.ours, name="cluster_2_vs_1") #table(results.ours$padj<0.01) 66 FALSE, 88 TRUE. When this is done with randomized cluster assignments get 107 True 47 False. After filtering to the main 20 ASVs this is 2 sig/20.
# 17 significantly different, 3 not for the Europe dataset. When it is randomized it is 16 No, 4 yes. I suppose I need an empirical estimate of the enrichment in differences. 

# What percentage of the ASVs that are significant between clusters have a significant interaction term


# Identify how many of them change with the treatment
results_df <- data.frame(gene=results.ours@rownames, 
                         moi=results.moi_prs[results.ours@rownames,]$log2FoldChange,
                         ours_1_3=results.ours[results.ours@rownames,]$log2FoldChange,
                         pval_treat = results.moi[results.ours@rownames,]$padj,
                         pval_joint=results.moi_joint[results.ours@rownames,]$padj,
                         pval_moi_prs=results.moi_prs[results.ours@rownames,]$padj,
                         pval_ours=results.ours[results.ours@rownames,]$padj)


# Now I want to put the taxonomy info on the 20 genes

asv_taxonomy <- data.frame(cbind(t(plant_clim$otu_table), plant_clim$tax_table))
asv_taxonomy$seqID <- rownames(plant_clim$tax_table)
results_df$seqID <- asv_taxonomy[results_df[,1],]$seqID
results_df$Genus <- asv_taxonomy[results_df[,1],]$Genus
write.table(results_df, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/moi_data.txt", sep= "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
                                                   

################
# Visualization of Moi data with pheatmap
################

dds = diagdds.moi_joint
rownames(dds) <- results_df$seqID
rld <- rlog(dds)
sampleDists <- dist( t( assay(rld) ) )
mat <- assay(rld)
mat <- mat - rowMeans(mat)
rownames(mat) <- results_df$seqID
mat_sub <- mat[sig_hps_lps,]
df <- as.data.frame(colData(rld)[,c("PRS","treatment.x")])
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/pheatmap_moi_comparisons.pdf", useDingbats = FALSE, 
    font = "ArialMT")
pheatmap(mat, scale = 'column', cluster_rows = F, annotation_col=df)
dev.off()           





sig_hps_lps <-  results_df$seqID[which (results_df$pval_moi_prs<=0.05)]
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/plotCounts.pdf", useDingbats = FALSE, 
    font = "ArialMT")
plotCounts(dds, gene="seq_12", intgroup=c("PRS"))
dev.off()

                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
###############
# Below code was used for sanity checkin
###############                                                   


                                                   
####I need to try a randomization of prs. This next code is checking the legitimacy
phylo_moi_random <- phylo_moi
set.seed(44452962) 
prs_ran <- sample(sample_data(phylo_moi)$PRS)
sample_data(phylo_moi_random)$PRS <- prs_ran
diagdds.moi_resamp = phyloseq_to_deseq2(phylo_moi_random, ~1+(PRS))
diagdds.moi_resamp = DESeq(diagdds.moi_resamp,test="Wald", fitType="parametric")           

diagdds.moi_score = phyloseq_to_deseq2(phylo_moi, ~1+(PRS))
diagdds.moi_score = DESeq(diagdds.moi_score,test="Wald", fitType="parametric")         
results.moi_resamp <- results(diagdds.moi_resamp, name="PRS_lps_vs_hps")
results.moi_score <- results(diagdds.moi_score, name="PRS_lps_vs_hps")
table(results.moi_resamp$padj<0.01)
table(results.moi_score$padj<0.01)
# Okay, maybe we are making progress. 20 ASVs pass all of the filtering. Two are significant 
                                          
#this analysis prior to filtering was really worrisome, because frequently about 1/4 of the genes were significant when PRS is permuted.  This blog suggests that many of these ASVs might not adhere to the negative binomial well. After filtering there were no significant ASVs in the permutation.
# https://towardsdatascience.com/deseq2-and-edger-should-no-longer-be-the-default-choice-for-large-sample-differential-gene-8fdf008deae9
dge = phyloseq_to_edgeR(phylo_moi, "PRS")
# Perform binary test
et = exactTest(dge)
# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.01 res
                                                   
sigtab = res[(res$FDR < alpha), ]
#sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), ], "matrix"))
dim(sigtab)
                                          
# Let's try Wilcoxan test read data
 #otu_table(st.phylo)[,which(colnames(otu_table(st.phylo)) %in% keep)]
readCount <- t(otu_table(phylo_moi))
conditions <- sample_data(phylo_moi)[,c("PRS")]
conditions <- factor(t(conditions))
# edgeR TMM normalize
y <- DGEList(counts = readCount, group = conditions)
## Remove rows conssitently have zero or very low counts
#keep <- which(y$samples$lib.size>=500)#
keep <- filterByExpr(y, min.count = 5000 )
#y <- y[keep, ]
y <- y[keep,keep.lib.sizes = FALSE]
## Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
# Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm), function(i){
data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
p <- wilcox.test(gene~conditions, data)$p.value
return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")
# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
write.table(outRst[outRst$FDR<fdrThres,], file = "examples/examples.WilcoxonTest.rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)
                                          

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/FC_moi_ours.pdf",useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3)
FC_moi <- ggplot(data=results_df, aes(x=moi, y=ours_1_3, col=joint_pval))+
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
FC_moi
 
dev.off()

####################################
# Are the same SNPs different between clusters and between droughts. 154
####################################
pval_ours <- results.ours$padj
pval_moi <- results.moi$padj
#this didn't work well. How about the relationship between latitude and the Fold Change of the Moi associated SNPs. Actually it worked better than I thought...

keep_dat <- sample_data(my_total)
keep_otu <- t(otu_table(my_total))
#keep_otu <- keep_otu[which(keep_dat$Lat<45 & keep_dat$Lat>30),]
                                                       
asv_relation <- function(asv){
  # rownames(keep_otu) <- keep_dat$id
  model1 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat[rownames(keep_otu),]$Lat)
  model2 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat[rownames(keep_otu),]$PDSI)
  model3 <- lm(as.numeric(keep_otu[,asv]) ~ keep_dat[rownames(keep_otu),]$PDSI + keep_dat[rownames(keep_otu),]$Lat)
  pval_lat <- summary(model1)$coefficients[2,4]
  pval_pds <- summary(model2)$coefficients[2,4]
  pval_lat_tog <- summary(model3)$coefficients[3,4]
  pval_pds_tog <- summary(model3)$coefficients[2,4]
  return(c(p.adjust(pval_lat, "BH"), p.adjust(pval_pds, "BH"), p.adjust(pval_lat_tog, "BH"), p.adjust(pval_pds_tog, "BH")))
}
N = length(colnames(keep_otu))
otu_sig <- data.frame(pval_lat = numeric(N), 
    pval_pds = numeric(N),
    pds_tog = numeric(N),
    lat_tog = numeric(N),
    lat_corr = numeric(N))
rownames(otu_sig) = colnames(keep_otu)

for(asv in colnames(keep_otu)){
  print(asv)
  pvals <- asv_relation(asv)
  otu_sig[asv,c(1:4)] = pvals
  otu_sig[asv, 5] = cor.test(as.numeric(keep_otu[,asv]),  keep_dat[rownames(keep_otu),]$Lat)$estimate
}


# Now let's look at the otu_sig and the ASVs that are shared.
keep_otu_sig <- otu_sig[rownames(results.moi),]
                                                       
# for some reason only about 10% of the ASVs are correlated in lm with latitude. Why so low?                                                       
#> table(otu_sig[,2]<0.01)

#FALSE  TRUE 
#  141    13 
                                                       
                                                       
                                                       
                                                       
                                                       



                                                       

kmeans <- MDS_plot_kmeans #readRDS( "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/MDS_plot_kmeans.rds")
kmeans_map <- readRDS( "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kmeans_map.rds")
all_MDS<- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/MDS_all.rds")
thal_cap <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/thal_cap_mds.rds" )

poop <- plot_grid(all_MDS, thal_cap)
p1 <- plot_grid( all_MDS + theme(legend.position = "none"), thal_cap + theme(legend.position = "none"),
                 kmeans, FC_moi + xlab("log2(FC Seasons") + ylab("log2(FC Clusters"), align = 'hv', ncol = 2)

p2 <- plot_grid(kmeans_map, proj, rel_widt = c(0.7, 0.3), align = 'hv', ncol =2)
fig2 <- plot_grid(p1, p2,
          #proj, FC_moi, 
          align = 'hv', nrow = 2)
  
# pdf(fig2,"/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/fig2_all.pdf", useDingbats = FALSE, 
#    font = "ArialMT", width = 3.5)


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/moi_comparisons_fig2.pdf", useDingbats = FALSE, 
    font = "ArialMT", width = 3.5, height  = 3.5/2)
plot_grid(proj, FC_moi)
dev.off()


