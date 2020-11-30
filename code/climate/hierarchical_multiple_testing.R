#!/usr/bin/env Rscript

library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)
library(ade4)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(genefilter)
library(DESeq2)
library(BiocParallel)
library(edgeR)
library(viridis)
#the goal of this script is to determine OTUs differentially abundant between  (1) A. thalaina and Capsella and differentially abundant between (2) Environmental Conditions.

#Much of the analysis comes from here: https://f1000research.com/articles/5-1492/v2

#######################################################################
# Read in Phyloseq object
#######################################################################
# Load the phyloseq object generated in after_dada2_make_otutable
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/GP50_subset_fin.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
#load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")



# #######################################################################
# # Prune phyloseq to only those samples with at least 1000 reads, and to those ASVs that are zero add 1
# #######################################################################
#subset to only my phylotypes of interest
GPmine <- prune_taxa(colnames(otu_table(GP_at15_all)), GP50)


# #######################################################################
# For differential analysis use limma-voom
# #######################################################################
# https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html
m = as(otu_table(GPmine), "matrix")

# Add one to protect against overflow, log(0) issues.
m = m + 1

# Define gene annotations (`genes`) as tax_table
taxonomy = tax_table(GPmine, errorIfNULL=FALSE)
if( !is.null(taxonomy) ){
  taxonomy = data.frame(as(taxonomy, "matrix"))
} 

# Now turn into a DGEList
d = DGEList(counts=t(m), genes=taxonomy, remove.zeros = TRUE)

# Calculate the normalization factors
z = calcNormFactors(d, method="RLE")

# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(z$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}

#plotMDS(z, col = as.numeric(factor(sample_data(GP_at15_all)$Species)), labels = as.numeric(factor(sample_data(GP_at15_all)$Species)))

# #######################################################################
# Test difference between soil and not soil
# #######################################################################

# Creat a model based on Slash_pile_number and depth
GP_mat <- data.frame(as(sample_data(GPmine), "matrix"))
mm <- model.matrix(~ 0 + Host_Species + Tour_ID, data= GP_mat) # specify model with no intercept for easier contrasts

#estimate dispersions
# , edgeR uses the Cox-Reid profile-adjusted
# likelihood (CR) method in estimating dispersions [22]. The CR method is derived to overcome
# the limitations of the qCML method as mentioned above. It takes care of multiple factors by
# fitting generalized linear models (GLM) with a design matrix.
dge <- estimateGLMCommonDisp(z, mm)
dge <- estimateGLMTagwiseDisp(dge, mm)
dge <- estimateGLMTrendedDisp(dge, mm)
plotBCV(dge)

# y <- voom(z, mm, plot = T)
# fit <- lmFit(dge, mm)
# head(coef(fit))
fit <- glmQLFit(dge, mm, robust = TRUE)
plotQLDisp(fit)
lrt_soil <- glmQLFTest( fit, contrast = c(1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
is.de <- decideTestsDGE(lrt_soil)
plotMD(lrt_soil, status=is.de, values =c(1,-1), col=c("red","blue"),legend="topright")
tmp2 <- topTags(lrt_soil, n = 1000000)$table
sigtab = cbind(as(tmp2, "data.frame"), as(tax_table(GPmine)[rownames(tmp2), ], "matrix"))
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen_soil <- sigtabgen
# # Phylum order
# x = tapply(sigtabgen_soil$logFC, sigtabgen_soil$Phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtabgen_soil$Phylum = factor(as.character(sigtabgen_soil$Phylum), levels = names(x))
# # # Genus order
# x = tapply(sigtabgen_soil$logFC, sigtabgen_soil$Genus, function(x) mean(x))
# x = sort(x, TRUE)
# sigtabgen_soil$Genus = factor(as.character(sigtabgen_soil$Genus), levels = names(x))

# Now prepare for plotting
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
colors <- viridis_pal()(10)
color_phylum_map <-viridis_pal()(length(levels(sigtabgen_soil$Phylum)))
names(color_phylum_map) <- levels(sigtabgen_soil$Phylum)

phylum_levels = levels(sigtabgen_soil$Phylum)



sigtabgen_soil$coloring <- color_phylum_map[sigtabgen_soil$Phylum]
sigtabgen_soil$sig <- sigtabgen_soil$FDR<0.01
sigtabgen_soil$sig[sigtabgen_soil$sig==TRUE]= 1
sigtabgen_soil$sig[sigtabgen_soil$sig==FALSE] = 0.25
#sigtabgen_soil <- sigtabgen_soi %>% filter(FDR<=0.01)
diff_soil_non_soil <- 
  ggplot(sigtabgen_soil, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=3, aes(alpha = sig))  + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + scale_color_viridis_d() + 
  theme(legend.position="bottom") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")


# #######################################################################
# Test difference thaliana and capsella
# #######################################################################
GP_plant <- prune_samples(sample_data(GPmine)$Sample_type == "M", GPmine)
good_pops <- sample_data(GP_plant)[which(sample_data(GP_plant)$Species == "Cap"),]$SiteID
cap_has <- prune_samples(sample_data(GP_plant)$SiteID %in%good_pops, GP_plant)
GP_plant <- cap_has
m.plant = otu_table(GP_plant)
d.plant = DGEList(counts=t(m.plant), genes=taxonomy, remove.zeros = TRUE)
z.plant = calcNormFactors(d.plant, method="TMM")

# Creat a model based on Slash_pile_number and depth
mm.random <- model.matrix(~ 0 + Species + TourID, data=data.frame(as(sample_data(GP_plant),"matrix"))) # specify model with no intercept for easier contrasts

#estimate dispersions
# , edgeR uses the Cox-Reid profile-adjusted
# likelihood (CR) method in estimating dispersions [22]. The CR method is derived to overcome
# the limitations of the qCML method as mentioned above. It takes care of multiple factors by
# fitting generalized linear models (GLM) with a design matrix.
dge <- estimateGLMCommonDisp(z.plant, mm.random)
dge <- estimateGLMTagwiseDisp(dge, mm.random)
dge <- estimateGLMTrendedDisp(dge, mm.random)
plotBCV(dge)
fit <- glmQLFit(dge, mm.random, robust = TRUE)
plotQLDisp(fit)
thal.caps <- makeContrasts(SpeciesAth-SpeciesCap, levels=mm.random)
lrt <- glmQLFTest( fit, contrast = thal.caps)#contrast = c(0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
is.de <- decideTestsDGE(lrt)
plotMD(lrt, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")
tmp3 <- topTags(lrt, n = 1000000)$table
sigtab = cbind(as(tmp3, "data.frame"), as(tax_table(GPmine)[rownames(tmp2), ], "matrix"))
theme_set(theme_bw())

sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = levels(sigtabgen_soil$Phylum))
sigtabgen$coloring <- color_phylum_map[sigtabgen$Phylum]
sigtabgen$sig <- sigtabgen$FDR<0.01
sigtabgen$sig[sigtabgen$sig==TRUE]= 1
sigtabgen$sig[sigtabgen$sig==FALSE] = 0.25
#sigtabgen <- sigtabgen %>% filter(FDR<=0.01)

diff_cap_ath <- ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=3, aes(alpha = sig)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + scale_color_viridis_d() + theme(legend.position="bottom") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grobs <- ggplotGrob(diff_cap_ath)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/diff_abundance_soil_capsella.pdf", height = 12, width = 15,
    useDingbats = FALSE, font = "ArialMT")
plot_grid(diff_soil_non_soil + ylab("logFC(Plant/Soil)") +
            theme(legend.position = "NULL"), 
          diff_cap_ath + ylab("logFC(A. thaliana/Other Brassicaceae)") +
          theme(legend.position = "NULL"), legend, rel_heights = c(1,1,.1), nrow = 3)
dev.off()




# # #######################################################################
# # Show a few specific
# # #######################################################################
keep = names(sort(rowSums(vm.plant$E), decreasing = TRUE))
first_seqs <- sigtabgen
first_seqs <- first_seqs[keep,]
first_seqs <- first_seqs[which(first_seqs$FDR < 0.01 & abs(first_seqs$logFC) > 2 ),][1:4,]
seqs.names <-rownames(first_seqs)[1:4]
vm.plant <- voom(d.plant, design = mm.random, plot = TRUE )
cpm <- data.frame(t(vm.plant$E[seqs.names,]))
colnames(cpm) <- paste(seqs.names, first_seqs$Genus, sep="_")
cpm$Species <- as.character(vm.plant$design[,"SpeciesAth"])
cpm[which(cpm$Species=="1"),]$Species = "A. thaliana"
cpm[which(cpm$Species=="0"),]$Species = "Other Brassicaceae"
#cpm[cpm$Species == "1",]$Species = "Ath" 

cpm.melt <- melt(cpm, id = "Species")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/spec_OTU_diff_abundance_soil_capsella.pdf",  width = 7, height = 3.5,
    useDingbats = FALSE, font = "ArialMT")
ggplot(data = cpm.melt, aes(x = Species, y = 2^value/10000, col = Species)) +
  geom_jitter(cex = 1.5, alpha = 0.5, width = 0.2, key_glyph = rectangle_key_glyph(fill = color)) +
  stat_summary(fun.ymin=median, fun.ymax=median, fun.y=median, geom="crossbar") +
  scale_y_continuous(trans = "log10", labels = comma) +
    facet_grid(. ~variable, ) +
  ylab("RA (%)") + scale_color_manual(values = wes_palette(n=3, name="Moonrise2")) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "bottom")
dev.off()

# # #######################################################################
# # Test difference thaliana and capsella using limma-voom
# # #######################################################################
# #vm <- voomWithQualityWeights(d, design=mm)
# #subset to sites in which have data for both
# 
# 
# mm.random <- model.matrix(~ 0 + Species + TourID, data=data.frame(as(sample_data(GP_plant),"matrix")))
# 
# # Now let's do voom which transforms count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.
# vm.plant <- voom(d.plant, design = mm.random, plot = TRUE )
# 
# #The corfit function is why we want to use voom for mixed effect models. We are using this with blocking on siteID
# plant_mat <- data.frame(as(sample_data(GP_plant),"matrix"))$SiteID
# 
# #corfit <- duplicateCorrelation(vm, mm.random, block = plant_mat)
# 
# #And now we use the correlation within block
# fit <- lmFit(vm.plant, mm.random) #m, block = data.frame(as(sample_data(GPmine),"matrix"))$SiteID, correlation=corfit$consensus)
# 
# # Comparison between cultivars A thaliana and Capsella with proper blocking
# contr <- makeContrasts(SpeciesAth - SpeciesCap,
#                        levels = colnames(coef(fit)))
# 
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp, robust = TRUE)
# tmp2 <- topTable(tmp, coef=1, sort.by = "P", n = Inf)
# tmp2$Taxa <- rownames(tmp2)
# tmp2 <- tmp2[,c("Taxa","logFC","AveExpr","P.Value","adj.P.Val")]
# length(which(tmp2$adj.P.Val < 0.01)) # number of DE genes
# 
# #tmp2 <- topTags(lrt, n = 1000000)$table
# sigtab = cbind(as(tmp2, "data.frame"), as(tax_table(GPmine)[rownames(tmp2), ], "matrix"))
# theme_set(theme_bw())
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }
# sigtabgen = subset(sigtab, !is.na(Genus))
# # Phylum order
# x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) mean(x))
# x = sort(x, TRUE)
# sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# sigtabgen$Phylum <- factor(sigtabgen$Phylum, levels = phylum_levels)
# #relevel(sigtabgen$Phylum) <- phylum_le
# # Genus order
# x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) mean(x))
# x = sort(x, TRUE)
# sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
# #sigtabgen <- sigtabgen %>% filter(FDR<=0.01)
# 
# diff_cap_ath_voom <- ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=3) + 
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + scale_color_viridis_d()
# 
# 
# 
# 
#









# #######################################################################
# # Convert Phyloseq to deseq
# #######################################################################
# ps_dds <- phyloseq_to_deseq2(GP_at15_all, ~ PlantID + Sample_type)
# 
# #######################################################################
# # Perform variance stabilizing transformation
# #######################################################################
# vsd <- varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "local")
# 
# 
# #######################################################################
# # Perform variance stabilizing transformation
# #######################################################################
# ps_dds <- estimateSizeFactors(ps_dds)
# ps_dds <- estimateDispersions(ps_dds)
# abund <- getVarianceStabilizedData(ps_dds)
# 
# # perform DEA
# dea_indp <- DESeq(dds_indp)
# dea_paired <- DESeq(ps_dds, parallel=TRUE, BPPARAM=MulticoreParam(8), quiet = FALSE)
# 
# # result analysis
# res_indp <- results(dea_indp)
# load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/dea_unpaired.rds")
# res <- results(dea_paired)
# 
# 
# t(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
# 
# # p-adjusted value
# padj_indp <- res_indp$padj
# padj_paired <- res_paired$padj
# 
# # gene list
# geneList <- rownames(control)

# significantly expressed genes
sig_indp <- geneList[which(padj_indp < 0.05)]
sig_paired <- geneList[which(padj_paired < 0.05)]


















####################################################################################
# Let's look at the  abundance of the microbes across environments
####################################################################################

GP_Pseud <- subset_taxa(GP_at15_all, Genus=="Pseudomonas")

GP_Sphing <- subset_taxa(GP_at15_all, Genus == "Sphingomonas")

GP_Burkholderia <- subset_taxa(GP_at15_all, Family == "Burkholderiaceae")

####################################################################################
# Show phylogeny 
####################################################################################
p.pseud = plot_tree(GP_Pseud, 
                    size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species")+
  scale_color_viridis_d() + ggtitle("Pseudomonas")

#+ coord_polar(theta="y")

p.sping = plot_tree(GP_Sphing, 
                    size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species") +
  scale_color_viridis_d() + ggtitle("Sphingomonas")

p.burk = plot_tree(GP_Burkholderia,
                   size = "abundance", plot.margin = 0.5, ladderize="left", color = "Species")+
  scale_color_viridis_d() + ggtitle("Burkholderiaceae")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/genus_level_phyl_abund.pdf", uuseDingbats = FALSE, fonts = "ArialMT", width = 12, height =12)
plot_grid(p.pseud, p.sping, p.burk, nrow = 1)
dev.off()


####################################################################################
# Compare abundance between environment
####################################################################################

