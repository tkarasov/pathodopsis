library(phyloseq)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(dplyr)

####################
# Load phyloseq object & pathogen identities
####################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds')

samp_data <- plant_clim$clim_data
rownames(samp_data) <- samp_data$Sequence_ID
GP_red <- phyloseq::prune_samples(rownames(samp_data), GP1000)
sample_data(GP_red) <- sample_data(plant_otu)

plant_otu <- phyloseq(sample_data(samp_data), 
                      otu_table = plant_clim$otu_table, 
                      phy_tree = plant_clim$phy_tree, 
                      refseq = plant_clim$refseq,
                      tax_table = plant_clim$tax_table)


####################
# Basic plot of abundant taxa
####################
physeq3 = transform_sample_counts(GP_red, function(x) x / sum(x) )
glom <- tax_glom(physeq3, taxrank = 'Phylum')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$class <- as.character(data_glom$Phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$class[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$class))
Count

spatial_plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=class)) + 
  facet_grid(~cluster, scales = "free")
spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "black")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))


plot_bar(plant_otu,"Family", fill = "Genus", facet_grid=~cluster)



# #######################################################################
# For differential analysis use edgeR
# #######################################################################
# https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html
#https://joey711.github.io/phyloseq-extensions/edgeR.html

phyloseq_to_edgeR <- function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}


#Run new function on plant_otu
dge = phyloseq_to_edgeR(plant_otu, group="cluster")
# Perform binary test
extract_results <- function(dge, phylo){
  et = exactTest(dge)
  # Extract values from test results
  tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
  res = tt@.Data[[1]]
  alpha = 0.01
  sigtab = res[(res$FDR < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo)[rownames(sigtab), ], "matrix"))
  dim(sigtab)
  return(sigtab)
}

#########
extract_ASV <- extract_results(dge, plant_otu)
sigtab <- extract_ASV

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))

# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))

# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))

# Plot
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/compare_taxa_clusters.pdf",  width = 7, height = 3.5,
    useDingbats = FALSE, font = "ArialMT")

ggplot(sigtabgen, aes(x = Genus, y = -1*logFC, color = Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  scale_color_viridis_d() +
  ylab("log2FC(Northern Cluster/Southern Cluster)")

dev.off()


##################
## Get a random representation of abundance


plant_otu <- phyloseq(sample_data(samp_data), 
                      otu_table = plant_clim$otu_table, 
                      phy_tree = plant_clim$phy_tree, 
                      refseq = plant_clim$refseq,
                      tax_table = plant_clim$tax_table)

# Choose some North and South
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP50.rds")

set.seed(4)
north = which(sample_data(plant_otu)$cluster==2)
south = which(sample_data(plant_otu)$cluster==1)
n20 <- sample_data(plant_otu)$Plant_ID[sample(north, 25)]
s20 <- sample_data(plant_otu)$Plant_ID[sample(south, 25)]
n_s <- as.factor(c(as.character(n20), as.character(s20)))
sub_sample <- subset_samples(GP1000, samples.out %in% n_s)
sample_data(sub_sample)$cluster <- sample_data(subset_samples(plant_otu, Sequence_ID  %in% n_s))$cluster


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/temp.pdf",
    useDingbats = FALSE, font = "ArialMT")
plot_bar(sub_sample, fill="Phylum" ) +
  facet_grid(~cluster, scale="free_x", drop=TRUE) + 
  scale_color_manual(c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744", "#AA4455"))


dev.off()
