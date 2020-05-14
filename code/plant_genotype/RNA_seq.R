
library(dplyr)
library(gstat)
library(ggplot2)
library(viridis)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(pegas)
library(sf)
library(rnaturalearth)
library(cowplot)
library(org.At.tair.db)


#############
# Load the RNAseq data
rna <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/RNAseq/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv", 
                  sep = "\t", header = TRUE, row.names = 1)
genot_list <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/1001_genomes_list.csv",
                         header = T, sep =",")
immunity_list <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/defense_genes_gander.txt", sep = "\t", header = T)
# pr1 <- rna[which(rownames(rna)=="AT1G64280"),]
# pr1 <- data.frame(t(pr1))
# pr1$genot = as.integer(as.character(gsub("X", "", rownames(pr1))))
# genot_list <-dplyr::left_join(pr1, genot_list, by = c("genot" = "pk"))


# all RNA
all_rna <- data.frame(t(rna))
all_rna$genot = as.integer(as.character(gsub("X", "", rownames(all_rna))))
genot_list <-dplyr::left_join(all_rna, genot_list, by = c("genot" = "pk"))

thousand_map <- ggplot(data = world) +
  geom_sf() +
  geom_jitter(data = genot_list, 
              aes(x = as.numeric(longitude), y = as.numeric(as.character(latitude)),  colour = (AT1G64280)), 
              size = 2,
              alpha =.5
  ) +
  coord_sf(xlim = c(-20, 40), ylim = c(30, 75), expand = FALSE) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(legend.position="bottom") +
  xlab("") + 
  ylab("")

genes <- genot_list[,grep("AT", colnames(genot_list))]
immunity <- genes #genes[,immunity_list[,1]]
corr_all <- data.frame(t(sapply(immunity, function(x)cor.test(x, genot_list$latitude, na.rm = T))))
hist(data.frame(corr_all$estimate))
corr_ <- as.numeric(data.frame(corr_all$estimate)[1,])
corr_mat <- data.frame(rho=as.numeric(data.frame(corr_all$estimate)[1,]), pval=as.numeric(data.frame(corr_all$p.value)[1,]))
corr_mat$fdr <- p.adjust(corr_mat$pval, method = "BH")

