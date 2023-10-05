# Horton analysis
# This script was originally written to run on TK's laptop. On 10/2023 she moved it to github for posterity.
# The goal of this script is to look at the SNPs that Matt found to be significant in structuring the microbiome and to look whether these show global fst 
library(dplyr)
library(ggplot2)

sig_snps <- read.table("~/Dropbox/pathodopsis/horton_analysis/horton_sig_snps_all.txt", header = TRUE, sep = "\t")
global_fst <- read.table("~/Dropbox/pathodopsis/horton_analysis/fst/global.fst.all.txt", header = TRUE, sep = ",")
sig_regions <- sig_snps[,c("chr", "pos", "name")]

there <- merge(sig_regions, global_fst, by= c("chr", "pos"))
sig_regions$chr <- paste("chr", sig_regions$chr, sep="")

# merge regions


# compare fst values between the groups
ggplot(data = there, aes(x = pos, y = fst, colour = chr)) +
  geom_point()


# ask question of whether snps for the loci of interest tend to have higher global fsts than the genome overall
df1 <- data.frame(global_fst$fst )
df1$condition <- "Global"
df_2 <- data.frame(there$fst)
df_2$condition <- "Significant"
colnames(df1) <- c("fst", "condition")
colnames(df_2) <- c("fst", "condition")
df3 <- rbind(df1,df_2)

pdf("/Users/talia/Dropbox/pathodopsis/horton_analysis/horton_fst.pdf", family = "ArialMT", useDingbats = F)
ggplot(data = df3, aes(x= condition, y = fst)) + geom_boxplot(outlier.shape = NA) + theme_bw()
dev.off()

global <- df3 %>% filter(condition == "Global")
significant <- df3 %>% filter(condition == "Significant")

#######################################################################
# Plot Fst along genome
#######################################################################

# These were the fst values calculated by me a while back
# fst_1_3 <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/fst_1_2.weir.fst",
#                      header = T)
fst_1_2 <- read.table("/Users/talia/Dropbox/pathodopsis/horton_analysis/fst_1_2.weir.fst", header = T)

fst_1_2 <- unique(fst_1_2)

fst_1_2$new_pos <- c(1:dim(fst_1_2)[1])
fst_1_2$what <- "all"
sig_regions$what <- "sig"
colnames(fst_1_2)[1:2] <- c("chr", "pos")
my_fst <- merge(fst_1_2, sig_regions, by = c("chr", "pos"))

there2 <- there
# there2$chr <- paste("chr", there$chr, sep="")
tog_fst <- merge(there2, my_fst)
cor.test(tog_fst[,4], tog_fst[,5])
plot(tog_fst[,4], tog_fst[,5])
# there is a very poor relationship between the global fst calculated by Matt and the Fst between populations 1 and 2that I calculated


#This graph shows the elevated Fst at ACD6
Fst_plot <- ggplot(fst_1_2, aes(x=new_pos, y=WEIR_AND_COCKERHAM_FST, col = chr)) + 
  geom_point() + 
  #facet_grid(rows = vars(chr), scales = "free_x", switch = "x") +
  theme_classic() +
  #geom_hline(aes(yintercept = fst99), lty = "dashed") +
  scale_color_viridis_d() +
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  ylab("Fst") +
  geom_vline(aes(xintercept = acd6 ),
             alpha = 0.5) +
  geom_text(aes(x=acd6, label="ACD6\n", y = 0.7), colour="blue", angle=90, text=element_text(size=11)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ylim(c(0,1)) 


#######################################################################
# Added on 10/03/2023 Find the SNPs studied by Horton et al. in the 1001 vcf files and calculate frequency of these SNPs with latitude.
#######################################################################
#these were the SNPs in A. thaliana that were significant in the Horton/Bergelson study
horton <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/horton_sig_snps_all.txt", header = T)
horton_chr_pos <- horton[, c("chr", "pos")]
