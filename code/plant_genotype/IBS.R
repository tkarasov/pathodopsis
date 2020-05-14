library(vegan)
library(SNPRelate)
library(gdsfmt)
library(microbiome)
library(phyloseq)
library(ade4)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
#library(genefilter)
library(limma)
library(edgeR)
library(sp)
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
library(dendextend)

#library(VariantAnnotation)

# The goal of this script is to associate microbial similarity with distance with genetic similarity
# basic info on handling vcf file and calculating ibs https://bioconductor.statistik.tu-dortmund.de/packages/3.2/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html

setwd("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/")
hue1_25 = c("#ada77c","#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77","#114477","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122","#DD7788", "darkgoldenrod1", "#771155", "#AA4488", "#CC99BB","#AA7744")

# Load in vcf and climate data
vcf = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolGVCF_gander2.vcf"
vcf_all = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolsGVCF.filtered_snps_final.PASS.bi.vcf"
genot_dir = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/"
output_direc = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin/"
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

# Load genome annotation information
gff = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/TAIR10_GFF3_genes.gff"
genemodels = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/ATH_GO_GOSLIM.txt"


#######################################################################
# VCF Fst analysis
#######################################################################
fst_1_3 <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/poolGVCF_gander_fst_1_3.weir.fst",
                      header = T)
fst_1_3_all<-read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/fst_1_3.weir.fst",
                        header = T)
fst_1_3$new_pos <- c(1:dim(fst_1_3)[1])

acd6 <- fst_1_3[which(fst_1_3$WEIR_AND_COCKERHAM_FST==max(fst_1_3$WEIR_AND_COCKERHAM_FST, na.rm =TRUE)),]$new_pos[1]

#This graph shows the elevated Fst at ACD6
Fst_plot <- ggplot(fst_1_3, aes(x=new_pos, y=WEIR_AND_COCKERHAM_FST, col = CHROM)) + 
  geom_point() + 
  #facet_grid(rows = vars(chr), scales = "free_x", switch = "x") +
  theme_bw() +
  #geom_hline(aes(yintercept = fst99), lty = "dashed") +
  scale_color_viridis_d() +
  xlab("Position") +
  ylab("Fst") +
  geom_vline(aes(xintercept = acd6 ),
             alpha = 0.5) +
  geom_text(aes(x=acd6, label="ACD6\n", y = 0.7), colour="blue", angle=90, text=element_text(size=11)) +
  ylim(c(0,1)) 
#facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x')

fst_1_3_acd6 <- fst_1_3_all %>% filter(CHROM=="chr4") %>% filter(POS>8283409 & POS<8300000)
txdb <- makeTxDbFromGFF(file="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/TAIR10_GFF3_genes.gff", format="gff3")
gene_track <- autoplot(txdb, which=GRanges("Chr4", IRanges(8283409, 8300000)), names.expr = "gene_id")+ theme_bw()

Fst_plot_on_acd6 <- ggplot(fst_1_3_acd6, aes(x=as.numeric(as.character(POS)), y=WEIR_AND_COCKERHAM_FST, col = CHROM)) + 
  geom_point() + 
  theme_bw() +
  scale_color_viridis_d() +
  xlab("Position") +
  ylab("Fst") +
  ylim(c(0,1)) 

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/fst_acd6_gene_model.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
tracks(Fst_plot_on_acd6, gene_track, heights = c(0.7, 0.3))
dev.off()


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/fst_acd6.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
Fst_plot
dev.off()


tped <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_patho.vcf.tped", 
                   header = FALSE, sep = "\t")
tfam <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_patho.vcf.tfam", 
                   header = FALSE, sep = "\t")
colnames(tped)[1:4] = c("CHROM", "p1", "p2", "POS" )
tped_red = tped[, c(1:4, seq(5, dim(tped)[2],2))]

# convert.snp.tped(tped="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tped",
# tfam="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tfam",
# out="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.raw",strand="+")

genot_list <- plant_clim$clim_data


# merge genot_list and tfam on pk and V1 respectively
tfam_info <- dplyr::left_join(tfam, genot_list, by = c("V1" = "Sequence_ID"))
tfam_info <- filter(tfam_info, V1 %in% genot_list$Sequence_ID)
tped_red <- tped_red[,c(1:4,(which(tfam_info$V1 %in% genot_list$Sequence_ID) +4))]
snp_8295146 <- tped_red[which(tped_red$POS==8295844),]
colnames(snp_8295146) = colnames(tped_red)
tfam_info$snp_8295146 <- t(snp_8295146[5:length(snp_8295146)])


#######################################################################
# Map of SNP location
#######################################################################
world <- ne_countries(scale = "medium", returnclass = "sf")

acd6_map <- ggplot(data = world) +
  geom_sf() +
  geom_jitter(data = tfam_info[which(tfam_info$snp_8295146!=0),], 
             aes(x = as.numeric(Long), y = as.numeric(as.character(Lat)),  colour = (snp_8295146)), 
             size = 2,
             alpha =.5
             ) +
  coord_sf(xlim = c(-20, 40), ylim = c(30, 75), expand = FALSE) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(legend.position="bottom") +
  xlab("") + 
  ylab("")


ggplot(data = tfam_info, aes(x = snp_8295146, y = as.numeric(as.character(Lat)))) +
  geom_boxplot() + geom_jitter()

#######################################################################
# Analysis of SNP frequency in 1001 genomes
#######################################################################
tped <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tped", 
                   header = FALSE, sep = "\t")
tfam <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tfam", 
                   header = FALSE, sep = "\t")
colnames(tped)[1:4] = c("CHROM", "p1", "p2", "POS" )
tped_red = tped[, c(1:4, seq(5, dim(tped)[2],2))]

# convert.snp.tped(tped="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tped",
                 # tfam="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.tfam",
                 # out="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_1001.raw",strand="+")

genot_list <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/1001_genomes_list.csv",
                         header = T, sep =",")


# merge genot_list and tfam on pk and V1 respectively
tfam_info <- dplyr::left_join(tfam, genot_list, by = c("V1" = "pk"))
snp_8295146 <- tped_red[which(tped_red$POS==8295844),]
tfam_info$snp_8295146 <- t(snp_8295146[5:length(snp_8295146)])

thousand_map <- ggplot(data = world) +
  geom_sf() +
  geom_jitter(data = tfam_info[which(tfam_info$snp_8295146!=0),], 
              aes(x = as.numeric(longitude), y = as.numeric(as.character(latitude)),  colour = (snp_8295146)), 
              size = 2,
              alpha =.5
  ) +
  coord_sf(xlim = c(-20, 40), ylim = c(30, 75), expand = FALSE) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(legend.position="bottom") +
  xlab("") + 
  ylab("")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/map_acd6.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
plot_grid(acd6_map + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()), 
          thousand_map +theme(axis.title.x=element_blank(),   
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()))

dev.off()

ggplot(data = tfam_info[tfam_info$snp_8295146!="0",], aes(x=longitude, y = latitude, col = snp_8295146)) +
  geom_point()

ggplot(data = tfam_info, aes(x = snp_8295146, y = latitude)) +
  geom_boxplot() +
  geom_jitter()

#######################################################################
# Plot pathodopsis acd6 on hierarchical cluster
#######################################################################
vcf_acd6 <- "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/merged_acd6.vcf"
snpgdsVCF2GDS(vcf_acd6, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.gds", method="biallelic.only")
genofile_acd6 <- snpgdsOpen("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.gds")

snp_table <- read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/merged_tables.txt", 
                        sep = "\t", header = T, row.names = 1)
# #first_match <- match(unique(snp_table$POS), snp_table$POS)
# snp_table <- snp_table[first_match,]
# rownames(snp_table) <- snp_table[,1]
# snp_table <- snp_table %>% select(-1)
# snp_table$poo <- "NA"
# snp_table <- data.frame(t(snp_table))[c(1:75),]
#sapply(snp_table, function(x){x <- gsub("./.",NA,(x))})

convert_snp<-function(x){
  lev <- levels(x)[1]
  xnew <-as.character(x)
  xnew[x!=lev]<-0
  xnew[x==lev]<-1
  return((as.numeric(xnew)))
}

#exclude low
snp_keep <- data.frame(t(snp_table))
snp_fin <- sapply(snp_keep, convert_snp)
good <- which(apply(snp_fin, 1, function(x) sum(is.na(x)))<70)
rownames(snp_fin) <- rownames(snp_keep)
snp_fin <- snp_fin[good,]
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) as.dist((1-cor(t(x), use="all.obs")/2))
d <- dist(snp_fin)#dist(snp_fin, method = "euclidean") #, metric = c("gower"))
hc <- hclustfunc(d) #, method = "complete")
dend <- as.dendrogram(hc)
colLab <- snp_fin[labels(dend),"X8297824"]
colLab[colLab==1]<-"RED"
colLab[colLab==0]<-"GREEN"
labels_colors(dend) <- colLab 
labels(dend) <- gsub(".GT","", labels(dend))
clust_dend<-plant_clim$clim_data[match( labels(dend),plant_clim$clim_data$Plant_ID),]$cluster
labels_colors(dend)<-clust_dend*3
plot(dend)



snp_table1 <- sapply(snp_table, function(x){x <- gsub("./.",NA,x)})
rownames(snp_table1) <- rownames(snp_table)
pca <- snpgdsPCA(genofile_acd6, num.thread = 4, maf=0.05)
#######################################################################
# Calculate kinship and pca
#######################################################################
#snpgdsVCF2GDS(vcf, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/pathodopsis.gds", method="biallelic.only")
#snpgdsVCF2GDS(vcf_all, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/pathodopsis_all.gds", method="biallelic.only")
clim = plant_clim$clim_data[, c("Plant_ID", "cluster")]
rownames(clim) = clim$Plant_ID
genofile <- snpgdsOpen("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/pathodopsis.gds")
genofile_all <- snpgdsOpen("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/pathodopsis_all.gds")
samp.id <- read.gdsn(index.gdsn(genofile_all, "sample.id"))
snp.id_all <- read.gdsn(index.gdsn(genofile_all, "snp.id"))
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.chromosome <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.position <- read.gdsn(index.gdsn(genofile, "snp.position"))
clim = clim[samp.id,]
clim = clim[!is.na(clim[,1]),]

ggplot(pca,aes(x=eigenvect[,1], y=eigenvect[,2])) + geom_point()
pca <- snpgdsPCA(genofile, num.thread = 4, maf=0.05)#, snp.id=snpset.id, num.thread=2)
pca$color <- grepl("PA", pca$sample.id)
fst <- snpgdsFst(genofile, sample.id = clim$Plant_ID, 
                 with.id = TRUE,
                 snp.id = snp.id, population = as.factor(clim$cluster), maf = 0.10)
fst_all <- snpgdsFst(genofile_all, sample.id = clim$Plant_ID, 
                 with.id = TRUE,
                 snp.id = snp.id, population = as.factor(clim$cluster), maf = 0.10)
ibs = snpgdsIBS(genofile, num.thread = 4)
kinship = ibs$ibs
snpgdsClose(genofile)
snpgdsClose(genofile_all)

#######################################################################
# Plot Fst along genome
#######################################################################
fst_mat = data.frame(chr = snp.chromosome[fst$snp.id], position = snp.position[fst$snp.id], fst = fst$FstSNP)
fst_quantile = quantile(fst_all$FstSNP, c(0,.999,1), na.rm =T)
fst99 <- as.numeric(fst_quantile[2])
fst_mat$new_pos <- c(1:dim(fst_mat)[1])
acd6 = fst_mat[which(fst_mat$fst==max(fst_mat$fst, na.rm=TRUE)),]$new_pos[1]

tsv = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.tsv", sep="\t")
hm=tsv[,c(28,1:20)]
mine=fst_mat[which(fst_mat$fst == max(fst_mat$fst, na.rm=T)),]


#This graph shows the elevated Fst at ACD6
Fst_plot <- ggplot(fst_mat, aes(x=new_pos, y=fst, col = chr)) + 
  geom_point() + 
  #facet_grid(rows = vars(chr), scales = "free_x", switch = "x") +
  theme_bw() +
  geom_hline(aes(yintercept = fst99), lty = "dashed") +
  scale_color_viridis_d() +
  xlab("Position") +
  ylab("Fst") +
  geom_vline(aes(xintercept = acd6 ),
             alpha = 0.1) +
  geom_text(aes(x=acd6, label="ACD6\n", y = 0.7), colour="blue", angle=90, text=element_text(size=11)) +
  ylim(c(0,1)) 
  #facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x')


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/fst_acd6.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
Fst_plot
dev.off()

#And just within ACD6
Fst_plot <- ggplot(fst_mat, aes(x=new_pos, y=fst, col = chr)) + 
  geom_point() + 
  #facet_grid(rows = vars(chr), scales = "free_x", switch = "x") +
  theme_bw() +
  geom_hline(aes(yintercept = fst99), lty = "dashed") +
  scale_color_viridis_d() +
  xlab("Position") +
  ylab("Fst") +
  geom_vline(aes(xintercept = acd6 ),
             alpha = 0.1) +
  geom_text(aes(x=acd6, label="ACD6\n", y = 0.7), colour="blue", angle=90, text=element_text(size=11)) +
  ylim(c(0,1)) 
#facet_grid(~chr, scales = 'free_x', space = 'free_x', switch = 'x')





#######################################################################
# Annotate vcf
#######################################################################
TAIR10 <- "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/TAIR10_GFF3_genes.gff"
TAIR10_gff <- import.gff(TAIR10) 
myGranges<-as(TAIR10_gff, "GRanges") 
TAIR10_genes <- myGranges[which(myGranges$type=="gene"),]
anno_track <- AnnotationTrack(myGranges, start = 1, 
                              end =1000, feature = "1",
                              stacking = "squish")
fst_df <- data.frame(chr=paste("Chr", fst_mat$chr, sep=""),
                   start=fst_mat$position,
                   end=fst_mat$position,
                   fst=fst_mat$fst,
                   #id=c('one', 'two', 'three'),
                   strand=c('+'))

fst_g <- with(fst_df, myGranges, IRanges(start, end), strand, id=id)

fst_track <- DataTrack(data = fst_mat$fst, start = fst_mat$position,
                    end = fst_mat$position,
                    chromosome = paste("Chr",fst_mat$chr, sep=""),
                    genome = genome(myGranges),
                    strand = c("+"))
  #fst_g, strand = c("+"))

annot_trac <- GeneRegionTrack(ranges = TAIR10_genes,
                             chromosome = "Chr1", start = 1000, end =1000000)

variants <-read.vcf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.recode.vcf")
variants2 <- VCFloci("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.recode.vcf")
variants3<-"/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/mapping/acd6_parsed.vcf"




snpgdsVCF2GDS(variants3, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.gds", method="biallelic.only")
genofile_acd6 <- snpgdsOpen("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6.gds")
samp.id <- read.gdsn(index.gdsn(genofile_acd6, "sample.id"))
snp.id <- read.gdsn(index.gdsn(genofile_acd6, "sample.id"))
fst <- snpgdsFst(genofile_acd6, sample.id = clim$Plant_ID, 
                 with.id = TRUE,
                 snp.id = snp.id, population = as.factor(clim$cluster), maf = 0.10)
colnames(variants) <- variants2$POS
mine <- variants[,which(colnames(variants)=="8295146")]
samps = as.character(plant_clim$clim_data$Sequence_ID)
plant_clim$clim_data$snp8295146 <- mine[samps,1]$`8295146`
#######################################################################
# Look at Variogram
#######################################################################

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")

plant_val = which(OTU_clim$clim_data$Host_Species=="Ath")

plant_clim <- list(otu_table = OTU_clim$otu_table[plant_val,], 
                   clim_data = OTU_clim$clim_data[plant_val,], 
                   tax_table = OTU_clim$tax_table, 
                   phy_tree = OTU_clim$phy_tree, 
                   refseq = OTU_clim$refseq)

plant_otu = data.frame(plant_clim$otu_table)
plant_otu$Latitude = as.numeric(as.character(plant_clim$clim_data$Lat))
plant_otu$Longitude = as.numeric(as.character(plant_clim$clim_data$Long))
plant_otu <- plant_otu %>% filter(is.na(Latitude)==FALSE)

set.seed(4)
plant_otu$Longitude = plant_otu$Latitude + runif(length(plant_otu$Latitude), 0, 0.0001)

# no trend:
my.transform <- "+proj=longlat +datum=WGS84 +no_defs"
sp::coordinates(plant_otu) = ~Longitude+Latitude
proj4string(plant_otu) = CRS(my.transform)
#proj4string(plant_otu) = CRS("+init=epsg:28992")

opti_var <- names(plant_otu)

# ################
# # LONG distance
# ################
# for(i in 1:20){
#   eq <- as.formula(paste(paste("log10(",opti_var[i],sep=""), "+0.01)~1", sep=""))
#   v = variogram(eq, plant_otu, cutoff = 1000, width = 20)
#   v.plot = ggplot(v, aes(x = dist, y = gamma)) + 
#     geom_point(alpha = 0.3, cex = 2) +
#     xlab("Distance (km)") +
#     ylab("Semivariance") +
#     annotate("Text", y= max(v$gamma), x =max(v$dist),label=opti_var[i],hjust=1) +
#     theme_bw()
#   assign(paste("plot", i, sep=""), v.plot)
#   myplots[[i]] = v.plot
# 
# }
# 
# 
# long <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, plot13, plot14, plot15, plot16, plot17, plot18, plot19, plot20)

my.sph<- vector('list', length(opti_var))
my.lin<- vector('list', length(opti_var))
################
# Look at variograms over 1000km
################
for(i in 1:dim(plant_clim$otu_table)[2]){
  my.label <- paste(plant_clim$tax_table[i], opti_var[i], sep = ",")[6]
  eq <- as.formula(paste(paste("log10(",opti_var[i],sep=""), "+0.01)~1", sep=""))
  v = variogram(eq, plant_otu, cutoff = 1000, width = 8)
  v.fit.lin = fit.variogram(v, vgm("Lin"))
  vgl = variogramLine(v.fit.lin, maxdist = 1000)
  v.plot = ggplot(v, aes(x = dist, y = gamma)) + 
    geom_point(alpha = 0.3, cex = 2) +
    geom_line(data = vgl, col = "RED") + 
    xlab("") +
    ylab("") +
    #ylim(c(0,3)) +
    labs(title = my.label, size = 4) +
    #annotate("Text", y = max(v$gamma), x =max(v$dist), label = my.label, hjust=1, size = 6 ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    assign(paste("plot", i, sep=""), v.plot)
  

 # v.fit.exp = fit.variogram(v, vgm(1, "exp"))
#  my.sph[[i]] = v.fit.sph
  my.lin[[i]] = v.fit.lin
 # my.exp[[i]] = v.fit.exp
  
}
names(my.lin) = opti_var

short <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, nrow = 4, label_x ="Distance (km)", label_y = "Semivariance")
#, plot11, plot12, plot13, plot14, plot15, plot16, plot17, plot18, plot19, plot20)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/variogram_100km.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
short
dev.off()

###################################################
# Calculate Decay rate
###################################################
my.nugget <- lapply(my.lin, `[[`, 2)
my.range <- lapply(lapply(my.lin, `[[`, 3), '[[', 2)
my.slope <- as.numeric(lapply(lapply(my.lin, `[[`, 2), '[[', 2))
names(my.slope) <- names(my.range)

my.range.1 <- as.numeric(as.character(my.range))
info_range_orig <- data.frame(range = my.range.1, genus = plant_clim$tax_table, slope = my.slope)
info_range_orig$max_abund <- apply(plant_clim$otu_table, 2, max, na.rm = TRUE)/1000
info_range <- info_range_orig %>% filter(slope!=0)
info_range_orig$max_abund <- apply(plant_clim$otu_table, 2, max, na.rm = TRUE)/1000
info_range_orig$mean_abund <- apply(plant_clim$otu_table, 2, mean, na.rm = TRUE)/1000

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/variogram_across_families.pdf", useDingbats = FALSE, font = "ArialMT", width  = 10, height = 4)
ggplot(info_range, aes(x = fct_reorder(genus.Family, range, .fun = median, .desc =TRUE), y = range)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point(alpha = 0.2) +
  ylab("Distance (km)") +
  xlab("")
#  facet_wrap(~genus.Family)
dev.off()

mean(info_range$range)
median(info_range$range)
sd(info_range$range)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/hist_autocorrelated_range.pdf", useDingbats = FALSE, font = "ArialMT", width  = 3.5, height = 3.5)
ggplot(info_range, aes(x = range)) + 
  geom_histogram(alpha = 0.3, color = "Grey20") +
  scale_x_continuous(trans = 'log10', breaks = c(1,10,100, 1000,4000)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Distance (km)") +
  ylab("Number of Phylotypes")

#  xlim(c(0, max(info_range$range)))
dev.off()


###################################################
# Relationship between OTUS and Latitude and Longitude
###################################################
#correlate every OTU with Latitude
y = data.frame(plant_clim$otu_table)
plant_clim$clim_data$Lat <- as.numeric(as.character(plant_clim$clim_data$Lat))
plant_clim$clim_data$Long <- as.numeric(as.character(plant_clim$clim_data$Long))
hm = apply(y, 2, function(col1) cor.test(col1, plant_clim$clim_data$Lat, na.rm =TRUE)[c("p.value", "estimate")])
hm.long = apply(y, 2, function(col1) cor.test(col1, plant_clim$clim_data$Long, na.rm =TRUE)[c("p.value", "estimate")])
cor.hm <- matrix(unlist(hm), ncol = 2, byrow = TRUE)
cor.hm.long <- matrix(unlist(hm.long), ncol = 2, byrow = TRUE)
colnames(cor.hm) <- c("pval.lat", "cor.lat")
colnames(cor.hm.long) <- c("pval.long", "cor.long")
info_range_orig<- cbind(info_range_orig, cor.hm)
info_range_orig<- cbind(info_range_orig, cor.hm.long)
info_range_orig$pval.adj <- p.adjust(info_range_orig$pval.lat, method = "BH")
info_range_orig$pval.adj.long <- p.adjust(info_range_orig$pval.long, method = "BH")

#number correlated with Latitude
length(which(info_range_orig$pval.adj<0.01))/length(info_range_orig$pval.adj)

#number correlated with Longitude
length(which(info_range_orig$pval.adj.long<0.01))/length(info_range_orig$pval.adj)

###################################################
# OTU5 picture
###################################################
hist(plant_clim$otu_table[,"seq_10"]/1000)

otu5_hist <- ggplot(data.frame(plant_clim$otu_table), aes(x = (seq_10 + 1)/10)) + 
  geom_histogram(alpha = 0.3, color = "Grey20") +
  scale_x_continuous(trans = 'log10', breaks = c(0,0.1,1, 5,10,20,50,70)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("OTU5 pathogen RA (%)") +
  ylab("Number of Plants")

otu5_spatial <- plant_clim$clim_data$otu5 <- as.numeric(plant_clim$otu_table[,"seq_10"][,1])
otu5_lat <- ggplot(data = plant_clim$clim_data, aes(x = Lat, y = otu5)) +
  geom_point() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 




i = 7
my.label <- paste(plant_clim$tax_table[i], opti_var[i], sep = ",")[6]
eq <- as.formula(paste(paste("log10(",opti_var[i],sep=""), "+0.01)/10~1", sep=""))
v = variogram(eq, plant_otu, cutoff = 1000, width = 8)
v.fit.lin = fit.variogram(v, vgm("Sph"))
vgl = variogramLine(v.fit.lin, maxdist = 1000)
v.plot = ggplot(v, aes(x = dist, y = gamma)) + 
  geom_point(alpha = 0.3, cex = 2) +
  geom_line(data = vgl, col = "RED") + 
  ylab(expression(paste("Semivariance (",paste(log[10], "(RA %))"), sep = ""))) +
  xlab("Distance (km)") +
  #ylim(c(0,3)) +
  #labs(title = my.label, size = 4) +
  #annotate("Text", y = max(v$gamma), x =max(v$dist), label = my.label, hjust=1, size = 6 ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

v.plot



load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_OTU5.rds")
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU5_all_subfigures.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
otu5_all <- grid.arrange(otu5_grid, plot_grid(otu5_hist, v.plot, ncol = 1, align = "hv"), ncol =2)
dev.off()
