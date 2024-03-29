#!/usr/bin/env Rscript

# Map the distributions of a list of pathogens
# tutorial on sf package: https://geocompr.robinlovelace.net/spatial-operations.html#spatial-ras
library(sf)
library(raster)
library(spData)
library(spDataLarge)
library(gstat)
#library(tmap)
library(ggplot2)
library(wesanderson)
library(spdep)
#library(proj4) #proj4 has major installation problems
library(smoothr)
library(tidyverse)
library(gstat)
library(automap)
data(world)
library(cowplot)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(showtext)
library(extrafont)
library(intrval)
Arial = font_add("Arial", "/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/pathodopsis/misc_tkarasov/fonts/Arial.ttf")
library(automap)
library(devtools)
library(dplyr)
library(phyloseq)
library(reshape2)

theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_rect(fill = 'light grey', colour = NA),
                       panel.border = element_rect(colour = "dark grey", fill=NA, size=.5),
                       #plot.background = element_rect(fill="light grey",
                       #size=1,linetype="solid",color="black"),
                       #axis.line = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.title = element_text(size=22),
                       aspect.ratio=1,
                       legend.justification=c(0,1), 
                       legend.position=c(0.75, 0.95),
                       legend.background = element_blank()))


my.transform <- "+proj=longlat +datum=WGS84 +no_defs"
#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#Laea for europe: https://spatialreference.org/ref/epsg/etrs89-etrs-laea/
laea <- "+init=epsg:3857" 

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


####################
# Load phyloseq object & pathogen identities
####################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds')
my_bac <- read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_bacteria.csv",
  header = TRUE, row.names = 1)
my_oom <- read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_oomycete.csv",
  header = TRUE, row.names = 1)
my_fungi = read.csv(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/meta_family_corrected_per_plant_v2_fungi.csv",
  header = TRUE, row.names = 1)

my_phylo <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/seqtab_final.rds")
my_phylo_tax <- readRDS('/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/tax_final.rds')

####################
# Build Lat/Long Dataframe
####################
clim <- plant_clim$clim_data
rownames(clim) <- clim$Sequence_ID
plant_phylo <- phyloseq(sample_data(clim), refseq = plant_clim$refseq, tax_table = plant_clim$tax_table, phy_tree = plant_clim$phy_tree, otu_table = plant_clim$otu_table )
sample_data(GP1000)$Lat <- as.numeric(as.character(sample_data(GP1000)$Lat))
sample_data(GP1000)$Long <- as.numeric(as.character(sample_data(GP1000)$Long))
#GP1000_dim <- subset_samples(GP1000, sample_data(GP1000)$samples.out %in% colnames(my_bac))
GP1000_dim <- subset_samples(plant_phylo, sample_data(plant_phylo)$Sequence_ID %in% colnames(my_bac))
reorder <- as.character(sample_data(GP1000_dim)$Sequence_ID)
my_bac <- my_bac[, reorder]
my_oom <- my_oom[, reorder]
my_fungi <- my_fungi[, reorder]
all_microb <- rbind(my_bac, my_oom, my_fungi)

load_tot <- data.frame(plant_ID = colnames(my_bac), 
                       oom_load = colSums(my_oom), 
                       fung_load = colSums(my_fungi), 
                       bac_load = colSums(my_bac, na.rm = TRUE),
                       hpa_load = t(my_oom["Peronosporaceae",]),
                       albug_load = t(my_oom["Albuginaceae",]))

sample_data(GP1000_dim)$oom_load = load_tot$oom_load
sample_data(GP1000_dim)$fung_load = load_tot$fung_load
sample_data(GP1000_dim)$bac_load = load_tot$bac_load
sample_data(GP1000_dim)$hpa_load = load_tot$Peronosporaceae
sample_data(GP1000_dim)$albug_load = load_tot$Albuginaceae
sample_data(GP1000_dim)$albug_load = load_tot$Albuginaceae
sample_data(GP1000_dim)$otu5 <- otu_table(GP1000_dim)[,"seq_10"]/1000

tot_frame <- data.frame(sample_data(GP1000_dim))
tot_frame$Lat <- as.numeric(as.character(tot_frame$Lat))


####################
# Correspondence between 16S and metagenome
####################

# aggregate ASV table at family
otu_table(GP1000_dim) <- otu_table(GP1000_dim)/1000
GP_fam <- tax_glom(GP1000_dim, "Family")
shared_fam <- c(tax_table(GP_fam)[,5][which(tax_table(GP_fam)[,5] %in%  rownames(my_bac))])
keep <-rownames(tax_table(GP_fam)[which(c(tax_table(GP_fam)[,5]) %in% shared_fam),])
GP_shared <- prune_taxa(keep, GP_fam)

my_bac_red <- t(my_bac[c(tax_table(GP_shared)[,5]),])
my_bac_red <- my_bac_red/(rowSums(my_bac_red, na.rm = TRUE) + 0.0000001)
GP_tab <- otu_table(GP_shared)
colnames(GP_tab) <- colnames(my_bac_red)
GP_tab <- GP_tab[,order(colSums(GP_tab), decreasing=TRUE)]
my_bac_red <- my_bac_red[,colnames(GP_tab)]
GP_10 <- melt(GP_tab[,c(1:10)], id=c("ID", "Family", "16S"))
my_bac_10 <- melt(my_bac_red[,c(1:10)])
#GP_10$Var3 <- "16S"
#my_bac_10$Var3 <- "Metagenome"
tog <- GP_10
tog$metagenome <- my_bac_10$value
colnames(tog) <- c("PlantID", "Family", "Amplicon", "Metagenome")

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/metagenome_16S.pdf", useDingbats = FALSE, fonts = "ArialMT", width = 3.5, height = 2)
ggplot(data = tog, aes(x = Amplicon*100, y = Metagenome*100)) +
  geom_point(aes(col = Family)) +
  scale_color_brewer(palette = "Paired", guide = guide_legend(override.aes = list(size = 3,
                                                                                  alpha = 1) ) ) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "Grey", linetype = "dashed") +
  xlab("16S Amplicon RA (%)") +
  ylab("Metagenomic RA (%)")

dev.off()

####################
# Build maps of abundance of pathogens
####################
####################################################################################
# Read in data and project data and europe to laea centered on europe
####################################################################################


#First add some wiggle to the OTU5 points then reproject the points
tot_frame$Long = tot_frame$Long + runif(n=length(tot_frame$Long), min=0, max=0.000001)
tot_frame$Lat = tot_frame$Lat + runif(n=length(tot_frame$Lat), min=0, max=0.000001)
tot_frame <- tot_frame %>% filter(is.na(Lat)==FALSE) %>% filter(is.na(Long)==FALSE)
#Now make sure everything has the same projection
coordinates(tot_frame) = c("Long", "Lat")
proj4string(tot_frame) = CRS(my.transform)
my.sf = st_as_sf(x = tot_frame, coords = c("Long", "Lat"), crs = my.transform)
my.sp = as_Spatial(my.sf)
my.sp = spTransform(my.sp, CRS(my.transform))
my.df = data.frame(my.sp)

# subset world map. Still need to remove canary islands and northern part of Norway
europe = world[world$continent == "Europe" | world$continent == "Asia" ,]
europe = as_Spatial(europe, )
crs(europe) = my.transform
europe.sf = st_as_sf(europe)

####################################################################################
# Smooth over the region with Kriging. 
####################################################################################
# Kriging functions
# take only the points falling in polygons. This is the step where we limit the infomration
krig_it <- function(kriging_result){
  Krig = kriging_result$krige_output[!is.na(over(kriging_result$krige_output,as(europe.laea,"SpatialPolygons"))),]  
  Krig.fin = spTransform(Krig, my.transform)
  Krig_df = as.data.frame(Krig.fin)
  names(Krig_df) = c(  "APPT_pred","APPT_var","APPT_stdev", "longitude","latitude")
  return(Krig_df)
}

#I'm going to use the function krige from gstat
#I need to create a grid of points I would like to predict the values for.
lon <- seq(extent(my.sp)[1] - 5, extent(my.sp)[2] + 5, length.out = 500)
lat <- seq(extent(my.sp)[3] -5, extent(my.sp)[4] + 5, length.out = 500)
grd <- expand.grid(lon = lon, lat = lat)
grd_sf  <-  st_as_sf(grd, coords = c("lon", "lat"), crs = my.transform, agr = "constant")
grd_sp <- as_Spatial(grd_sf)
crs(grd_sp) = my.transform

#automap seems a little faster
sp.laea <- spTransform(my.sp, CRS(laea))
europe.laea <-spTransform(europe, CRS(laea))
grd.laea <- spTransform(grd_sp, CRS(laea))

#Run Kriging interpolation for features of interest
kriging_result_hpa = autoKrige(formula = log10(hpa_load+0.0001)~1, input_data = sp.laea, new_data = grd.laea)
kriging_result_albug = autoKrige(formula = log10(albug_load+0.0001)~1, input_data = sp.laea, new_data = grd.laea)
kriging_result_tot = autoKrige(formula = log10(bac_load+0.0001)~1, input_data = sp.laea, new_data = grd.laea)
kriging_result_otu5 =  autoKrige(formula = otu5~1, input_data = sp.laea, new_data = grd.laea)

#Now pull out the relevant dataframe
my.df = as.data.frame(my.sp)
Krig_df_tot <- krig_it(kriging_result_tot)
Krig_df_hpa <- krig_it(kriging_result_hpa)
Krig_df_albug <- krig_it(kriging_result_albug)
Krig_df_otu5 <- krig_it(kriging_result_otu5)


####################################################################################
# Plot world map with data points
####################################################################################
xmin_plot = extent(my.sp)[1] -5
xmax_plot = extent(my.sp)[2] +5
ymin_plot = extent(my.sp)[3] -5
ymax_plot = extent(my.sp)[4] +5

pal <- wes_palette("Zissou1", 100, type = "continuous")

base_europe_df <-
  ggplot() +
  geom_raster(data = Krig_df_tot, 
              aes(x = Krig_df$longitude, 
                  y = Krig_df$latitude, fill = NULL
              ), 
              inherit.aes = T)

#lets just make a raster of the Krig_df
tot_point <-
  base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = log10((bac_load)/10)), 
              width = .5, cex = 2) +
  scale_colour_gradientn(colours = pal,
                         values = c(0, 0.1, 0.4, 1), name = "log10(RA(%))") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_plot, ymax_plot)) +
  scale_x_continuous(expand = c(0,0), limits = c(xmin_plot, xmax_plot)) +
  theme_opts + 
  #theme_pubs +
  xlab(label = NULL) +
  ylab(label = NULL)


####################################################################################
# Plot Kriging Interpolation
####################################################################################
plot_krig <- function(Krig_df, lab = "RA (%)", breaks = seq(0,1,.05)){
  ggplot(data = my.df, aes(x=coords.x2, y=coords.x1)) +
  geom_raster(data = Krig_df, 
              aes(x = Krig_df$longitude, 
                  y = Krig_df$latitude, 
                  fill = Krig_df$APPT_pred), 
              inherit.aes = T) + 
  scale_y_continuous(expand = c(0,0), 
                     limits = c(ymin_plot, ymax_plot)) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(xmin_plot, xmax_plot)) +
  scale_fill_distiller(palette = "Spectral", 
                       labs(size= lab), 
                       breaks = breaks) +
  theme_opts +
  #theme_pubs  +
  xlab(label = NULL) +
  ylab(label = NULL) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

plot_tot <- plot_krig(Krig_df_tot, "log10(Bacterial Load)", 
                      breaks = seq(-4,2, .4))
plot_albug <- plot_krig(Krig_df_albug, "log10(Load Albugo)",
                        breaks = seq(-4,2, .4))
plot_hpa <- plot_krig(Krig_df_hpa, "Load HpA")
plot_otu5 <- plot_krig(Krig_df_otu5, "RA(%) OTU5" )

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_all_path.pdf", family = "ArialMT", useDingbats = F)
all_grid <-plot_grid(plot_tot, plot_albug, plot_hpa, plot_otu5, nrow = 2)
all_grid
dev.off()

save(otu5_grid, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_OTU5.rds" )

####################################################################################
# What environmental factor relates to disease index?
####################################################################################
lm_rec<- data.frame(otu = colnames(plant_clim$otu_table), pval = c(NA)*575)

for(i in 1:dim(plant_clim$otu_table)[2]){
  df = data.frame(OTU = plant_clim$otu_table[,i], Lat =  plant_clim$clim_data$Disease_RL ) 
  colnames(df) = c("OTU", "Lat")
  model <- lm( OTU~Lat, data = df)
  lm_rec[i,2] = summary(model)$coefficients[2,4]
  lm_rec$freq[i] = sum(plant_clim$otu_table[,i])/(1000*dim(plant_clim$otu_table)[1])
}

lm_rec$FDR <- p.adjust(lm_rec$pval, method = "BH")

lm_rec$sig_0.01 <- lm_rec$FDR<0.01

lm_rec <- cbind(lm_rec, plant_clim$tax_table )

lm_fam <- lm_rec %>% filter(Family %in% names(which(table(lm_rec$Family)>=5)))

# Most significant associations with Lat come from Sphingomonadaceae

pval_plot <- ggplot(aes(x = Family, y = -log10(FDR), col = freq*100), data = lm_fam, col = freq) + 
  geom_point(aes(size = freq))+ scale_color_viridis_c( direction = -1 ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 2, col = "Red", lty = "dashed", alpha = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "top")


pval_plot <- pval_plot + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/disease_otu_association.pdf", 
    family = "ArialMT", useDingbats = F, width = 3.5)
pval_plot
dev.off()

####################################################################################
# albugo load relationship to disease index?
####################################################################################
shared_bac <- my_bac[,which(colnames(my_bac) %in% plant_clim$clim_data$Plant_ID)]
shared_oom <- my_oom[,which(colnames(my_oom) %in% plant_clim$clim_data$Plant_ID)]
shared_org <- t(rbind(my_bac, my_oom[,colnames(my_bac)]))
plant_shared <- plant_clim$clim_data[which(plant_clim$clim_data$Plant_ID %in% rownames(shared_org)),]
shared_org <- shared_org[as.character(plant_shared$Plant_ID),]
shared_org <- shared_org[,which(is.na(colSums(shared_org))==FALSE)]

lm_rec<- data.frame(otu = colnames(shared_org), pval = c(NA)*395)
plant_shared$disease = c(plant_shared$Disease_RL %in% c(5))
for(i in 1:dim(shared_org)[2]){
  df = data.frame(OTU = shared_org[,i], Lat =  plant_shared$disease) 
  colnames(df) = c("OTU", "Lat")
  model <- lm( OTU~Lat, data = df)
  lm_rec[i,2] = summary(model)$coefficients[2,4]
  lm_rec$freq[i] = sum(shared_org[,i])/395
}

lm_rec$FDR <- p.adjust(lm_rec$pval, method = "BH")

lm_rec$sig_0.01 <- lm_rec$FDR<0.01

lm_rec <- lm_rec %>% filter(freq > 0.001)



# Most significant associations with Lat come from Sphingomonadaceae

pval_plot <- ggplot(aes(x = otu, y = -log10(FDR), col = freq*100), data = lm_rec, col = freq) + 
  geom_point(aes(size = freq))+ scale_color_viridis_c( direction = -1 ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 2, col = "Red", lty = "dashed", alpha = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "top")


albugo_load <- pval_plot + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)




####################################################################################
# HpA load relationship to disease index?
####################################################################################
shared_bac <- my_bac[,which(colnames(my_bac) %in% plant_clim$clim_data$Plant_ID)]
shared_oom <- my_oom[,which(colnames(my_oom) %in% plant_clim$clim_data$Plant_ID)]
shared_org <- t(rbind(my_bac, my_oom[,colnames(my_bac)]))
plant_shared <- plant_clim$clim_data[which(plant_clim$clim_data$Plant_ID %in% rownames(shared_org)),]
shared_org <- shared_org[as.character(plant_shared$Plant_ID),]
shared_org <- shared_org[,which(is.na(colSums(shared_org))==FALSE)]

lm_rec<- data.frame(otu = colnames(shared_org), pval = c(NA)*395)
plant_shared$disease = c(plant_shared$Disease_RL %in% c(2))
for(i in 1:dim(shared_org)[2]){
  df = data.frame(OTU = shared_org[,i], Lat =  plant_shared$disease) 
  colnames(df) = c("OTU", "Lat")
  model <- lm( OTU~Lat, data = df)
  lm_rec[i,2] = summary(model)$coefficients[2,4]
  lm_rec$freq[i] = sum(shared_org[,i])/395
}

lm_rec$FDR <- p.adjust(lm_rec$pval, method = "BH")

lm_rec$sig_0.01 <- lm_rec$FDR<0.01

lm_rec <- lm_rec %>% filter(freq > 0.001)



# Most significant associations with Lat come from Sphingomonadaceae

pval_plot <- ggplot(aes(x = otu, y = -log10(FDR), col = freq*100), data = lm_rec, col = freq) + 
  geom_point(aes(size = freq))+ scale_color_viridis_c( direction = -1 ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 2, col = "Red", lty = "dashed", alpha = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "top")


hpa_load <- pval_plot + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/albugo_hpa_association.pdf", 
    family = "ArialMT", useDingbats = F, width = 11)
plot_grid(albugo_load, hpa_load, ncol = 1)
dev.off()
