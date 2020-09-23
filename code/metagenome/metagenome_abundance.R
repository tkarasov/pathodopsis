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
load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds')
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
sample_data(GP1000)$Lat <- as.numeric(as.character(sample_data(GP1000)$Lat))
sample_data(GP1000)$Long <- as.numeric(as.character(sample_data(GP1000)$Long))
GP1000_dim <- subset_samples(GP1000, sample_data(GP1000)$samples.out %in% colnames(my_bac))
reorder <- as.character(sample_data(GP1000_dim)$samples.out)
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