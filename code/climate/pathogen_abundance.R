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
my_phylo <- readRDS("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/seqtab_final.rds")
my_phylo_tax <- readRDS('/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/tax_final.rds')

pathogen <- read.table(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/pathogen_list_bartoli.txt",
  header = TRUE)

####################
# subset to ASVs corresponding to pathogen list starting with genus
####################
OTU5_match=read.table(
  "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/OTU5_ident.blast.out",
  header =FALSE)
more99 = OTU5_match[OTU5_match$V3>=99.0,]
keep_asvs <- more99$V2
genera <- as.character(unique(pathogen$Genus))
species <- as.character(unique(pathogen$Species))
pres_species <- tax_table(GP1000)[,7][which(tax_table(GP1000)[,7] %in% species)]
pres_genera <- tax_table(GP1000)[,6][which(tax_table(GP1000)[,6] %in% genera)]
xanth <- subset_taxa(GP1000, Genus=="Xanthomonas")
ralst <- subset_taxa(GP1000, Genus=="Ralstonia")
otu5 <- subset_taxa(GP1000, rownames(tax_table(GP1000))=="seq_10")
all_otu5 <- subset_taxa(GP1000, rownames(tax_table(GP1000))%in%keep_asvs)

#subset phyloseq to those identified as pathogens
pathogens_ph <- subset_taxa(GP1000, rownames(tax_table(GP1000)) %in% rownames(pres_species))


####################
# Build Lat/Long Dataframe
####################
sample_data(GP1000)$Lat <- as.numeric(as.character(sample_data(GP1000)$Lat))
sample_data(GP1000)$Long <- as.numeric(as.character(sample_data(GP1000)$Long))

spec_asv <- data.frame(sample_data(GP1000))
spec_asv$Xanthomonas = c(otu_table(xanth)[,1])
spec_asv$Ralstonia = c(otu_table(ralst)[,1])
spec_asv$otu5 = c(otu_table(otu5)[,1])
spec_asv$all_otu5 = rowSums(otu_table(all_otu5))

plant_asv = spec_asv %>% filter(Host_Species=="Ath")
soil_asv = spec_asv %>% filter(Used_for == "SOIL")


ggplot(data=plant_asv, aes(x=Long, y=Lat, col=log10(otu5+1)))+
  geom_point()+
  scale_color_viridis_c()
####################
# Build maps of abundance of pathogens
####################
####################################################################################
# Read in data and project data and europe to laea centered on europe
####################################################################################


# path = "/ebio"
# OTU5_tab = read.table(paste(path, 
#                             "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/OTU5_RA_table.txt", sep =""), 
#                       sep="\t", 
#                       header = T )
# OTU5_tab$OTU5 = as.numeric(as.character(OTU5_tab$OTU5))
# 
# #remove NA's and Capsella
# OTU5_tab = OTU5_tab[!grepl("PC", rownames(OTU5_tab)),]
# 
# # Percentage of individuals that have OTU5:
# has_OTU5 = length(which(OTU5_tab$OTU5!=0))/length(OTU5_tab$OTU5)
OTU5_tab <- plant_asv

#First add some wiggle to the OTU5 points then reproject the points
OTU5_tab$Long = OTU5_tab$Long + runif(n=length(OTU5_tab$Long), min=0, max=0.000001)
OTU5_tab$Lat = OTU5_tab$Lat + runif(n=length(OTU5_tab$Lat), min=0, max=0.000001)
OTU5_tab <- OTU5_tab %>% filter(is.na(Lat)==FALSE) %>% filter(is.na(Long)==FALSE)
#Now make sure everything has the same projection
coordinates(OTU5_tab) = c("Long", "Lat")
proj4string(OTU5_tab) = CRS(my.transform)
my.sf = st_as_sf(x = OTU5_tab, coords = c("Long", "Lat"), crs = my.transform)
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

krig_it<-function(vkriging_result = autoKrige(formula = otu5~1, input_data = sp.laea, new_data = grd.laea)ar){
  
  #kriging_result1 = autoKrige(formula = otu5~1, input_data = sp.laea, new_data = europe.laea)
  #Krig = kriging_result$krige_output
  
  # take only the points falling in polygons. This is the step where we limit the infomration
  Krig = kriging_result$krige_output[!is.na(over(kriging_result$krige_output,as(europe.laea,"SpatialPolygons"))),]  
  Krig.fin = spTransform(Krig, my.transform)
  Krig_df = as.data.frame(Krig.fin)
  names(Krig_df) = c(  "APPT_pred","APPT_var","APPT_stdev", "longitude","latitude")
  my.df = as.data.frame(my.sp)
  
  #Make a dataframe of Kriging in which the poorly interpolated points are just white
  the_mode = getmode(kriging_result$krige_output@data$var1.pred)
  Krig.fin = kriging_result
  Krig.fin$krige_output@data$var1.pred[which(Krig.fin$krige_output@data$var1.pred==the_mode)] = NA
  
}
kriging_result = autoKrige(formula = otu5~1, input_data = sp.laea, new_data = grd.laea)
#kriging_result1 = autoKrige(formula = otu5~1, input_data = sp.laea, new_data = europe.laea)
#Krig = kriging_result$krige_output

# take only the points falling in polygons. This is the step where we limit the infomration
Krig = kriging_result$krige_output[!is.na(over(kriging_result$krige_output,as(europe.laea,"SpatialPolygons"))),]  
Krig.fin = spTransform(Krig, my.transform)
Krig_df = as.data.frame(Krig.fin)
names(Krig_df) = c(  "APPT_pred","APPT_var","APPT_stdev", "longitude","latitude")
my.df = as.data.frame(my.sp)

#Make a dataframe of Kriging in which the poorly interpolated points are just white
the_mode = getmode(kriging_result$krige_output@data$var1.pred)
Krig.fin = kriging_result
Krig.fin$krige_output@data$var1.pred[which(Krig.fin$krige_output@data$var1.pred==the_mode)] = NA

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
  geom_raster(data = Krig_df, 
              aes(x = Krig_df$longitude, 
                  y = Krig_df$latitude, fill = NULL
              ), 
              inherit.aes = T)

#lets just make a raster of the Krig_df
OTU5_point <-
  base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = log10((otu5+1)/10)), 
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


krig_OTU5 <-
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
                       labs(size="OTU5 RA(%)"), 
                       breaks = seq(0,1,.05)) +
  theme_opts +
  #theme_pubs  +
  xlab(label = NULL) +
  ylab(label = NULL)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_OTU5.pdf", family = "ArialMT", useDingbats = F)
otu5_grid <-plot_grid(OTU5_point, krig_OTU5, nrow = 2)
otu5_grid
dev.off()

save(otu5_grid, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_OTU5.rds" )