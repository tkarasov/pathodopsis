# This script does Kriging for the abundance of OTU5 in my dataset


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
library(proj4)
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
library(inrval)
Arial = font_add("Arial", "/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/pathodopsis/misc_tkarasov/fonts/Arial.ttf")



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

####################################################################################
# Read in data and project data and europe to laea centered on europe
####################################################################################

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
pseud_sphing<-plant_clim$clim_data
pseud_sphing$sphing <- plant_clim$otu_table[,1]/10
pseud_sphing$pseud <- plant_clim$otu_table[,7]/10


# Percentage of individuals that have OTU5:
# has_OTU5 = length(which(OTU5_tab$OTU5!=0))/length(OTU5_tab$OTU5)

#First add some wiggle to the OTU5 points then reproject the points
pseud_sphing$Longitude = pseud_sphing$Long + runif(n=length(pseud_sphing$Long), min=0, max=0.000001)
pseud_sphing$Latitude = pseud_sphing$Lat + runif(n=length(pseud_sphing$Lat), min=0, max=0.000001)

#Now make sure everything has the same projection
coordinates(pseud_sphing) = c("Longitude", "Latitude")
proj4string(pseud_sphing) = CRS(my.transform)
my.sf = st_as_sf(x = pseud_sphing, coords = c("Longitude", "Latitude"), crs = my.transform)
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
kriging_result = autoKrige(formula = pseud~1, input_data = sp.laea, new_data = grd.laea)
kriging_result1 = autoKrige(formula = pseud~1, input_data = sp.laea, new_data = europe.laea)
Krig = kriging_result$krige_output

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
pseud_point <-
  base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = log10(pseud)), 
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


krig_pseud <-
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
otu5_grid <-plot_grid(pseud_point, krig_pseud, nrow = 2)
otu5_grid
dev.off()

save(otu5_grid, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_OTU5.rds" )


####################################################################################
# Plot Variogram 
####################################################################################  

vgOK <- kriging_result$exp_var
vgOK$distance <- vgOK$dist/1000
fitted.variogram <- kriging_result$var_model
variogram.line <- variogramLine(fitted.variogram, max(vgOK$dist)) 

vg <- ggplot() +
  geom_point(data = vgOK, aes(x = (dist/1000), y = gamma), cex = 3) +
  geom_line(data = variogram.line, aes(x = (dist/1000), y = gamma, colour = "RED")) +
  xlab("Distance between points (km)") + 
  scale_x_log10() +
  scale_y_log10() +
  ylab("Semivariance") +
  theme_bw() +
  theme(legend.position = "none")

                              
plot(vgOK, fitted.variogram, xlab = "Meters")


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTU5_variogram.pdf", family = "ArialMT", useDingbats = F)
vg + theme_pubs
dev.off()

