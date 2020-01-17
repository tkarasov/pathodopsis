# This script does Kriging for the abundance of myvar in my dataset


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
library(intrval)
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
path = "/ebio"
myvar_tab = read.table(paste(path, "/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/myvar_RA_table.txt", sep =""), 
                      sep="\t", 
                      header = T )
# met_load = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_9_2018_reads.txt", sep ="\t", header = T)
# met_load$Total_load = rowSums(met_load[,c(15:dim(met_load)[2])], na.rm = T)
met_load = read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/Total_load_v2.csv", sep = ",", header = T, row.names = 1)
myvar_tab = met_load
myvar_tab$myvar = met_load$Total_load

myvar_tab$myvar = as.numeric(as.character(myvar_tab$myvar))

# Remove top and bottom 5% of values
top5 = sort(myvar_tab$myvar)[round(length(myvar_tab$myvar)*.95)]
bottom5 = sort(myvar_tab$myvar)[round(length(myvar_tab$myvar)*.05)]
myvar_tab = filter(myvar_tab, myvar<top5 & myvar>bottom5)


#remove NA's and Capsella
myvar_tab = myvar_tab[!grepl("PC", rownames(myvar_tab)),]

# Percentage of individuals that have myvar:
has_myvar = length(which(myvar_tab$myvar!=0))/length(myvar_tab$myvar)

#First add some wiggle to the myvar points then reproject the points
myvar_tab$Longitude = myvar_tab$Longitude + runif(n=length(myvar_tab$Longitude), min=0, max=0.000001)
myvar_tab$Latitude = myvar_tab$Latitude + runif(n=length(myvar_tab$Latitude), min=0, max=0.000001)

#Now make sure everything has the same projection
coordinates(myvar_tab) = c("Longitude", "Latitude")
proj4string(myvar_tab) = CRS(my.transform)
my.sf = st_as_sf(x = myvar_tab, coords = c("Longitude", "Latitude"), crs = my.transform)
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
kriging_result = autoKrige(formula = myvar~1, input_data = sp.laea, new_data = grd.laea)
kriging_result1 = autoKrige(formula = myvar~1, input_data = sp.laea, new_data = europe.laea)
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
myvar_point <-
  base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = log10(myvar*100)), 
              width = .5, cex = 2) +
  scale_colour_gradientn(colours = pal,
                         values = c(0, 0.1, 0.4, 1), name = "log10(RA(%))") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_plot, ymax_plot)) +
  scale_x_continuous(expand = c(0,0), limits = c(xmin_plot, xmax_plot)) +
  theme_opts + 
#  theme_pubs +
  xlab(label = NULL) +
  ylab(label = NULL)


####################################################################################
# Plot Kriging Interpolation
####################################################################################

krig_myvar <-
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
                       labs(size="myvar RA(%)"), 
                       breaks = seq(0,1,.05)) +
  theme_opts +
#  theme_pubs  +
  xlab(label = NULL) +
  ylab(label = NULL)

file_out = paste(paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/", myvar, sep=""), ".pdf", sep="")
pdf(file_out, family = "ArialMT", useDingbats = F)
plot_grid(myvar_point, krig_myvar)
dev.off()


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


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/myvar_variogram.pdf", family = "ArialMT", useDingbats = F)
vg + theme_pubs
dev.off()



# 
# ####################################################################################
# # Smooth over the region with kernel smoothing
# ####################################################################################
# my.polygons <- rasterToPolygons(my.raster)
# 
# p_smooth_ksmooth <- smoothr::smooth(my.polygons, method = "ksmooth", smoothness = 1)
# 


# ####################################################################################
# # Smooth over the region with nearest neighbor
# ####################################################################################
# # Tutorial on regional smoothing: https://pudding.cool/process/regional_smoothing/
# hex_points <- spsample(as_Spatial(europe), n = 100,type = "hexagonal", cellsize = 150)
# 
# #convert to a SPDF
# hex_polygons <- HexPoints2SpatialPolygons(hex_points) %>%
#   #clip the hexagon grid to the edges of Ireland
#   gIntersection(ireland_gdal, byid = TRUE) %>%
#   st_as_sf(crs = st_crs(ireland_sf)) %>%
#   st_transform(crs = st_crs(ireland_sf))
# 
# #find which small area shapefile (ireland_sf) fall within each hexagon
# bin_list <- st_intersects(hex_polygons, ireland_sf)
# 
# coords <- coordinates(my.sp)
# shapefile_df <- as(my.sp, "data.frame")
# IDs <- row.names(shapefile_df)
# knn50 <- knn2nb(knearneigh(coords, k = 50), row.names = IDs)
# knn50 <- include.self(knn50)   
# 
# # Creating the localG statistic for each of counties, with a k-nearest neighbor value of 5, and round this to 3 decimal places
# localGvalues <- localG(x = as.numeric(shapefile_df$myvar), listw = nb2listw(knn50, style = "B"), zero.policy = TRUE)
# localGvalues <- round(localGvalues,3)
# 
# # Create a new data frame that only includes the county fips codes and the G scores
# new_df <- data.frame(shapefile_df$GEOID)
# new_df$values <- localGvalues
# 
# #Huzzah! We are now ready to export this CSV, and visualize away!
# write.table(row.names = FALSE,new_df, file = "smooooothstuff.csv", sep=",")







####################################################################################
# Old Kriging text
####################################################################################
# #First estimate the thin plane spline
# my.Krig <- Krig(x = as.matrix(my.df[,c("coords.x1", "coords.x2")]), Y = as.numeric(my.df$myvar))
# set.panel(2,2)
# plot(my.Krig)
# set.panel()
# 
# #Rasterize the original dataset
# r <- raster(my.sp, res=100)
# 
# #Make a rasterlayer with interpolated values using the fitted model object from the thin plate spline
# #interpolate the abundance over r
# krig.int <- raster::interpolate(r, my.Krig)
# 
# #Predict the values using my.Krig
# krig.pred <- predictSurface(my.Krig)
# 
# #Mask
# krig.mask <- mask(krig.int, my.sp)
# krig.mask.df <- fortify(tps_mask_df)
# plot(krig.mask)
# 
# obs.pred <- data.frame(y = )
# #RMSE
# rmse(tps.pred$tps, tps.pred$y)


#Now I need to create a variogram. Like in lm values can be dependent on some feature like distance from the river (meuse dataset) and if there is no such variable, use 1
#The variogram function calculates the semivariance of the trait over space
# dt.vgm <- variogram(myvar~1, locations = coordinates(my.sp), my.sp)
# class(dt.vgm)
# 
# #Now build a model to fit the variogram. The first parameter is a sample variogram. The second is a model with parameters to be fit to the sample variogram.
# dt.fit <-fit.variogram(dt.vgm, model = vgm(model = c("Exp", "Mat", "Sph")),fit.kappa = T) # fit model
# 
# # vgm() list of models
# plot(dt.vgm, dt.fit)
# 
# #This is just nearest neighbor estimation. This doesn't seem to work Looks super weird.
# lzn.kriged0 <- krige((myvar) ~ 1, my.sp, grd_sp, set = list(method = "med"), nmax = 11)
# 
# # Now I can perform kriging and plot the result. But this takes many minutes and I got impatient.
# lzn.kriged1 <- krige((myvar) ~ 1, my.sp, grd_sp, model = dt.fit, maxdist = 1000)
#lzn.kriged %>% as.data.frame %>% rename(lon=coords.x1, lat=coords.x2) %>% 
#  ggplot(aes(x=lon, y=lat)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
#  scale_fill_gradient2(low="green", mid = "yellow",  high="red",midpoint = 15) +
#  theme_bw()+
#  geom_point(data=dt, aes(color=z), size=10)+
#  geom_text(data=dt, aes(label=z), color="white")

  
  # # create an empty raster object to the extent of the points
  # e = extent(my.sp)
  # r <- raster(e, ncol=100, nrow=100)
  # 
  # # rasterize your irregular points 
  # my.raster <- rasterize(my.sp@coords[,c(1,2)], r, my.sp@data$myvar, fun = mean)
  

#But weirdly this base map of europe turns out weird when I graph next to the myvar distribution.
# Create a base map of europe
# base_europe <- ggplot() + 
#   geom_sf(data = europe.sf, 
#           fill="white") +
#   xlim(xmin_plot, xmax_plot) +
#   ylim (ymin_plot, ymax_plot ) +
#   theme_opts   
# 
# 
# #Add myvar onto the european basemap
# myvar_point <- base_europe +
#   geom_jitter(data = my.df, 
#               aes(coords.x1, coords.x2, color = log10(myvar*100)), 
#               width = .5, cex = 2) +
#   scale_colour_gradientn(colours = pal,
#                          values = c(0, 0.1, 0.4, 1), name = "log10(RA(%))") +
#   theme_opts
# 
# myvar_point

  