library(automap)
library(fpc)
library(dplyr)
library(sf)
library(cluster)
library(vegan)
library(pvclust)
library(rgdal)
library(raster)
library(spData)
library(spDataLarge)
library(wesanderson)
library(cowplot)
# library(NbClust)
library(dendextend)


#The goal of this script is to identify clusters in A. thaliana data. With this script we identified two major OTU clusters via kmeans clustering


load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

#Website is amazing for options for determining the number of clusters: https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters

######################
# How many ASVs to consider (to reduce complexity). Consider as many as constitute an average of 50% of microbiome
######################

p = plant_clim$otu_table/1000
p_order = p[,order(colSums(-p, na.rm=T))]


plot_dist <- function(p_order, n){
  keep = p_order[,c(1:n)]
  mean_ = mean(rowSums(keep))
  return(mean_)
}

dist_vector <- c(dim(p_order)[1])

for(i in 1:dim(p_order)[2])
  dist_vector[i] = plot_dist(p_order, i)

for(i in 1:length(dist_vector)){
  print(i)
  if(dist_vector[i]<0.5){
  print("pass")
  }
  if(dist_vector[i]>=0.5){
  print(paste("The 50% value is ", i, sep=""))
  break
}
}

#subset the otu_table to the ASVs that constitute 50%
#plant_clim$otu_table <- plant_clim$otu_table[,colnames(p_order)[1:i]]

######################
# What is the optimal number of clusters?
######################
set.seed(16)

otu_scale1 <- plant_clim$otu_table

#now we do a hellinger transformation of the scaled data
otu_scale <- sqrt(otu_scale1/1000)

######################
# Elbow plot (does not coverge)
######################

wss <- (nrow(otu_scale)-1)*sum(apply(otu_scale,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(otu_scale,
                                     centers=i)$withinss)
plot(1:50, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

######################
# Silhouette says two clusters
######################

#Silouette: 2 clusters
pamk.best <- pamk(otu_scale)
my.pam <- pam(data.frame((otu_scale)), pamk.best$nc)
#For some reason the following throws an error
#plot(pam(data.frame((otu_scale)), pamk.best$nc))
asw <- numeric(20)

for (k in 2:20)
  asw[[k]] <- pam(otu_scale, k) $ silinfo $ avg.width

pamkzk.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", pamkzk.best, "\n")

# ######################
# # Calinsky says 2 
# ######################
# fit <- cascadeKM(scale(otu_scale, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)
# calinski.best <- as.numeric(which.max(fit$results[2,]))
# cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")


######################
# So let's do k-means clustering with two clsuters
######################
# nc <- NbClust(otu_scale, diss=NULL, distance = "euclidean", min.nc=2, max.nc=6, 
#       method = "kmeans", index = "silhouette")
# fviz_nbclust(otu_scale, kmeans, method = c("silhouette", "wss",
#    "gap_stat"))


set.seed(16)
clusters <- kmeans(otu_scale, 2, nstart = 100)

plant_clim$clim_data$cluster <- clusters$cluster
save(plant_clim, file = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

# and now look at the silhouette score for the two clusters
dis = dist(otu_scale)^2 # k-means clustering uses squared euclidean distances
sil = silhouette(clusters$cluster, dis)
pdf("//ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/silhouette_scored_poo2.pdf",
    useDingbats = FALSE, font = "ArialMT", width = 7, height  = 6)
plot(sil) 
dev.off()

######################
# How about hierarchical clustering
######################
d <- (otu_scale) %>% dist()
hc <- hclust(d)
hc$labels <- clusters$cluster
cl_members <- cutree(tree = hc, k = 3)
#plot(hc, hang = -1, cex = 0.6,  leaflab = "none")
plot(x = hc, labels =  row.names(hc), cex = 0.5)
rect.hclust(tree = hc, k = 3, which = 1:3, border = 1:2, cluster = cl_members, )
plant_clim$clim_data$hc_cuttree2 = cl_members
dend = as.dendrogram(hc) %>% set("leaves_colors", hc$labels) %>% set("leaves_pch", 16) 
#labels_colors(dend) <- hc$labels

ggplot(dend, labels = FALSE) 
#choosing optimal clusters supported by data
#https://www.datanovia.com/en/lessons/computing-p-value-for-hierarchical-clustering/
# pv <- parPvclust(cl=NULL, t(otu_scale), method.hclust = "average",
#              method.dist = "correlation", nboot = 100, 
#              iseed = NULL)

# # Default plot
# plot(pv, hang = -1, cex = 0.5)
# pvrect(pv)
# clusters <- pvpick(pv)
# 
# # average silhouette score for hclust
# h.cut <- hcut((otu_scale), k = 2, hc_method = "complete")
# fviz_silhouette(h.cut)
# fviz_cluster(h.cut)


######################
# And Let's plot the MDS with the clusters
######################

MDS <- sqrt(plant_clim$otu_table/1000) %>% dist() %>% cmdscale(eig = TRUE,  k = (dim(plant_clim$otu_table)[1]-1))
MDS.points <- data.frame(MDS$points)
colnames(MDS.points)[c(1:2)] = c("MDS1", "MDS2")
exp3 <-  ((MDS$eig) / sum(MDS$eig))[1]*100
exp4 <-  ((MDS$eig) / sum(MDS$eig))[2]*100

MDS_plot_kmeans <- 
  ggplot(data = MDS.points, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(col=as.factor(clusters$cluster)), cex = 3, alpha = 0.5) +
  scale_colour_manual(name = "Cluster", values = c('#d95f02', '#1b9e77')) +
  theme_bw() +
  xlab(paste(paste("MDS1 (", round(exp3), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp4), sep=""),"%)",sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey70") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey70") +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.9),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )

MDS_plot_hclust <- 
  ggplot(data = MDS.points, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(col=as.factor(cl_members)), cex = 2.5, alpha = 0.5) +
  scale_colour_brewer(name = "Cluster", palette = c('#d95f02', '#1b9e77')) +
  theme_bw() +
  xlab(paste(paste("MDS1 (", round(exp3), sep=""),"%)",sep="")) +
  ylab(paste(paste("MDS2 (", round(exp4), sep=""),"%)",sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey70") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey70") +
  theme(legend.justification=c(0,0), 
        legend.position=c(.7,.9),
        legend.title = element_blank(),
        legend.text.align = 0,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.box.background = element_rect(colour = "black")
  )

####################################################################################
# Read in data and project data and europe to laea centered on europe
####################################################################################
my.transform <- "+proj=longlat +datum=WGS84 +no_defs"
laea <- "+init=epsg:3857" 

myvar_tab = data.frame(myvar = clusters$cluster, myvar2 = cl_members, Longitude =  as.numeric(as.character(plant_clim$clim_data$Long)), 
                       Latitude = as.numeric(as.character(plant_clim$clim_data$Lat)))
myvar_tab = myvar_tab[complete.cases(myvar_tab),]

# Percentage of individuals that have myvar:
#has_myvar = length(which(myvar_tab$myvar!=0))/length(myvar_tab$myvar)

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
Krig = kriging_result$krige_output

# take only the points falling in polygons. This is the step where we limit the infomration
Krig = kriging_result$krige_output[!is.na(over(kriging_result$krige_output,as(europe.laea,"SpatialPolygons"))),]  
Krig.fin = spTransform(Krig, my.transform)
Krig_df = as.data.frame(Krig.fin)

####################################################################################
# Plot world map with data points
####################################################################################
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

xmin_plot = extent(my.sp)[1] -5
xmax_plot = extent(my.sp)[2] +5
ymin_plot = extent(my.sp)[3] -5
ymax_plot = extent(my.sp)[4] +5

pal <- wes_palette("Zissou1", 5, type = "discrete")

base_europe_df <-
  ggplot() +
  geom_raster(data = Krig_df, 
              aes(x = Krig_df$coords.x1, 
                  y = Krig_df$coords.x2), fill = "White", 
              inherit.aes = T)

#lets just make a raster of the Krig_df
myvar_point_kmeans <-base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = as.factor(myvar)), 
              width = .5, cex = 1.5, alpha = 0.5) +
  scale_colour_manual(name = "Cluster", values = c('#d95f02', '#1b9e77')) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_plot, ymax_plot)) +
  scale_x_continuous(expand = c(0,0), limits = c(xmin_plot, xmax_plot)) +
  theme_opts + 
  #  theme_pubs +
  xlab(label = NULL) +
  ylab(label = NULL)


myvar_point_hclust <-base_europe_df + 
  geom_jitter(data = my.df, 
              aes(coords.x1, coords.x2, colour = as.factor(myvar2)), 
              width = .5, cex = 1.5, alpha = 0.5) +
  scale_colour_brewer(name = "Cluster", palette = "Dark2") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_plot, ymax_plot)) +
  scale_x_continuous(expand = c(0,0), limits = c(xmin_plot, xmax_plot)) +
  theme_opts + 
  #  theme_pubs +
  xlab(label = NULL) +
  ylab(label = NULL)

kmeans.p <- plot_grid(MDS_plot_kmeans + theme(legend.position = "none"), 
                    myvar_point_kmeans + theme(legend.position = "none"))

hclust.p <- plot_grid(MDS_plot_hclust + theme(legend.position = "none"), 
                    myvar_point_hclust + theme(legend.position = "none"))

cluster_plot <- kmeans.p #plot_grid(kmeans.p, hclust.p, labels = c("k-means", "hclust")
# , nrow = 2)
                          

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/map_clusters_poo.pdf", useDingbats = FALSE, font = "ArialMT", width = 7, height  = 3)

cluster_plot

dev.off()


