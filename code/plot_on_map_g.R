##### This is an R script to plot the GPS points on the map############################
#original code written by Gautam (copied by Talia 1/2019) then modified 

library(maps)
library(mapdata)
library(scales)
library(mapproj)
library(intrval)
colfunc<-colorRampPalette(c("dodgerblue2","khaki","orangered")) # create colours
mycol <- function(x, myrange, n=100) round( 1+(x-myrange[1])/diff(myrange) * (n-1))

###whole range map#####
pdf("~/Dropbox/pathodopsis/pathodopsismicrobiome/figures/pathodopsis_collections_locations.pdf", width = 10 , height=8)
map(database = "world",xlim = c(-15, 50), ylim = c(35, 65),col="#f0e4cd", fill= TRUE, resolution = 0, lty=1, lwd=0.2)
map.axes()
#map.scale()
##############
#import a data frame from .csv file
test5 <- read.table("~/Dropbox/pathodopsis/pathodopsismicrobiome/data/v1_Master_Pathodopsis_M_2018-07-20.csv", header = T, sep=",")
#("#003893","#fcd116","#ce1126","#0b8a46","#1ec5e4","#edbd39","#2162c8","#f3a782","#14b4c6","#145bc6","#14c67f","#c6145b","#767676","#d82129","#db977f","#1f174b","#dfa196","#ae030e","#8a6f7e","#d69031")
##############
###Cluster color is modified here, you can comment out next line if you don't want to customize your color i.e colors will be selected from the pallete above)
#cluster_color_colors <- c("#CC0744", "#111111", "orangered", "#00846F", "#7B4F4B", "purple4","#A1C299" ,"deepskyblue4", '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
cluster_color_colors <- c("#900c3f", "#204838", "#1f78b4", "#834848", "#505e5e", "#e31a1c", "#6a3d9a", "#b15928", "#ff7f00", "#35978f", "#d53e4f", "#bf812d", "#1a1a1a", "#4575b4", "#5aae61", "#dd3497", "#63531a", "#08306b", "#8c6bb1", "#d7301f")
points(x=test5$V3,y=test5$V2, pch=19,cex=1, col=(cluster_color_colors))
#points(x=test5$V3,y=test5$V2, pch=as.integer(test5$V6),cex=1, col=(cluster_color_colors))

#adjust transparency###
colorList <- adjustcolor(cluster_color_colors[test5$V1], alpha.f=0.9)

###add points from all the sampled sites######
sampled <- data.frame(cbind(cbind(test5$Lat, test5$Long), as.character(test5$Plant_ID)))
colnames(sampled)=c("Latitude", "Longitude", "Plant_ID")
###Transparency#####
colorGray <- adjustcolor("black", alpha.f=0.6)
###Points#####
points(x=sampled$Longitude  , y=sampled$Latitude, cex=1, pch=20, bg=colorGray, col="dodgerblue4", lwd=1)

dev.off()

###Read in table
metagenome=read.table("~/work_main/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/2018_9_metagenome_reads/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)

map(database = "world",xlim = c(-15, 50), ylim = c(35, 65),col="#f0e4cd", fill= TRUE, resolution = 0, lty=1, lwd=0.2)
map.axes()
Pseudo=data.frame(t(metagenome[c("Pseudomonadaceae"),]))
#Pseudo$Latitude=sampled[sampled$Plant_ID==rownames(Pseudo),]$Latitude
rownames(sampled)=sampled$Plant_ID
Pseudo$Latitude=as.numeric(as.character(sampled[rownames(Pseudo), ]$Latitude))
Pseudo$Longitude=as.numeric(as.character(sampled[rownames(Pseudo), ]$Longitude))
#points(x=Pseudo$Longitude  , y=Pseudo$Latitude, cex=1, pch=16, bg=colorGray, col = colfunc(12000)[round(Pseudo$Pseudomonadaceae)*2000+1], lwd=1)
points(x=Pseudo$Longitude  , y=Pseudo$Latitude, cex=1, pch=16, bg=colorGray, lwd=1)
legend("topleft",title="",legend=c(0:5),col = colfunc(12000)[c(1:6)*2000-1999], pch=20)
