#!/usr/bin/env Rscript


library(maps)
library(mapdata)
library(plotrix)
library(ggmap)
library(RColorBrewer)
#Sites from which we have microbiome samples (AM_b)	
#	Different color if PP_n <20
# Different symbol by PP_type

locations = read.table("~/Dropbox/pathodopsis/maps/22.5.2018.sampling_info.txt", header=T, sep="\t")
locations$microb_pop = locations$PP_n>0
#sb = scalebar(9.05,48.45, 0.005, , 0.0143, "km" )


map <-  ggmap(get_googlemap(c(0,44), zoom = 3,  style = 'feature:administrative.country|element:labels|visibility:off'))

#microbiome samples
microb = locations[locations$AM_b>0,]

#map + geom_point(data=locations, aes(x=Longitude, y= Latitude, colour="red")) + ylim(20,70) + xlim(-30,50)

#plot size of population information
pdf("~/Dropbox/pathodopsis/maps/ps_population_size.pdf", width=8, height=8)
map + geom_point(data=locations, aes(x=Longitude, y= Latitude, colour=microb_pop)) + ylim(20,70) + xlim(-30,50) + scale_color_brewer(palette="Dark2")
dev.off()

#plot the type of the population
pdf("~/Dropbox/pathodopsis/maps/population_type.pdf", width=8, height=8)
map + geom_point(data=locations, aes(x=Longitude, y= Latitude, colour=as.factor(PP_type))) + ylim(20,70) + xlim(-30,50) + scale_color_brewer(palette="Dark2")
dev.off()

#plot whether hpa is present
pdf("~/Dropbox/pathodopsis/maps/hpa_presence.pdf", width=8, height=8)
map + geom_point(data=locations, aes(x=Longitude, y= Latitude, colour=as.factor(H_b))) + ylim(20,70) + xlim(-30,50) + scale_color_brewer(palette="Dark2")
dev.off()

#pop size with altitude
pdf("~/Dropbox/pathodopsis/maps/pop_size.pdf", width=8, height=8)
map + geom_point(data=locations, aes(x=Longitude, y= Latitude, colour=log10(Pop_size))) + ylim(20,70) + xlim(-30,50) + scale_color_gradientn(colours = topo.colors(10))
dev.off()



