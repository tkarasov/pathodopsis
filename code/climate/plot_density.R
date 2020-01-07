################################################################################
# Script to generate a denstiy map from geographic points
# author: moisesexpositoalonso@gmail.com
################################################################################

library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(raster)
library(devtools)
library(RColorBrewer)
devtools::install_github('MoisesExpositoAlonso/moiR')
library(moiR)


#### load bioclim ####


Euroclim<-make_Euroclim()
fieldenvs<-make_fieldenvs()




#### get the arabidopsis data ####

data(gbif)
gbif<- filter(gbif,countrycode !='IE') # Ireland
gbif<- filter(gbif,countrycode !='RU') # Russia
gbif<- filter(gbif,countrycode !='BY') # Belarus
gbif<- filter(gbif,countrycode !='UA') # Ukraine
gbif<- filter(gbif,countrycode !='KZ') # Kazakhstan

plotbaseraster()
points(gbif$latitude ~ gbif$longitude,col='white')

data(acc)
acc<-filter(acc, country !='RUS')

coords<- rbind(gbif[,c('longitude','latitude')], acc[,c('longitude','latitude')])
coords$longitude <- fn(coords$longitude)
coords$latitude <- fn(coords$latitude)


#### estimate density ####

exampleraster=Euroclim[[1]]
baseraster=Euroclim[[1]]
saveRDS(file="dataint/baseraster.rda",baseraster)
baseraster<-makebaseraster(exampleraster)
# moiR::envirplot(baseraster,discrete = T)

density<-getdensityraster(coords$longitude, coords$latitude, refraster=baseraster,dilute = 50,method='bilinear')
moiR::envirplot(density)

densityall <- focal(density, w=matrix(1, 31, 31), mean)
moiR::envirplot(densityall)

density1001<- focal(getdensityraster(fn(acc$longitude), fn(acc$latitude), refraster =baseraster,dilute = 50,method='bilinear'), w=matrix(1, 11, 11), mean) 
density1001<-maskpresence(density1001,baseraster,1)

examplemasked<-maskpresence(density1001,densityall, quantile(values(densityall),0.5,na.rm=T) )
pdf('figs/geographic_density_map.pdf')
plotbaseraster()
moiR::envirplot(examplemasked,add=T,vecol = mypalettes('jet'))
points(dry$longitude , dry$latitude,col=('white'),cex=0.2)
dev.off()

saveRDS(file = 'maps/aramask.rda',examplemasked)