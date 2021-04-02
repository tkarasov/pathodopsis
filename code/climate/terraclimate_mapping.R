#Plotting PDSI (Figure 4)
# This script takes terraclimate PDSI data and graphs. Used for making Figure 4
#got script from here

library(rgdal)
library(raster)
library(ggplot2)
library(forcats)
library(nnet)
a <- "/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/TerraClimate_PDSI_2018.nc"
pdsI <- stack(a)
pdsI <- subset(pdsI, 1:4)
pdsI_mean <- calc(pdsI, # RasterStack object
                 fun = mean, # Function to apply across the layers
                 na.rm = TRUE)

ext <- extent(c(xmin = -15, xmax = 54, 
                ymin = 32, ymax = 68))

pdsI_crop <- crop(x = pdsI_mean, y = ext )


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/pdsI_map_january_April_2018.pdf",family = "ArialMT", useDingbats = F)
plot(pdsI_crop, xlim = c(-15,54), 
     axes = FALSE, 
     ext = ext)
dev.off()

# Box plots for PDSI and cluster assignment
load("~/work_main/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
clim_data <- plant_clim$clim_data
clim_data$cluster <- as.factor(clim_data$cluster)
clim_data$cluster <- fct_relevel(clim_data$cluster, "2", after = 2)
clim_data$Lat <- as.numeric(as.character(clim_data$Lat))

pdf("~/work_main/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/pdsI_boxplot.pdf",
    family = "ArialMT", useDingbats = F, width = 3.5, height = 3.5)
ggplot(data = clim_data, aes(x = cluster, y = PDSI)) +
  #geom_boxplot() +
  geom_jitter(color=clim_data$cluster, size=0.4, alpha=0.9,width = .2) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  theme_classic() +
  xlab("Cluster")
  dev.off()

model <- multinom(cluster ~ Lat + PDSI, data = clim_data, family = )
zvalues <- summary(model)$coefficients / summary(model)$standard.errors
