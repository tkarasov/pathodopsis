library(phyloseq)
library(cowplot)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyr)
library(OpenStreetMap)

# The goal of this script is to provide basic measures of plant status
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")

#########################################
# Plant health traits
#########################################
clim <- plant_clim$clim_data
clim$Lat <- as.numeric(as.character(clim$Lat))
clim$Long <- as.numeric(as.character(clim$Long))

world <- ne_countries(scale = "medium", returnclass = "sf")

base_map <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = c(-20, 45), ylim = c(30, 70), expand = FALSE) +
    theme_bw() +
    theme(legend.position="bottom") +
    xlab("") +
    ylab("") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
	  legend.position = c(0.2, 0.8),
	  legend.title.align=0.5,
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank())

# original script had something called clim_df. identity lost
clim_df <- clim
# full_base <- base_map +
#   geom_jitter(data = clim_df, 
#               aes(x = as.numeric(Long), y = as.numeric(as.character(Lat))), 
#               size = .75,
#               alpha =.5, colour = "RED") 



######## Better looking map
map <- openproj(openmap( upperLeft = c(70, -15), lowerRight = c(30, 50), zoom = 2, type = "nps"))
plot1 <- autoplot(map) 


plot2 <- plot1 + geom_jitter(data = clim_df, 
            aes(x = as.numeric(Long), y = as.numeric(as.character(Lat))), 
            size = .75,
            alpha =.5, colour = "steelblue4") +
            xlab("Longitude") +
            ylab("Latitude")


herbivory <- base_map + 
  geom_jitter(data = clim_df, 
              aes(x = as.numeric(Long), y = as.numeric(as.character(Lat)),  colour = Herbivory), 
              size = 2,
              alpha =.5) +
  scale_color_viridis_c("Herbivory State") 

disease <-  base_map + 
  geom_jitter(data = clim_df, 
              aes(x = as.numeric(Long), y = as.numeric(as.character(Lat)),  colour = Disease_RL), 
              size = 2,
              alpha =.5) +
  scale_color_viridis_c("Disease State") 

size <- base_map + 
  geom_jitter(data = clim_df, 
              aes(x = as.numeric(Long), y = as.numeric(as.character(Lat)),  colour = as.numeric(as.character(R_diameter))), 
              size = 2,
              alpha =.5) +
  scale_color_viridis_c("Rosette Diameter (cm)") 

dev <- base_map + 
  geom_jitter(data = clim_df, 
              aes(x = as.numeric(Long), y = as.numeric(as.character(Lat)),  colour = as.numeric(as.character(Developmental_state))), 
              size = 2,
              alpha =.5) +
  scale_color_viridis_c(name="Developmental State") 

variables <- clim_df %>% dplyr::select(c(Developmental_state, Disease_RL, Herbivory, R_diameter, Lat))
variables_gather <- reshape2::melt(variables, id.vars="Lat")
# Remove NA
variables_gather <- variables_gather %>% dplyr::filter(value!="NA")
  
pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/map_plant_health.pdf", useDingbats = FALSE, font = "ArialMT", width  = 7.2)
plot_grid(herbivory, disease, size, dev)
dev.off()

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/points_map.pdf", useDingbats = FALSE, family  = "ArialMT", width  = 2.8, height = 2.8)
#full_base + coord_sf(xlim = c(-15, 40), ylim = c(30, 65), expand = FALSE)
plot2 + theme_bw() +  theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            # Remove panel background
                            panel.background = element_blank())
dev.off()

n_fun <- function(x){
  val = stats::quantile(x)[4] + 1.25
  return(data.frame(y = as.numeric(val), 
                    label = paste0("n = ",length(x))))
}


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/barplot_plant_health.pdf", useDingbats = FALSE, family  = "ArialMT", width  = 4.6, height=4.5)
ggplot(variables_gather, aes(x=value, y =Lat)) +
  geom_boxplot(colour="steelblue4") +
  facet_wrap(~variable) +
  theme_bw() + geom_jitter(alpha = 0.1, width=0.1, color = "steelblue4") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank()) +
    stat_summary(
      geom = "text",
      fun.data = n_fun,
      family = "ArialMT",
      size = 2) +
  ylab("Latitude")
dev.off()

#########################################
# Now OTUs and plant health
########################################
# OTU <- (as.data.frame(plant_clim$otu_table))
#OTU$disease <- as.numeric(disease)
# disease <- plant_clim$clim_data$Disease_RL
# corr <- function(x){return(cor.test(x, OTU$disease, na.rm =TRUE))}
# corr_vec <- data.frame(disease = numeric(dim(OTU)[1]), herbivory = numeric(dim(OTU)[1]), development = numeric(dim(OTU)[1]))#, herbivory = c(), development = c())
# for(i in 1:dim(OTU)[2]){
#  corr_vec$disease[i] <- cor.test(log10(OTU[,i] + 1), plant_clim$clim_data$Disease_RL, na.rm=TRUE)$estimate
#  corr_vec$herbivory[i] <- cor.test(log10(OTU[,i] + 1), plant_clim$clim_data$Herbivory, na.rm=TRUE)$estimate
#  corr_vec$development[i] <- cor.test(log10(OTU[,i] + 1), plant_clim$clim_data$Developmental_state, na.rm=TRUE)$estimate
# }
# OTU_corr <-OTU %>% mutate_all(corr)
  #sapply(OTU, function(x)cor.test(c(x), disease, na.rm = TRUE))

# plot_bar(plant_clim, "Family", fill="Genus", facet_grid=~Disease_RL)OA
