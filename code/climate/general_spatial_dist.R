
library(phyloseq)
library(dplyr)
library(reshape2)

# The goal of this script is to take a phyloseq object and a list of microbes and to assess the occupancy-abundance relationships
#########################################
# Functions
#########################################

per_site <- function(GP, seq, comp){
  #Compiles data for comparisons between soil and species
  # seq must be a vector
  if(comp == "between"){
    ath = subset_samples(GP, Species == "Ath")
    ath = prune_taxa(c(seq), ath)
    ath = make_replicate(ath)
    fin = cbind(sample_data(ath), otu_table(ath))
    fin = reshape2::dcast(fin, Latitude ~ rep, value.var = eval(seq))
  }
  if(comp == "capsella"){
    ath = subset_samples(GP, Sample_type == "M")
    ath = prune_taxa(c(seq), ath)
    fin = cbind(sample_data(ath), otu_table(ath))
    fin = reshape2::dcast(fin, Latitude ~ Species, mean, value.value = eval(seq))
  }
  if(comp == "soil"){
    ath = subset_samples(GP, Species != "Cap")
    ath = prune_taxa(c(seq), ath)
    fin = cbind(sample_data(ath), otu_table(ath))
    fin = reshape2::dcast(fin, Latitude ~ Sample_type, mean, value.var = eval(seq))
  }
  return(fin)
}


make_replicate <- function(GP){
  #adds sample data column
  unique_loc = unique(sample_data(GP)$Latitude)
  sample_data(GP)$rep = c(0)
  for(val in unique_loc){
    chosen = which(sample_data(GP)$Latitude==val)
    sample_data(GP)$rep[chosen] = paste("rep", c(1:length(chosen)), sep="")
  }
  return(GP)
}

#########################################
# Pathogen lables
#########################################
OTU5 <- c("seq_10")
Xan <- "seq_49"
sphing <- "seq_1"


#########################################
# Load ITS and 16S data
#########################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")


#########################################
# Correlation between samples from same site
#########################################
ath = per_site(GP_at15_all, Xan, "between")
cap = per_site(GP_at15_all, Xan, "capsella")
soil = per_site(GP_at15_all, Xan, "soil")

tax_names <- colnames(otu_table(GP_at15_all))
ath_all <-sapply(tax_names, function(x) per_site(GP_at15_all, x, "between"))

p1 <- ggplot(data = ath, aes(x = (rep1/10 +1) , y = (rep2/10 +1))) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()
  

#########################################
# Correlation between A. thaliana and Capsella
#########################################


#########################################
# Correlation between A. thaliana and Capsella
#########################################
