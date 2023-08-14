# The goal of this script is to take the eigenvectors output from gcta and run anovas in vegan

library(vegan)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(lme4)

eigen_vec <- read.table("test.eigenvec", row.names =1 )
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
clim_data <- plant_clim$clim_data
otu_data <- plant_clim$otu_table
clim_keep <- clim_data[clim_data$Sequence_ID %in% rownames(eigen_vec),]
otu_keep <- otu_data[clim_data$Sequence_ID %in% rownames(eigen_vec),]
clim_keep <- clim_keep[clim_keep$Sequence_ID %in% rownames(phenotype3),]
otu_keep <- otu_keep[clim_keep$Sequence_ID %in% rownames(phenotype3),]
eigen_vec1 <- eigen_vec[order(clim_keep$Sequence_ID),]
eigen_vec1 <- cbind(eigen_vec1, clim_keep[, c("Lat", "Date", "PDSI", "cluster")])
colnames(eigen_vec1)[1]="Sequence_ID"
phenotype3 <- read.table("phenotype3.txt", row.names =1)
pheno <- phenotype3[order(clim_keep$Sequence_ID),]
colnames(pheno) <- c("Sequence_ID", "PCo1_loading")
PCo1_loading <- pheno$PCo1_loading
fin <- cbind(eigen_vec1,PCo1_loading)
fin <- na.omit(fin)

# Now do a lmer
mo1 <- lm(cluster~ V3+V4+V5+V6+V7+Lat+PDSI), dat=fin)
capture.output(summary(mo1))

##data(dune,dune.env) 
mod<-cca(PCo1_loading~V3+V4+V5+V6+V7+Lat+PDSI,data =fin) 
##overalltest anova(mod) 
##testsforindividualterms anova(mod,by="term") anova(mod,by="margin") 
##sequentialtestforcontrasts anova(mod,by="onedf") 
##testforaddingallenvironmentalvariables anova(mod,cca(dune~.,dune.env))

# what I want to do is to measure the relative contributions of latitude, and other variables on the PCo1 and PCo2 loadings
# next step is to do the bray curtis distance. then do adonis anova.



distance_matrix <- dist(otu_clim$)
ad <- adonis(PCo1_loading~V3+V4+V5+V6+V7+Lat+PDSI,data =fin, permutations = 999, method = "bray",
       strata = NULL, contr.unordered = "contr.sum",
       contr.ordered = "contr.poly"
