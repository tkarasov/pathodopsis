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
#remove all NAs
fin <- na.omit(fin)
otu_keep <- otu_keep[rownames(fin),]

# Now do a lmer
mo1 <- glm(as.factor(cluster)~ V3+V4+V5+V6+V7+Lat+PDSI, family=binomial(link='logit'), dat=fin)
capture.output(summary(mo1))

#[3] "glm(formula = as.factor(cluster) ~ V3 + V4 + V5 + V6 + V7 + Lat + "
# [4] "    PDSI, family = binomial(link = \"logit\"), data = fin)"        
# [5] ""                                                                  
# [6] "Deviance Residuals: "                                              
# [7] "    Min       1Q   Median       3Q      Max  "                     
# [8] "-3.0734  -0.2942  -0.1408   0.3164   2.9147  "                     
# [9] ""                                                                  
#[10] "Coefficients:"                                                     
#[11] "            Estimate Std. Error z value Pr(>|z|)    "              
#[12] "(Intercept) 13.65865    2.88770   4.730 2.25e-06 ***"              
#[13] "V3          -7.06332    9.21068  -0.767  0.44316    "              
#[14] "V4           2.69063   11.75597   0.229  0.81897    "              
#[15] "V5          13.66086   16.65738   0.820  0.41215    "              
#[16] "V6           3.91046    5.67368   0.689  0.49068    "              
#[17] "V7          -0.56140    6.04698  -0.093  0.92603    "              
#[18] "Lat         -0.30019    0.05711  -5.256 1.47e-07 ***"              
#[19] "PDSI        -0.30717    0.11153  -2.754  0.00589 ** "              
#[20] "---"                                                               
#[21] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"    
#                                                  
#[23] "(Dispersion parameter for binomial family taken to be 1)"                                                                        
#[25] "    Null deviance: 480.04  on 357  degrees of freedom"             
#[26] "Residual deviance: 186.78  on 350  degrees of freedom"             
#[27] "AIC: 202.78"                                                                                                                       
#[29] "Number of Fisher Scoring iterations: 7"                            
                                      





##data(dune,dune.env) 
mod<-cca(PCo1_loading~V3+V4+V5+V6+V7+Lat+PDSI,data =fin) 
##overalltest anova(mod) 
##testsforindividualterms anova(mod,by="term") anova(mod,by="margin") 
##sequentialtestforcontrasts anova(mod,by="onedf") 
##testforaddingallenvironmentalvariables anova(mod,cca(dune~.,dune.env))

# what I want to do is to measure the relative contributions of latitude, and other variables on the PCo1 and PCo2 loadings
# next step is to do the bray curtis distance. then do adonis anova.
# convert OTU keep to RA
otu_RA <- otu_keep/rowSums(otu_keep)
otu_RA <- otu_RA[rownames(fin),]
otu_dist <- (vegdist(otu_RA, method = "bray"))
ad <- adonis(otu_dist~V3+V4+V5+V6+V7+Lat+PDSI,data =fin, permutations = 999, method = "bray",
       strata = NULL, contr.unordered = "contr.sum",
       contr.ordered = "contr.poly")
