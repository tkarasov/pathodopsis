# The goal of this script is to take the eigenvectors output from gcta and run anovas in vegan

library(vegan)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(lme4)
library(car)

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
all_pheno <- read.table("phen_MDS1_MDS2.txt", row.names = 1) 
pheno <- all_pheno[order(clim_keep$Sequence_ID),]
colnames(pheno) <- c("Sequence_ID", "PCo1_loading", "PCo2_loading")
PCo1_loading <- pheno$PCo1_loading
PCo2_loading <- pheno$PCo2_loading
fin <- cbind(eigen_vec1,PCo1_loading)
fin <- cbind(fin, PCo2_loading)
#remove all NAs
fin <- na.omit(fin)
fin2 <- fin[!rownames(fin)%in%c("PA0763", "PA0805"),] #these two samples were not found in the OTU file
otu_keep <- otu_keep[rownames(fin2),]

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
                                      

mo_MDS1 <- lm((PCo1_loading)~ V3+V4+V5+V6+V7+Lat+PDSI, dat=fin)
mo_MDS2 <- lm((PCo2_loading)~ V3+V4+V5+V6+V7+Lat+PDSI, dat=fin)

# > Anova(mo1)
# Analysis of Deviance Table (Type II tests)
# 
# Response: as.factor(cluster)
# LR Chisq Df Pr(>Chisq)    
# V3      0.653  1   0.418944    
# V4      0.064  1   0.800149    
# V5      0.670  1   0.413121    
# V6      0.468  1   0.494094    
# V7      0.009  1   0.926298    
# Lat    32.004  1  1.538e-08 ***
#   PDSI    7.775  1   0.005298 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > Anova(mo_MDS1)
# Anova Table (Type II tests)
# 
# Response: (PCo1_loading)
# Sum Sq  Df F value    Pr(>F)    
# V3        0.0169   1  0.9824 0.3222973    
# V4        0.0972   1  5.6382 0.0181120 *  
#   V5        0.2232   1 12.9508 0.0003659 ***
#   V6        0.0000   1  0.0003 0.9859697    
# V7        0.0099   1  0.5759 0.4484379    
# Lat       0.2260   1 13.1138 0.0003364 ***
#   PDSI      0.1461   1  8.4807 0.0038196 ** 
#   Residuals 6.0314 350                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > Anova(mo_MDS2)
# Anova Table (Type II tests)
# 
# Response: (PCo2_loading)
# Sum Sq  Df F value    Pr(>F)    
# V3         1.6712   1 57.9316 2.543e-13 ***
#   V4         0.9664   1 33.4999 1.584e-08 ***
#   V5         0.0547   1  1.8947   0.16955    
# V6         0.1292   1  4.4779   0.03504 *  
#   V7         0.1336   1  4.6303   0.03210 *  
#   Lat        0.0623   1  2.1586   0.14267    
# PDSI       0.0000   1  0.0000   0.99481    
# Residuals 10.0969 350                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# what I want to do is to measure the relative contributions of latitude, and other variables on the PCo1 and PCo2 loadings
# next step is to do the bray curtis distance. then do adonis anova.
# convert OTU keep to RA
otu_RA <- otu_keep/rowSums(otu_keep)
otu_RA <- otu_RA[rownames(fin2),]
otu_dist <- (vegdist(otu_RA, method = "bray"))
ad <- adonis(otu_dist~V3+V4+V5+V6+V7+Lat+PDSI,data =fin, permutations = 999, method = "bray",
       strata = NULL, contr.unordered = "contr.sum",
       contr.ordered = "contr.poly")
