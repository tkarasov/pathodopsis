library(phyloseq)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(leaps)

# The goal of this script is to Perform model comparisons to look at whether PDSI is important within region. This also produces figures of the relationship between latitude and abundance.


# https://uc-r.github.io/model_selection#compare

####################
# Load phyloseq object & pathogen identities
####################
 load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_clim.rds")
 load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/OTU_clim.rds")
 load('/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000.rds')

#load("~/Documents/Google Drive/pathodopsis/temp_stuff/plant_clim.rds")
#load("~/Documents/Google Drive/pathodopsis/temp_stuff/OTU_clim.rds")
#load('~/Documents/Google Drive/pathodopsis/temp_stuff/OTUtab_GP1000.rds')

samp_data <- plant_clim$clim_data
rownames(samp_data) <- samp_data$Sequence_ID
GP_red <- phyloseq::prune_samples(rownames(samp_data), GP1000)


plant_otu <- phyloseq(sample_data(samp_data), 
                      otu_table = plant_clim$otu_table, 
                      phy_tree = plant_clim$phy_tree, 
                      refseq = plant_clim$refseq,
                      tax_table = plant_clim$tax_table)

sample_data(GP_red) <- sample_data(plant_otu)
# recalculate date:
new_date <- c(plant_clim$clim_data$Date) - 20180200
i = 0
for(rec in new_date){
  i=i+1
  hm = rec
  if(hm<100) rec1= rec
  if(hm>=100 & hm < 200) rec1= 28 +  rec - 100
  if(hm>=200 & hm < 300) rec1= 28 + 31 +  rec - 200
  if(hm>=300 & hm < 400) rec1= 28 + 31 + 30 + rec - 300
  new_date[i]=rec1
}
plant_clim$clim_data$Date <- new_date

c("PDSI", "srad", "vpd", "Elevation", "ws", "pet", "Lat", "Date" )

list_facs <- c("Lat", "Long", "Soil_temp", 
               "Soil_hum", "Humidity_ground.1", 
               "Air_temp", "Air_hum", "R_diameter")


# vars <- c("PDSI", "srad", "vpd", "Elevation", "ws", "pet", 
#          "Soil_temp", "tmax", "Air_temp", "soil", "Ath.Ath", "def", 
#          "Strata_herbs", "Air_hum", "Soil_hum", "Date")
  
vars <- c("PDSI", "Date", "Developmental_state", "R_diameter", "Herbivory", "Disease_RL", 
          "Tour_ID",  "cluster",  "Lat")

poo <- plant_clim$clim_data[,vars]
poo[poo=="n/a"]<-NA
plant_red <-data.frame(lapply(poo[,c(vars )], as.numeric))
plant_red$Sequence_ID <- plant_clim$clim_data$Sequence_ID
#plant_red <- scale(plant_red)
#plant_red[,vars] <- scale(apply(plant_red[,vars], function(x) as.numeric(as.character(x))) )

####################
# Ordination
####################
num_eig <- dim((plant_clim$otu_table))[1]-1
my.responseorig <- data.frame(sqrt(plant_clim$otu_table/1000)) %>% 
  dist() %>% cmdscale(eig = F, k = num_eig) %>% data.frame()

colnames(my.responseorig)[1:2] = c("MDS1", "MDS2")

col1 = plant_clim$clim_data$Tour_ID

#################################
# Set response variable as MDS1 or as cluster
#################################
my.response1 <- my.responseorig$MDS1 #my.responseorig[,1][match(data_frame_predictors$Sequence_ID, rownames(my.responseorig))] 

my.response2 <- make.names(as.factor(plant_clim$clim_data$cluster))

my.response3 <- my.responseorig$MDS2

#################################
# Setting up the model
#################################
plant_red <- as.data.frame(plant_red)
plant_red$MDS1 <- my.response1
plant_red$MDS2 <- my.response3
plant_red$Tour_ID <- as.factor(plant_red$Tour_ID)
 
#################################
# Model comparisons including Lat
#################################
plant_red <- plant_red %>% filter(is.na(PDSI)==FALSE)
m1 <- lm(MDS1 ~ PDSI +  Lat, data = plant_red)
m2  <- lm(MDS2 ~ PDSI  + Lat, data = plant_red)

m1_noLat <- lm(MDS1 ~ Lat, data = plant_red)
AIC(m1, m1_noLat)

anova(m1, m1_noLat, test = "F")
#

mc1 <- glm(as.factor(cluster) ~ PDSI + 
               Lat, data = plant_red, family = "binomial")

mc1_noP <- glm(as.factor(cluster) ~  
             Lat, data = plant_red, family = "binomial")

AIC(mc1, mc1_noP)

anova(mc1, mc1_noP, test = "Chisq")

#################################
# Model comparisons including Tour ID
#################################

m1 <- lmer(MDS1 ~ PDSI +  1|Tour_ID, data = plant_red)
m2  <- lmer(MDS2 ~ PDSI  + 1|Tour_ID, data = plant_red)

m1_noP <-  lmer(MDS1 ~  1|Tour_ID, data = plant_red)
m2_noP <-  lmer(MDS2 ~  1|Tour_ID, data = plant_red)

check_model(m2)

anova(m1_noP, m1)
anova(m2_noP, m2)
anova(m1, m1_noP, test = "F")

mc1 <- glmer(as.factor(cluster) ~ PDSI +  1|Tour_ID, data = plant_red, family = "binomial")
mc2  <- glmer(as.factor(cluster) ~ PDSI  + 1|Tour_ID, data = plant_red, family = "binomial")

mc1_noP <-  glmer(as.factor(cluster) ~  1|Tour_ID, data = plant_red, family = "binomial")
mc2_noP <-  lmer(MDS2 ~  1|Tour_ID, data = plant_red)

anova(mc1_noP, mc1)

#################################
# Which microbes associated with PDSI
#################################
lm_rec<- data.frame(otu = colnames(plant_clim$otu_table), pval = c(NA)*575)

for(i in 1:dim(plant_clim$otu_table)[2]){
  df = data.frame(OTU = plant_clim$otu_table[,i], Lat =  plant_clim$clim_data$Lat ) 
  colnames(df) = c("OTU", "Lat")
  model <- lm( OTU~Lat, data = df)
  lm_rec[i,2] = summary(model)$coefficients[2,4]
  lm_rec$freq[i] = sum(plant_clim$otu_table[,i])/(1000*dim(plant_clim$otu_table)[1])
}

lm_rec$FDR <- p.adjust(lm_rec$pval, method = "BH")

lm_rec$sig_0.01 <- lm_rec$FDR<0.01

lm_rec <- cbind(lm_rec, plant_clim$tax_table )

lm_fam <- lm_rec %>% filter(Family %in% names(which(table(lm_rec$Family)>=5)))

# Most significant associations with Lat come from Sphingomonadaceae

pval_plot <- ggplot(aes(x = Family, y = -log10(FDR), col = freq*100), data = lm_fam, col = freq) + 
  geom_point(aes(size = freq))+ scale_color_viridis_c( direction = -1 ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 2, col = "Red", lty = "dashed", alpha = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
theme(legend.position = "top")
  

pval_plot <- pval_plot + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)




sphingo <- data.frame(plant_clim$otu_table)
sphingo <- sphingo[,which(plant_clim$tax_table[,5] == "Sphingomonadaceae")]
sphingo$Lat <- plant_clim$clim_data$Lat

sph_gather <- sphingo %>% gather(key = "key", value = "value", -Lat )


sp_plot <- ggplot(data = sph_gather, aes(x = Lat, y = value/10, col = key)) + 
 # geom_point()
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  ylim(c(0,15)) +
  xlab("Latitude") +
  ylab("RA (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None") + 
  scale_color_brewer(palette = "Set1") +
  annotate("text",  x=60, y = 15, label = "Sphingomonadaceae", vjust=1, hjust=1)
  



pseudo <- data.frame(plant_clim$otu_table)
pseudo <- pseudo[,which(plant_clim$tax_table[,5] == "Pseudomonadaceae")]
pseudo$Lat <- plant_clim$clim_data$Lat

ps_gather <- pseudo %>% gather(key = "key", value = "value", -Lat )


ps_plot <- ggplot(data = ps_gather, aes(x = Lat, y = value/10, col = key)) + 
  # geom_point()
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  ylim(c(0,1)) +
  xlab("Latitude") +
  ylab("RA (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "None" ) + 
  scale_color_brewer(palette = "Set1") +
  annotate("text",  x=60, y = 1, label = "Pseudomonadaceae", vjust=1, hjust=1)


p2 <- plot_grid(nrow = 2, sp_plot, ps_plot, align = "hv")


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/pval_lat_asv.pdf", useDingbats = FALSE, family = "ArialMT", width  = 7, height=4)
lat = plot_grid(pval_plot + theme(axis.text.x=element_text(size=rel(0.5))), p2, ncol = 2, labels = c("A", "B", "C"))
lat
dev.off()

saveRDS(lat, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/lat_fdr.rds")

load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/kriging_pseud_sping.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/Fst_plot.rds")

p1 <- plot_grid(pval_plot + theme(axis.text.x=element_text(size=rel(0.5))), p2, ncol = 2, rel_widths = c(2,1))
p3 <- plot_grid(pseud_sphing_grid, Fst_plot)

pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/first_part_lat.pdf", 
    useDingbats = FALSE, family = "ArialMT", width  = 7, height = 4)
p1
dev.off()
# 
# forward <- regsubsets(MDS1 ~ ., plant_red, nvmax = 19, method = "forward")
# backward <- regsubsets(MDS1 ~ ., plant_red, nvmax = 19, method = "backward")
# 
# 
# # create training - testing data
# set.seed(1)
# sample <- sample(c(TRUE, FALSE), nrow(plant_red), replace = T, prob = c(0.6,0.4))
# train <- plant_red[sample, ]
# test <- plant_red[!sample, ]
# 
# # perform best subset selection
# best_subset <- regsubsets(MDS1 ~ ., plant_red, nvmax = 12)
# results <- summary(best_subset)
# 
# # extract and plot results
# tibble(predictors = 1:12,
#        adj_R2 = results$adjr2,
#        Cp = results$cp,
#        BIC = results$bic) %>%
#   gather(statistic, value, -predictors) %>%
#   ggplot(aes(predictors, value, color = statistic)) +
#   geom_line(show.legend = F) +
#   geom_point(show.legend = F) +
#   facet_wrap(~ statistic, scales = "free")
# 
# #################################
# # Directly Esimating Test error
# #################################
# 
# test_m <- model.matrix(MDS1 ~ ., data = plant_red)
# # create empty vector to fill with error values
# validation_errors <- vector("double", length = 19)
# 
# for(i in 1:12) {
#   coef_x <- coef(best_subset, id = i)                     # extract coefficients for model size i
#   pred_x <- test_m[ , names(coef_x)] %*% coef_x           # predict salary using matrix algebra
#   validation_errors[i] <- mean((test$Salary - pred_x)^2)  # compute test error btwn actual & predicted salary
# }
# 
# # plot validation errors
# plot(validation_errors, type = "b")
# 
# 
# #cross validation
# predict.regsubsets <- function(object, newdata, id ,...) {
#   form <- as.formula(object$call[[2]]) 
#   mat <- model.matrix(form, newdata)
#   coefi <- coef(object, id = id)
#   xvars <- names(coefi)
#   mat[, xvars] %*% coefi
# }
# 
# k <- 10
# set.seed(1)
# folds <- sample(1:k, nrow(hitters), replace = TRUE)
# cv_errors <- matrix(NA, k, 19, dimnames = list(NULL, paste(1:19)))
# 
# for(j in 1:k) {
#   
#   # perform best subset on rows not equal to j
#   best_subset <- regsubsets(Salary ~ ., hitters[folds != j, ], nvmax = 19)
#   
#   # perform cross-validation
#   for( i in 1:19) {
#     pred_x <- predict.regsubsets(best_subset, hitters[folds == j, ], id = i)
#     cv_errors[j, i] <- mean((hitters$Salary[folds == j] - pred_x)^2)
#   }
# }
# 
# 
# mean_cv_errors <- colMeans(cv_errors)
# 
# plot(mean_cv_errors, type = "b")
# 
# final_best <- regsubsets(Salary ~ ., data = hitters , nvmax = 19)
# coef(final_best, 11)
