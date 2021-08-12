
library(phyloseq)
library(dplyr)
library(reshape2)
devtools::install_git(url = 'https://github.com/tkarasov/taliaRgeneral.git')
library(taliaRgeneral)
library(parallel)
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(viridis)
#githubURL <- 
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/color_pal.rds")
load("/ebio/abt6_projects9/pathodopsis_microbiomes/taliaRgeneral/R/theme_pubs.rds")
#download.file(githubURL, c("color_vector.rds","theme_pubs.rds"), method="curl")

# The goal of this script is to take a phyloseq object and to compare abundance between soil and capsella and A. thaliana
#########################################
# Functions
#########################################

per_site <- function(GP, seq, comp1="between", comp2="capsella", comp3="soil"){
  #Compiles data for comparisons between soil and species
  # seq must be a vector
  if(comp1 == "between"){
    ath1 = subset_samples(GP, Host_Species == "Ath")
    ath1 = prune_taxa(c(seq), ath1)
    ath1 = make_replicate(ath1)
    fin1 = cbind(sample_data(ath1), otu_table(ath1))
    #hm = fin1$Host_Species == "Soil"
    #fin1$Used_for == hm
    fin1 = reshape2::dcast(fin1, Site_ID ~ rep, value.var = eval(seq))
  }
  if(comp2 == "capsella"){
    ath2 = subset_samples(GP, Host_Species == "Cap")
    ath2 = prune_taxa(c(seq), ath2)
    fin2 = cbind(sample_data(ath2), otu_table(ath2))
    hm = fin2$Host_Species == "Cap"
    fin2$Used_for = hm
    fin2 = reshape2::dcast(fin2, Site_ID ~ Host_Species, mean, value.value = eval(seq))
  }
  if(comp3 == "soil"){
    ath3 = subset_samples(GP, Host_Species == "Soil")
    ath3 = prune_taxa(c(seq), ath3)
    fin3 = cbind(sample_data(ath3), otu_table(ath3))
    hm = fin3$Host_Species == "Soil"
    fin3$Used_for = hm
    fin3 = reshape2::dcast(fin3, Site_ID ~ Used_for, mean, value.var = eval(seq))
  }
  
  #merge all
  all_comp <- left_join(fin1, fin2)
  all_comp <- left_join(all_comp, fin3)
  
  return(data.frame(all_comp))
}


make_replicate <- function(GP){
  #adds sample data column
  unique_loc = unique(sample_data(GP)$Site_ID)
  sample_data(GP)$rep = c(0)
  for(val in unique_loc){
    chosen = which(sample_data(GP)$Site_ID==val)
    sample_data(GP)$rep[chosen] = paste("rep", c(1:length(chosen)), sep="")
  }
  return(GP)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(#a = format(unname(coef(m)[1]), digits = 2),
                        #b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

run_lm<-function(per_site_object, variable1, variable2, colname=""){
  per_site_object = data.frame(per_site_object)
  #remove NAs
  #per_site_object = per_site_object %>% filter(!is.na(Cap)) %>% filter(!is.na(SOIL))
  per_site_object$mean_rep <- apply(per_site_object, 1, function(x) mean(c(x["rep1"],x["rep2"]), na.rm = TRUE))
  equation = paste(
    paste("log10(", 
          paste(eval(variable1),"+1)~",sep=""), 
          paste("log10(", eval(variable2), sep=""), sep=""),"+1)", sep="")
  my_lm <- summary(lm(data=per_site_object, 
              as.formula(equation)))
  to_return <-c(my_lm$r.squared, my_lm$coefficients[,4][c(2)])
  names(to_return) = paste( c("R2", "pval"),colname, sep="_")
  
  return(c(my_lm$r.squared, my_lm$coefficients[,4][c(2)]))
}

run_log<-function(per_site_object){
  #print("yes")
  per_site_object = data.frame(per_site_object)
  per_site_object[per_site_object$rep1 > 0,]$rep1 <-1
  per_site_object[which(per_site_object$rep2 > 0),]$rep2 <-1
  per_site_object$combined_ath <- rowSums(per_site_object[,c("rep1", "rep2")], na.rm=TRUE)
  per_site_object[per_site_object$combined_ath>0,]$combined_ath=1
  #per_site_object[which(per_site_object$Cap > 0),]$Cap<-1
  if(var(per_site_object$combined_ath, na.rm = TRUE)==0)
    return("NA")
  if(var(per_site_object$Cap, na.rm = TRUE) == 0)
    return("NA")
  else{
    glm.fit <- summary(glm(combined_ath ~ Cap, 
                 data = per_site_object, family = binomial))
  return(glm.fit$coefficients[2,4])
  }
}

#########################################
# Pathogen labels
#########################################
OTU5 <- c("seq_10")
Xan <- "seq_49"
sphing <- "seq_1"


#########################################
# Load ITS and 16S data
#########################################
load("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/OTUtab_GP1000_at15.rds")

#########################################
# Build data frame that has Capsella, Athalian and Soil per ASV
#########################################

tax_names <- colnames(otu_table(GP_at15_all))

ath_all <-mcsapply(tax_names, function(x) per_site(GP_at15_all, x), mc.cores = 16)





#########################################
# Linear regression to determine whether surrounding species is a good predictor
#########################################
ath_lm <- data.frame(t(apply((ath_all), 2, 
                             function(x) run_lm(x, "rep1", "rep2", colname="ath" ))))
ath_cap <- data.frame(t(apply((ath_all), 2, 
                              function(x) run_lm(x, "rep1", "Cap", colname="capsella" ))) )
ath_soil <- data.frame(t(apply((ath_all), 2, 
                               function(x) run_lm(x, "rep1", "TRUE", colname="soil" ))))
ath_lm_all <- cbind(ath_lm, ath_cap, ath_soil)

colnames(ath_lm_all) <- c("pval_rep", "R2_rep", "pval_cap", "R2_cap", "pval_soil", "R2_soil")

#########################################
# Logistic regression for relationship with other species
#########################################
ath_log <- data.frame(t(apply((ath_all), 2, 
                             function(x) run_log(x))))


#########################################
# Tree
#########################################
# https://yulab-smu.github.io/treedata-book/chapter7.html
tree <- phy_tree(GP_at15_all)
circ <-ggtree(phy_tree(GP_at15_all), layout = "circular", branch.length = "none")
rownames(ath_lm_all) = tree$tip.label
df <- ath_lm_all[,c(1,3,5)]

df_fam <-data.frame(tax_table(GP_at15_all)[,3])
rownames(df_fam) = tree$tip.label

set.seed(24)
p1 <- gheatmap(circ, df_fam, offset = 0.2, width=0.2, colnames = FALSE,
               color = NULL) +
  scale_fill_manual(values=sample(col_vector, 15))
  #scale_fill_viridis(discrete = TRUE)

p2 <- p1 + new_scale_fill()
  
p3 <- gheatmap(p2, df, offset=10, width=0.5,
               colnames_angle=95, colnames_offset_y = .25,
               color = NULL,
               font.size = 7) +
  scale_fill_viridis(option="A", name="continuous\nvalue")  

colnames(ath_lm) <- c("pval_capsella", "pval_soil")
ath_lm_cap <- ath_lm#filter(ath_lm, !is.na(pval_capsella))
ath_lm_cap$fdr_cap <- p.adjust(ath_lm_cap$pval_capsella)


pdf("/ebio/abt6_projects9/pathodopsis_microbiomes/data/figures_misc/circos_r2_ASVs.pdf", useDingbats = FALSE, family = "ArialMT")
p3
dev.off()


p1 <- ggplot(data = ath, aes(x = (rep1/10 +1) , y = (rep2/10 +1))) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()
  

#########################################
# Correlation between A. thaliana and Capsella
#########################################
ath_cap <- t(data.frame(ath_all))
cap <- data.frame(ath_cap[,5])
ath <- data.frame(ath_cap[,6])
for(i in c(1:10)){
  temp = data.frame(cbind(ath[,i], cap[,i]))
  colnames(temp) = c("athal", "cap")
  m = lm(log10(temp$cap+1) ~ log10(temp$ath + 1))
  r2 = format(summary(m)$r.squared, digits = 2)
  pval = format(summary(m)$coefficients[2,4], digits = 2)
  temp2 = ggplot(data = temp, aes(x = (athal + 1)/10, y = (cap + 1)/10)) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
                  geom_point(cex = 2, alpha = 0.5) +
                   scale_x_log10()+
                   scale_y_log10()+
                  theme_pubs 
    #geom_text(x = 25, y = 300, label = paste(r2, pval, sep=" "), parse = TRUE)
  assign(paste("plot", i, sep=""), temp2)
}

full_plot <- plot_grid(plot1, plot2, plot3, 
                       plot4, plot5, plot6, 
                       plot7, plot8, plot9)
#########################################
# Correlation between A. thaliana and Capsella
#########################################
