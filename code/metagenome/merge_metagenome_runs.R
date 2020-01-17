#!/usr/bin/env Rscript
#the only purpose of this script is to take the output from the centrifuge runs on different datasets and merge them

library(dplyr)
library(tibble)

combine_tables <- function(df1, df2){
  #combine coverage per genome together
  bind_rows(df1 %>% tibble::rownames_to_column(), 
            df2 %>% tibble::rownames_to_column()) %>% 
    # evaluate following calls for each value in the rowname column
    group_by(rowname) %>% 
    # add all non-grouping variables
    summarise_all(sum) -> combined_df
  return(combined_df)
}


merge_by_name <- function(name){
  #merge data from all relevant folders. Just need to provide name of desired taxonomic group
  df1 = read.csv(paste(folders[1], name, sep = ""), header = T, row.names = 1)
  
  for(folder in folders[1:length(folders)[1]]){
    df2 = read.csv(paste(folder, name, sep = ""), header = T, row.names = 1)
    name_join = combine_tables(df1, df2)
    df1 = data.frame(name_join, row.names = 1)
  }
  return(df1/length(folders)) # need to get the average coverage per genome
}



folders = c("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run116_2018_9_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run104_2018_7_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run119_2018_9_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run144_2019_5_metagenome_reads/"
            )

fungi = "meta_family_corrected_per_plant_v2_fungi.csv"
bac = "meta_family_corrected_per_plant_v2_fungi.csv"
oom = "meta_family_corrected_per_plant_v2_fungi.csv"

# Merge the tables
fung_all = merge_by_name(fungi)
bac_all = merge_by_name(bac)
oom_all = merge_by_name(oom)

write.table(fung_all, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/fungi_all_runs_meta_family_corrected_per_plant_v2.csv", sep = ",", header = T, quote = F)

write.table(bac_all, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/bac_all_runs_meta_family_corrected_per_plant_v2.csv", sep = ",", header = T, quote = F)

write.table(oom_all, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/oom_all_runs_meta_family_corrected_per_plant_v2.csv", sep = ",", header = T, quote = F)

#write out summed load file

fung_sum = colSums(fung_all)
bac_sum = colSums(bac_all)
oom_sum = colSums(oom_all)

all_load = data.frame(fung_load = colSums(fung_all))
all_load$bac_load = bac_sum[rownames(all_load)]
all_load$oom_load = oom_sum[rownames(all_load)]
all_load$Total_load = rowSums(all_load)

write.table(fung_all, "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/Total_load_v2.csv", sep = ",", header = T, quote = F)






# folder = "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run116_2018_9_metagenome_reads/"
# fung = read.csv(paste(folder, fungi, sep = ""), header = T, row.names = 1)
# 
# fung_main = data.frame()
# for(folder in folders){
#   fung_temp = read.csv(paste(folder, fungi, sep = ""), header = T, row.names = 1)
#   fung_join = join(fung_main, fung_temp)
#   fung_main = fung_join
#   
# }





