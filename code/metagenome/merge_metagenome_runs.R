#!/usr/bin/env Rscript
#the only purpose of this script is to take the different metagenome runs and merge the data. I had some issues with the capsella merging and one of the controls merging so have excluded them here in 4/2020

library(dplyr)
library(tibble)

setwd("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs")
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


merge_microbe_reads_by_name <- function(){
  #list all samples across all folders. Names should end in R1R2.fq.gz
  read_list = c()
  for(fold in folders){
    my_files = list.files(path = fold, pattern = "R1.fq")
    my_files = sapply(strsplit(my_files, "_"), "[", 1)
    my_files = sapply(strsplit(my_files, "Metagenomic"), "[", 1)
    my_files2 = gsub("R1", "R2", my_files)
    my_files = gsub(".R1", "*.R1",my_files)
    my_files2 = gsub(".R2", "*.R2", my_files2)
    print(my_files[1])
    read_list = append(read_list, my_files)
    read_list = append(read_list, my_files2)
  }
  for(read in unique(read_list)){
    for(fold in folders){
      samp=paste(fold, read, sep="")
      samp_name=gsub("\\*", "", read)
      dest = paste(paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/", samp_name, sep=""), "", sep="")
      cmd = paste(paste(paste("cat", samp, sep= " "), ">>", sep=""),dest, sep="")
      com = system(cmd)
    }
  }
  
}

merge_plant_reads_by_name <- function(){
  #list all samples across all folders. Names should end in R1R2.fq.gz
  read_list = c()
  for(fold in folders){
    my_files = list.files(path = fold, pattern = ".bam")
    my_files2 = grep(my_files, pattern = "unmapped", inv = T, value = T)
    my_files = sapply(strsplit(my_files2, "_"), "[", 1)
    my_files = gsub(".bam", "", my_files)
    #my_files = sapply(strsplit(my_files, "Metagenomic"), "[", 1)
    #my_files = sapply(strsplit(my_files, ".fq.gz"), "[",1)
    print(my_files[1])
    read_list = append(read_list, my_files)
  }
  read_list = unique(read_list)
  for(read in read_list){
    cmd = paste("touch",
              paste(
                paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/", read, sep=""), 
                ".bam", sep = ""), sep = " ")
    com = system(cmd)
      }

  
  for(read in read_list[2:length(read_list)]){
    tot_files = c()
    dest =paste(paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/", 
                      read, sep=""), ".bam", sep = "")
    for(fold in folders){
      samp = list.files(path=fold, pattern = read )
      samp = grep(samp, pattern = ".bam", value = T)
      samp = grep(samp, pattern = "unmapped", inv = T, value = T)
      my_files2 = grep(my_files, pattern = "unmapped", inv = T, value = T)
      #print(my_files2)
      if(length(samp)==0)
        print("NA")
      else{
        samp_name=paste(fold, samp, sep="")
        #samp_name=gsub("\\*", "", read)
        #dest = paste(paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/", samp_name, sep=""), "", sep="")
        tot_files = append(tot_files, samp_name)
      }
        print(tot_files)
        tot_char = stri_paste(tot_files, collapse=' ')
        cmd = paste(paste("samtools merge ", dest, sep = " " ), tot_char, sep = " ")
        print(cmd)
        com = system(cmd)
        #cmd2 = paste("mv temp ", dest, sep = " ")
        #com2 = system(cmd2)
      }
    }
  }
  



merge_original_reads <- function(){
  dir="/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/metagenome"
  file_list <-  list.dirs(dir)
  samp <- sapply(strsplit(sapply(strsplit(file_list, "SampleId"), "[",2), "_"), "[",1)
  samp <- unique(samp)
  
  #Exclude the capsella samples and the control sample that has size problems
  samp <- grep("PC", samp, inv = TRUE, value = TRUE)
  samp <- grep("control.", samp, inv = TRUE, value =TRUE)
  
  for(sample in samp){
    mv_file_list <-  grep(sample, file_list, value = TRUE)
    mv_file_list1 <- paste(mv_file_list, "/*R1*", sep="")
    mv_file_list2 <- paste(mv_file_list, "/*R2*", sep="")
    #mv_file_list2 <- grep("R2", grep(sample, file_list, value = TRUE))
    mv_file_list1 <- stri_paste(mv_file_list1, collapse = " ")
    mv_file_list2 <- stri_paste(mv_file_list2, collapse = " ")
    new_name1 <- paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/raw_concat/", 
                       paste(sample, ".R1.fq.gz", sep = ""), sep = "")
    new_name2 <- paste("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/raw_concat/", 
                       paste(sample, ".R2.fq.gz", sep = ""), sep = "")
    cmd1 = paste(paste(paste("zcat ", mv_file_list1, sep = " "), ">>", sep = ""), new_name1, sep = "")
    cmd2 = paste(paste(paste("zcat ", mv_file_list2, sep = " "), ">>", sep = ""), new_name2, sep = "")
    print(cmd1)
    com1 <- system(cmd1)
    com2 <- system(cmd2)
  }
  
  

  }





folders = c("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run116_2018_9_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run104_2018_7_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run119_2018_9_metagenome_reads/",
            "/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run144_2019_5_metagenome_reads/"
            )



# This command writes concatenated text to the combined_runs file
#merge_reads_by_name()
#merge_plant_reads_by_name()

merge_original_reads()




