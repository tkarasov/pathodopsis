#the goal of this script is to take the full metagenome table and subset to only those htat have the highest variance
library("rwunderground")

#read in info for all plants
#plant_info=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/v1_22_5_merged.txt", sep=",", header=T)
#plant_info=subset(plant_info, select=-c(Notes))
plant_info=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course2.txt", sep="\t", header=T)

#for all latitude and longitudes, build info on temperature and humidity
#history_daily(location, date = "20150101", use_metric = FALSE, key = get_api_key(), raw = FALSE, message = TRUE)

#read in metatable
meta=read.table("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/2018_9_metagenome_reads/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)
temp=gsub("..centrifuge_output.","", colnames(meta))
temp=gsub(".R1.fq.report","",temp)
colnames(meta)=temp


#full data
meta_t_all=as.data.frame(t(meta))
sorted_meta = names(sort(apply(meta_t_all, 2 ,var, na.rm=T), decreasing = T)[1:25])
meta_t_all$Plant_ID=rownames(meta_t_all)
meta_fin_all=merge(plant_info, meta_t_all)
write.table(meta_fin_all,"/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/all_metagenome_metadata_9_2018_reads.txt", quote = F, row.names = F, sep="\t")

#top 25 of variance
meta_red=meta[sorted_meta,]
meta_t=as.data.frame(t(meta_red))
meta_t$Plant_ID=rownames(meta_t)
meta_fin=merge(plant_info, meta_t)
write.table(meta_fin,"/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/metagenome_metadata_9_2018_reads.txt", quote = F, row.names = F, sep="\t")
