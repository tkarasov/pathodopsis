#library("ClimClass")
library('kgc')

data = read.table("~/work_main/abt6_projects9/pathodopsis_microbiomes/data/meta_data/Pathodopsis_site_metadata_20190808.txt", sep="\t", fill = TRUE, header = T)

data <- data.frame(data,
                   rndCoord.lon = RoundCoordinates(data$Long),
                   rndCoord.lat = RoundCoordinates(data$Lat))

data <- data.frame(data,ClimateZ=LookupCZ(data, res="course"))
#look up climate zone

write.table(data, "~/work_main/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/data/Pathodopsis_site_metadata_20190808_CZ_course.txt", quote = FALSE, row.names = FALSE, sep="\t")
