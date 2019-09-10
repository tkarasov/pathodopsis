#!/bin/sh
#this script demultiplexes the non_flashed reads. This needs to be used with the dada2 pipeline

#cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads
read_direc=$1  #/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/16S/16S_12_2018
output_direc=$2

#Now iterate through the folders in the raw read directory
for direc in `ls $read_direc`; do
	name=`ls $read_direc/$direc/| grep L001_R1 | cut -f1 -d '_'`;
	R1=`ls $read_direc/$direc/*R1*`;
	R2=`ls $read_direc/$direc/*R2*`;
	qsub -v name=$name,output_direc=$output_direc,R1=$R1,R2=$R2 /ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/demultiplex_python_call.sh;
done

