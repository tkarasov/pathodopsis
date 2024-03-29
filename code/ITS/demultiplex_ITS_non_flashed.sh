#!/bin/sh
#this script demultiplexes the non_flashed reads. This needs to be used with the dada2 pipeline

#cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads
read_direc=$1
output_direc=$2


#Now iterate through the folders in the raw read directory
for direc in `ls $read_direc`; do
	name=`ls $read_direc/$direc/| grep L001_R1 | cut -f1 -d '_'`;
	R1=`ls $read_direc/$direc/*R1*`;
	R2=`ls $read_direc/$direc/*R2*`;
	qsub -v name=$name,output_direc=$output_direc,R1=$R1,R2=$R2 /ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/ITS/demultiplex_python_call_ITS.sh;
done

#once all of these files are done then call clip_barcodes_after_python.sh to remove barcodes and rename the files
/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/clip_barcodes_after_python_ITS.sh
