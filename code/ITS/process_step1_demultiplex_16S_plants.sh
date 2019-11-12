#!/bin/bash
#$ -N run_dada2
#
#$ -cwd
#
# Request nine cores
#$ -pe parallel 9
#
# Request 8Gb per process
#$ -l h_vmem=8G
#
# Merge standard error and output
#$ -j yes
#$ -o $HOME/tmp/stdout_of_job_demultiplex
#The goal of this script is to take the reads from a multiplexed 16S sequencing and demultiplex/clip them of their barcodes and rename the files

#location of scripts
conda activate pathodopsis

S16="/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S"
amp_general="/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/amplicon_general"

input1=/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/16S/run141_16S_12_2018
input2=/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/16S/run159_16S_6_2019
input3=/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/16S/run163_16S_6_2019
output_direc=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all

silva=/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database


############
mkdir $output_direc/run141_16S_12_2018
mkdir $output_direc/run159_16S_6_2019
mkdir $output_direc/run163_16S_6_2019

mkdir $output_direc/run141_16S_12_2018/demult_python
mkdir $output_direc/run159_16S_6_2019/demult_python
mkdir $output_direc/run163_16S_6_2019/demult_python
############

#Now we demultiplex the three separate runs
$S16/demultiplex_16S_non_flashed.sh $input1 $output_direc/run141_16S_12_2018

$S16/demultiplex_16S_non_flashed.sh $input2 $output_direc/run159_16S_6_2019

$S16/demultiplex_16S_non_flashed.sh $input3 $output_direc/run163_16S_6_2019

#sleep for an hour to let shell scripts finish
sleep 1h

#Once the demultiplexing is done, we need to cat the reads together from all three runs. First we will find the union of all barcodes observed across the runs.
cd $output_direc
mkdir all_runs

#list the barcodes observed per run
ls $output_direc/run141_16S_12_2018/demult_python/ > run141.txt
ls $output_direc/run159_16S_6_2019/demult_python/ > run159.txt
ls $output_direc/run163_16S_6_2019/demult_python/ > run163.txt

#Find the union of the barcodes across runs
cat run141.txt run159.txt run163.txt | sort | uniq -d > union.txt

rm run141.txt run159.txt run163.txt

#and now cat together all observed files
for file1 in `cat union.txt`; do
    # cat them together
    echo $file1
    cat $output_direc/run141_16S_12_2018/demult_python/$file1 $output_direc/run159_16S_6_2019/demult_python/$file1 $output_direc/run163_16S_6_2019/demult_python/$file1> $output_direc/all_runs/all.$file1
done

#now we need to remove the barcodes from the files and rename them. The resulting files go into demult_python

demultiplexing_file=/ebio/abt6_projects9/pathodopsis_microbiomes/data/
$S16/clip_barcodes_after_python.sh $output_direc/ $demultiplexing_file




#info=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/16S_plate_locations.csv

#read through every line of info file, and get renaming
#cd $output_direc/16S_all/all_runs/
#while IFS= read -r line;
#do
#    read -r f1 f2 f3 f4 <<<"$line"
#    plate_pos=$f3
#   barcodes=$f4
#    demult_file_R1=`ls | grep $plate_pos.fastq | grep $barcodes | grep R1`
#    demult_file_R2=`ls | grep $plate_pos | grep $barcodes | grep R2`
#    echo $demult_file_R1
#    cp $demult_file_R1 ${f1}_16S_R1.fastq
#    cp $demult_file_R2 ${f1}_16S_R2.fastq
# done<$info


#and remove all of the intermediate fastqs
rm $output_direc/all_runs/strip*
rm $output_direc/all_runs/all.*

