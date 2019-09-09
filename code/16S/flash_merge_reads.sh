#!/bin/sh
#
#  Reserve 1 CPUs for this job
#  Request 10G of RAM
#$ -l h_vmem=10G
#
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
# -N $name
#
#  The path used for the standard output stream of the job
# -o
#
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
# -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#  Run job from current working directory
#$ -cwd
#

cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads

R1=$R1
R2=$R2
merged=$merged

echo "R1 is " $R1
echo "R2 is " $R2
echo "merged is " $merged



flash=/ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/bin/flash

#the goal of this script is to merge 16S reads together using FLASH
flash $R1 $R2 -m 30 -M 300 -o "$merged" -r 300 -f 310
