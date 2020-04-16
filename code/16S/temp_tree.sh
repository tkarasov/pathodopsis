#!/bin/sh
#
#  Reserve 16 CPUs for this job
#$ -pe parallel 16
#  Request 64G of RAM per core
#$ -l h_vmem=64G
#$ -o $HOME/tmp/stdout_of_job_dada2
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N tree
#  Run job from current working directory
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#

#location of scripts
conda activate pathodopsis

S16="/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S"
amp_general="/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/amplicon_general"

#output_direc=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all

output_direc=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

silva=/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database

#mkdir /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

#cp /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/* /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

#cp /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/soil_16S_ITS/run174_16S_9_2019_soil_16S/demult_python/* /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

cd $output_direc
echo "Making tree"
# Create phylogeny from table
$amp_general/after_dada2_make_tree.R
