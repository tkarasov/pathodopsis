#!/bin/sh
#
#  Reserve 1 CPUs for this job
#  Request 10G of RAM
#$ -l h_vmem=10G
#
# The name shown in the qstat output and in the output file(s). The
# default is to use the script name.
# -N $name
#
# The path used for the standard output stream of the job
# -o
#
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
# -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
# Run job from current working directory
#$ -cwd

name=$name
R1=$R1
R2=$R2


python /ebio/abt6_projects9/pathodopsis_microbiomes/scripts/ITS/demultiplex_ITS_non_flashed.py $name $R1 $R2


