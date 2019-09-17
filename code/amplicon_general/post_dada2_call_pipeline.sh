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
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
# -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
# Run job from current working directory
#$ -cwd

source /ebio/abt6/tkarasov/.bashrc
conda activate pathodopsis

