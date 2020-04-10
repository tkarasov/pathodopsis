#!/bin/sh
#
#  Reserve 4 CPUs for this job
#$ -pe parallel 4
#  Request 128G of RAM total (32per)
#$ -l h_vmem=64G
#$ -o $HOME/tmp/ITS_stdout_of_job_dada2
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N ITS_dada2
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

output_direc=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/all_runs/demult_python
#/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/

#tax_ref=/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/sh_general_release_04.02.2020
tax_ref=/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database

#mkdir /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

#cp /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python/* /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

#cp /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/soil_16S_ITS/run174_16S_9_2019_soil_16S/demult_python/* /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S_soil_phyllo_fin

cd $output_direc

# # first filter with dada2
# $amp_general/dada2_filter_tlk.R \
#        -f $output_direc --truncLen 270,200 \
#        --maxEE 3,7 --truncQ 2 --threads 9 --f_match _R1.*fastq \
#        --r_match _R2.*fastq --verbose TRUE
# 
# # then infer
# $amp_general/dada2_inference_tlk.R \
#        -f $output_direc/filtered_fastqs --seed 1659 -t 9 --verbose --plot_errors

# Remove chimeric variants and assign taxonomy
$amp_general/dada2_chimera_taxa.R -i $output_direc/seqtab.rds \
 -r $tax_ref/fin_oomycete_ITS.fasta  \
 -t 9 --skip_species --verbose TRUE \
 --tax_out $output_direc/tax_final.rds

echo "Done with taxonomy assignments, moving onto converting dada table"

# generate table
 $amp_general/convert_dada2_out.R -i \
       $output_direc/seqtab_final.rds -b $output_direc/seqtab.biom.tsv -f $output_direc/seqtab.fasta -S-taxa_in $output_direc/tax_final.rds

echo "Done with converting dada2 moving on to biom conversion"

# Create phylogeny from table. This function takes a lot of memory also
# $amp_general/after_dada2_make_tree_ITS.R


# convert table to biom
biom convert -i $output_direc/filtered_fastqs/seqtab.biom.tsv -o $output_direc/filtered_fastqs/seqtab.biom --to-hdf5

# add metadata to biom
biom add-metadata -i seqtab.biom -o seqtab_tax.biom --observation-metadata-fp  taxa_metadata.txt \
       --observation-header OTUID,taxonomy --sc-separated taxonomy

#let's look at the table for how to rarefy
biom summarize-table -i seqtab.biom -o seqtab_summary.txt

#if I want to rarefy???
# single_rarefaction.py -i seqtab_tax.biom -o seqtab_tax_rarified.biom -d 4000





















# 
# 
# 
# # first filter with dada2
#  # $amp_general/dada2_filter_tlk.R \
#  #       -f /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_all/all_runs/demult_python --truncLen 270,200 \
#  #       --maxEE 3,7 --truncQ 2 --threads 9 --f_match _R1.*fastq \
#  #       --r_match _R2.*fastq
# 
# # then infer
#  # $amp_general/dada2_inference_tlk.R \
#  #       -f $output_direc/all_runs/demult_python/filtered_fastqs --seed 1659 -qt 9 --verbose --plot_errors
# 
# 
# # Remove chimeric variants and assign taxonomy
# $amp_general/dada2_chimera_taxa.R -i $output_direc/seqtab.rds \
#  -r $silva/silva_nr_v132_train_set.fa.gz \
#  -t 9 --skip_species --verbose TRUE \
#  --tax_out $output_direc/tax_final.rds
# 
# echo "Done with taxonomy assignments, moving onto converting dada table"
# 
# # generate table
#  $amp_general/convert_dada2_out.R -i \
#        $output_direc/seqtab_final.rds -b $output_direc/seqtab.biom.tsv -f $output_direc/seqtab.fasta --taxa_in $output_direc/tax_final.rds
# 
# echo "Done with converting dada2 moving on to biom conversion"
# 
# # convert table to biom
# biom convert -i $output_direc/all_runs/demult_python/filtered_fastqs/seqtab.biom.tsv -o $output_direc/all_runs/demult_python/filtered_fastqs/seqtab.biom --to-hdf5
# 
# # add metadata to biom
# biom add-metadata -i seqtab.biom -o seqtab_tax.biom --observation-metadata-fp  taxa_metadata.txt \
#        --observation-header OTUID,taxonomy --sc-separated taxonomy
# 
# #let's look at the table for how to rarefy
# biom summarize-table -i seqtab.biom -o seqtab_summary.txt
# 
# #if I want to rarefy???
# # single_rarefaction.py -i seqtab_tax.biom -o seqtab_tax_rarified.biom -d 4000
