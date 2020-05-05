#!/bin/bash
#
# Use bash as default shell so that I can load conda
#$ -S /bin/bash
#
#  Reserve 6 threads
#$ -pe parallel 6
#
#  Request 16G of RAM per core
#$ -l h_vmem=32G
#
#$ -N freebayes
#
#
#$ -S /bin/bash



#this script maps all metagenome reads to acd6 gene
cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/mapping
reference="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/TAIR10_chr_all.fas"
raw_reads="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/raw_concat"

############################################################
# Call SNPs with freebayes
############################################################
#make sure the reference index is gone
rm /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_TAIR10.fasta.fai

#now run freebayes
/ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/pathodopsis/bin/freebayes -f $reference --gvcf  *.bam > pathodopsis_freebayes.vcf
#sed 's/.R1..bam//g' acd6.vcf  | sed 's/.R1.//g'> acd6_parsed.vcf

#Pull ACD6 region from vcf

#Make phylogeny?
# Take only the region around acd6
# vcftools --vcf poolsGVCF.filtered_snps_final.PASS.bi.vcf \
# 	--chr chr4 --from-bp 8294164 --to-bp 8299195 --recode \
#        	--out acd6

# vcftools --gzvcf acd6.recode.vcf.gz --freq --out allel_freq
# #output tsv (will be used for list of samples below)
# vk vcf2tsv long acd6.recode.vcf.gz > acd6.tsv

# #convert into fasta with N's       	
# vk phylo fasta acd6.recode.vcf.gz > acd6.samples.fasta

