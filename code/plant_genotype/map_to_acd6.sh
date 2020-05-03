#!/bin/bash

#this script maps all metagenome reads to acd6 gene
cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/mapping
reference="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_TAIR10.fasta"
raw_reads="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/raw_concat"

bwa index ${reference}
total_files=`ls $raw_reads | grep R1.fq.gz | wc -l`
arr=( $(ls $raw_reads | grep R1.fq.gz) )
echo "mapping started" >> map.log
echo "---------------" >> map.log



############################################################
# Map reads with bwa
############################################################

for ((i=0; i<$total_files; i+=1))
{
sample_name1=`echo ${arr[$i]} | awk -F "_" '{print $1}'`
sample_name2=`echo ${arr[$i]} | sed 's/R1/R2/g' | awk -F "_" '{print $1}'`
echo "[mapping running for] $sample_name1"
printf "\n"
echo "bwa mem -t 12 ${reference} ${arr[$i]} ${arr[$i+1]} > $sample_name1.sam" >> map.
read1=$raw_reads/$sample_name1
read2=$raw_reads/$sample_name2
sample=`echo $sample_name1 | sed 's/fq.gz//g'`
bwa mem -R "@RG\tID:$sample\tSM:$sample" -t 12 ${reference} ${read1} ${read2}  |\
    /ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/samtools/samtools view -b -F 4 |
    /ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/samtools/samtools sort -m 970m - > "$sample".bam 
}

############################################################
# Call SNPs with freebayes
############################################################
#make sure the reference index is gone
rm /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_TAIR10.fasta.fai

#now run freebayes
freebayes -f $reference --gvcf  *..bam > acd6.vcf