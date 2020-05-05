#!/bin/bash
#
#  Reserve 16 CPUs for this job
#$ -pe parallel 12
#
#  Request 16G of RAM per core
#$ -l h_vmem=16G
#
#$ -N acd6_map
#
# The following sets the task as an array job which range from 1 to 1179 with a step size of 1, where it will submit total of 4 jobs, one at a time
#$ -t 1-589:1
#
#$ -S /bin/bash



#this script maps all metagenome reads to acd6 gene
cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/mapping
reference="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/acd6_mapping/todesco_major_alleles.fasta"
raw_reads="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/combine_runs/raw_concat"

# The reference needs to be indexed but if indexed already don't run
#bwa index ${reference}
total_files=`ls $raw_reads | grep R1.fq.gz | wc -l`
echo "mapping started" >> map.log
echo "---------------" >> map.log



############################################################
# Map reads with bwa
############################################################
arr=($(ls $raw_reads | grep R1.fq.gz))
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"




#for ((i=0; i<$total_files; i+=1))
#{
sample_name1=`echo ${arr[$SGE_TASK_ID - 1]} | awk -F "_" '{print $1}'`
sample_name2=`echo ${arr[$SGE_TASK_ID - 1]} | sed 's/R1/R2/g' | awk -F "_" '{print $1}'`
sample=`echo $sample_name1 | sed 's/.R1.fq.gz//g'`
echo "[mapping running for] $sample_name1"
printf "\n"
echo "bwa mem -t 12 ${reference} ${arr[$SGE_TASK_ID - 1]} ${sample_name2} > $sample.sam" >> "$sample".bam 
read1=$raw_reads/$sample_name1
read2=$raw_reads/$sample_name2
#sample=`echo $sample_name1 | sed 's/fq.gz//g'`
mkdir $sample.logs

bwa mem -R "@RG\tID:$sample\tSM:$sample" -t 12 ${reference} ${read1} ${read2} 2> $sample.logs/bwa.err|\
    /ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/samtools/samtools view -b -F 4 |
    /ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/samtools/samtools sort -m 970m - > "$sample".bam 2> $sample.logs/samsort.err




