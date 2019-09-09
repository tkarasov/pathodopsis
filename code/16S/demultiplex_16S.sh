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
#
#  Run job from current working directory
#$ -cwd
#
#  Send email when the job begins, ends, aborts, or is suspended
##$ -m beas

merged=$merged


usearch=/ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/bin/usearch11.0.667_i86linux32

cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/demultiplexed_reads

#341 f1=19, f2=20, f3=21, f4=22, f5=23, f6=24
#806 f1=22, f2=23, f3=24, f4=25, f5=26, f6=27

#Plate1: Ts341F_f2 + 806R_f5
#341_f2 with #806_f5
egrep "^TGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAGTC$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 20 -stripright 26  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F2.806F5.fastq
rm "$merged"_temp*

#Plate2: Ts341F_f3 + 806R_f4
#341_f3 with #806_f4
egrep "^CTGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAGT$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 21 -stripright 25  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F3.806F4.fastq
rm "$merged"_temp*

#Plate3: Ts341F_f4 + 806R_f3
#341_f4 with #806_f3
egrep "^ACTGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAG$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 22 -stripright 24  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F4.806F3.fastq
rm "$merged"_temp*

#Plate4: Ts341F_f5 + 806R_f2
#341_f5 with #806_f2
egrep "^GACTGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTA$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 23 -stripright 23  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F5.806F2.fastq
rm "$merged"_temp*

#Plate5: Ts341F_f6 + 806R_f1
#341_f6 with #806_f1
egrep "^TGACTGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGT$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 24 -stripright 22  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F6.806F1.fastq
rm "$merged"_temp*

#Plate7: Ts341F_f1 + 806R_f5
#341_f1 with #806_f6*******f5
egrep "^GACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAGTC$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 19 -stripright 26  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F1.806F5.fastq
rm "$merged"_temp*

#Plate8: Ts341F_f1 + 806R_f6
#341_f1 with #806_f6
egrep "^GACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAGTCA$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 19 -stripright 27  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F1.806F6.fastq
rm "$merged"_temp*

#8single:Ts341F_f2 + 806R_f4
#341_f2 with #806_f4
egrep "^TGACCTACGGGAGGCAGCAG|ATTAGA[ATCG]ACCC[ATCG][ATCG]GTAGTCCGTAGT$" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/flashed_reads/$mergedSSS -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 20 -stripright 25  -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_341F2.806F4.fastq
rm "$merged"_temp*
