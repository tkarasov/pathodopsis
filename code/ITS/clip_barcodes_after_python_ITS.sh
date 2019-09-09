#!/bin/sh

#pull the bad part from every read and rename read
usearch=/ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/bin/usearch11.0.667_i86linux32
for var1 in 1 2 3 4 5 6;
do
    for var2 in 1 2 3 4 5 6;
    do
        my_files=`ls | grep AgITSF${var1}_AgITSR${var2}_R1 | grep -v strip`
        for file in $my_files;
        do
            forward_strip=$(($var1+23))
            reverse_strip=$(($var2+21))
            echo $file
            reverse=`echo $file | sed -e 's/_R1_/_R2_/g'`
            $usearch -fastx_truncate $file -stripleft $forward_strip -fastqout strip_$file
            $usearch -fastx_truncate $reverse -stripleft $reverse_strip -fastqout strip_$reverse
        done
    done
done

$usearch -fastx_truncate "$merged"_temp.fq -stripleft 20 -stripright 26  -fastqout "$merged"_temp.fq.stripped

#once these are done rename all of the files
info=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/ITS_plate_locations.txt

cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python
#read through every line of info file, and get renaming
while IFS=  read -r line;
do
    read -r f1 f2 f3 f4 <<<"$line"
    platepos=$f3
    barcodes=$f4
    demult_file_R1=`ls | grep strip | grep "R1_$platepos.fastq" | grep $barcodes`
    demult_file_R2=`ls | grep strip | grep "R2_$platepos.fastq" | grep $barcodes`
    echo $demult_file_R2
    echo "That was supposed to be the demult"
    echo $platepos, $barcodes
    cp $demult_file_R1  /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/${f1}_ITS_R1.fastq
    cp $demult_file_R2  /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/ITS_1_2019/demult_python/${f1}_ITS_R2.fastq
done<$info

