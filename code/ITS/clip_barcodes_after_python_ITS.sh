#!/bin/sh

#pull the bad part from every read and rename read
usearch=/ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/bin/usearch11.0.667_i86linux32
output_direc=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/all_runs 

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
info=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/plate_location_map_16S_ITS1.txt

cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/all_runs
mkdir /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/ITS/all_runs/demult_python

#read through every line of info file, and get renaming
while IFS= read -r line;
do
    read -r f1 f2 f3 f4 f5 f6 f7 <<<"$line"
    wrong_platepos=$f4
    correct_platepos=$f3
    barcodeF=`echo $f7 | cut -f1 -d '_'`
    barcodeR=`echo $f7 | cut -f2 -d '_'`
    demult_file_R1=`ls | grep strip | grep "R1_$wrong_platepos.fastq" | grep "$barcodeF" | grep "$barcodeR"`
    demult_file_R2=`ls | grep strip | grep "R2_$wrong_platepos.fastq" | grep "$barcodeF" | grep "$barcodeR"`
    echo "output is:"$output_direc/demult_python/${f1}_${f2}_${f3}_${f5}_ITS_R1.fastq
    echo "demult_file_R1:"$demult_file_R1
    #if [ -z "$var" ]
    #then
    #    break
    #fi
    cp $demult_file_R1 $output_direc/demult_python/${f1}_${f2}_${f3}_${f5}_16S_R1.fastq
    cp $demult_file_R2  $output_direc/demult_python/${f1}_${f2}_${f3}_${f5}_16S_R2.fastq
done< $info