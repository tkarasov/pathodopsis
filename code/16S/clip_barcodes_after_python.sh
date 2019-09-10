#pull the bad part from every read and rename read
output_direc=$1

usearch=/ebio/abt6_projects9/metagenomic_controlled/Programs/metagenomics_pipeline_software/bin/usearch11.0.667_i86linux32
for var1 in 1 2 3 4 5 6;
do
    for var2 in 1 2 3 4 5 6;
    do
        my_files=`ls | grep 341F${var1}_806F${var2}_R1`
        for file in $my_files;
        do
            forward_strip=$(($var1+18))
            reverse_strip=$(($var2+21))
            echo $file
            reverse=`echo $file | sed -e 's/R1/R2/g'`
            $usearch -fastx_truncate $file -stripleft $forward_strip -fastqout strip_$file
            $usearch -fastx_truncate $reverse -stripleft $reverse_strip -fastqout strip_$reverse
        done
    done
done


#$usearch -fastx_truncate "$merged"_temp.fq -stripleft 20 -stripright 26  -fastqout "$merged"_temp.fq.stripped

#once these are done rename all of the files
info=/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/16S/16S_12_2018/16S_plate_locations.csv

mkdir $output_direc/demult_python
#cd $output_direc/demult_python
#read through every line of info file, and get renaming
while IFS= read -r line;
do
    read -r f1 f2 f3 f4 <<<"$line"
    platepos=$f3
    barcodeF=`echo $f4 | cut -f1 -d '.'`
    barcodeR=`echo $f4 | cut -f2 -d '.'`
    demult_file_R1=`ls | grep strip | grep "R1_$platepos.fastq" | grep "$barcodeF" | grep "$barcodeR"`
    demult_file_R2=`ls | grep strip | grep "R2_$platepos.fastq" | grep "$barcodeF" | grep "$barcodeR"`
    echo $demult_file_R2
    cp $demult_file_R1 $output_direc/demult_python/${f1}_16S_R1.fastq
    cp $demult_file_R2  $output_direc/demult_python/${f1}_16S_R2.fastq
done<$info

