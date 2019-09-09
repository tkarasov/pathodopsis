#cd ~/links/my_small_projects/MiSeq_083/ampliconsequencing
for i in illumina_M00333_flowcellA_SampleIdH*
do
	echo $i
	cd $i
	plexseq.py -r *R1* --r2 *R2* -i idxfile.txt -o demultiplexed --nostrict > demultiplexed.stdout 2> demultiplexed.stderr
	echo "demultiplexing finished"
	cd demultiplexed
	gzip *.fq
	echo "compression finished"
	for j in *
	do
		foldername=$(basename -s .fq.gz $j)
		mkdir $foldername
		mv $j $foldername
		cd $foldername
		/ebio/abt6_projects9/microbiome_analysis/data/software/seqtk-master/seqtk seq -1 $j | gzip -c > "$foldername"_R1.fq.gz
		/ebio/abt6_projects9/microbiome_analysis/data/software/seqtk-master/seqtk seq -2 $j | gzip -c > "$foldername"_R2.fq.gz
		cd ..
		done
	echo "seqtk finished"
	for l in {01,02,03,04,05,06,07}
	do
		echo $l
		filename=$(echo *$l)
		echo $filename
		cd $filename
		/ebio/abt6_projects9/microbiome_analysis/data/software/bwa-0.7.15/bwa mem -t 4 /ebio/abt6_projects7/small_projects/esymeonidi/MiSeq_010_crispr/bwainddexes/gICS1/gICS1 *R1.fq.gz *R2.fq.gz 2> bwa.stderr > "$filename".sam
		echo "bwa finish"
		/ebio/abt6/esymeonidi/Software/bin/samtools view -bS "$filename".sam > "$filename".bam
		echo "view finish"
		rm "$filename".sam
		echo "sam remove"
		/ebio/abt6/esymeonidi/Software/bin/samtools sort -@ 4 "$filename".bam > S"$filename".bam 2> samsort.stderr
		rm "$filename".bam
		echo "sort finished"
		/ebio/abt6/esymeonidi/Software/bin/samtools index S"$filename".bam 2> samindex.stderr
		/ebio/abt6/esymeonidi/Software/bin/samtools flagstat S"$filename".bam > S"$filename".stat 2> samflags.stderr
		/ebio/abt6_projects9/microbiome_analysis/data/software/freebayes_1.1/bin/freebayes -f /ebio/abt6_projects7/small_projects/esymeonidi/MiSeq_010_crispr/bwainddexes/gICS1.fasta S"$filename".bam > S"$filename".vcf
		echo "freebayes finish"
		cd ..
		done
	for k in {08,09,10}
	do
		filename=$(echo *$k)
		echo $filename
		cd $filename
		/ebio/abt6_projects9/microbiome_analysis/data/software/bwa-0.7.15/bwa mem -t 4 /ebio/abt6_projects7/small_projects/esymeonidi/MiSeq_010_crispr/bwainddexes/gFLC/gFLC *R1.fq.gz *R2.fq.gz 2> bwa.stderr > "$filename".sam
		echo "bwa finished"
		/ebio/abt6/esymeonidi/Software/bin/samtools view -bS "$filename".sam > "$filename".bam
		echo "view finished"
		rm "$filename".sam
		echo "sam remove"
		/ebio/abt6/esymeonidi/Software/bin/samtools sort -@ 4 "$filename".bam > S"$filename".bam 2> samsort.stderr
		rm "$filename".bam
		echo "sort finished"
		 /ebio/abt6/esymeonidi/Software/bin/samtools index S"$filename".bam 2> samindex.stderr
		 /ebio/abt6/esymeonidi/Software/bin/samtools flagstat S"$filename".bam > S"$filename".stat 2> samflags.stderr
		/ebio/abt6_projects9/microbiome_analysis/data/software/freebayes_1.1/bin/freebayes -f /ebio/abt6_projects7/small_projects/esymeonidi/MiSeq_010_crispr/bwainddexes/gFLC.fasta S"$filename".bam > S"$filename".vcf
		echo "freebayes finished"
		cd ..
		done
cd ../..
done 
