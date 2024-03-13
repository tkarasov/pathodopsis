# The goal of this script is to take the vcf files and get fasta sequences from them. It's laborious and weirdly hard but the key is using a mix of bcf and vcftools
#AT4G14400 is ACD6
# position of most significant SNPs are 8295146, then also 8295218 and 8295231 and 8295845
cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/
############################################################
# Sorting references
############################################################
#rename TAIR10 chromosomes 
cat TAIR10_chr_all.fas | sed 's/>/>chr/g' > TAIR10_chr_all_recode.fas

# run snpEff in the snpEff.jar directory
java -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar testAthalianaTair10 $output_dir/poolsGVCF_lowCov.filtered_snps_final.PASS.bi.vcf.gz > $output_dir/poolsGVCF_SNPeff_annotation.vcf

#convert gff to bed
sortBed -i TAIR10_GFF3_genes.gff | gff2bed > TAIR10_GFF3_genes.bed



############################################################
# Start subsetting the files to the regions of interest
############################################################
# Keep only gander genes (immunity genes)
grep -f defense_genes_gander.txt TAIR10_GFF3_genes.bed > TAIR10_gander.bed

# replace chromosome name
sed -i 's/Chr/chr/g' TAIR10_gander.bed
awk -v OFS='\t' -v s=20000 'BEGIN {FS="\t"}{print $1, $2-s, $3+s, $4, $5, $6, $7, $8, $9, $10}' TAIR10_gander.bed \
| awk -v OFS='\t' 'BEGIN {FS="\t"}{$2=($2<0)?1:2}1' > TAIR10_gander_20kb.bed

# Select vcf file for only those genes in gander and 20kb around those genes
intersectBed -a TAIR10_gander_20kb.bed -b poolsGVCF_SNPeff_annotation.vcf > poolGVCF_gander.bed

intersectBed -a poolsGVCF_lowCov.filtered_snps_final.PASS.bi.vcf.gz -b TAIR10_gander.bed -header \
> poolGVCF_gander2.vcf

echo -e "chr4\t8294164\t8299195" >acd6.bed
#echo -e "chr4\t1\t8315146" >acd6.bed

# Take only the region around acd6
bcftools index poolsGVCF_lowCov.filtered_snps_final.PASS.bi.vcf.gz
bcftools view \
	--regions chr4:8294164-8299195  \
       	--output-file acd6.vcf poolsGVCF_lowCov.filtered_snps_final.PASS.bi.vcf.gz

############################################################
# Perform fst analysis
############################################################
vcftools --gzvcf poolsGVCF_lowCov.filtered_snps_final.PASS.bi.vcf.gz \
	--weir-fst-pop cluster1x.pop --weir-fst-pop cluster2x.pop --out fst_1_2

vcftools --vcf poolGVCF_gander2.vcf \
	--weir-fst-pop cluster1x.pop --weir-fst-pop cluster2x.pop --out poolGVCF_gander_fst_1_2

############################################################
# Convert acd6 to tped and tfam
############################################################

vcftools --plink-tped --out acd6_patho.vcf --vcf acd6.vcf


vcftools --vcf acd6.vcf --freq --out allel_freq
#output tsv (will be used for list of samples below)
vk vcf2tsv long acd6.recode.vcf.gz > acd6.tsv

#convert into fasta with N's       	
vk phylo fasta acd6.recode.vcf.gz > acd6.samples.fasta

# get this same region from the TAIR10 geneom
#bedtools getfasta -fi TAIR10_chr_all_recode.fas -bed acd6.bed  | sed "s/chr4/chr4\t/g"> acd6_TAIR10.fasta

# Now prepare sample for bedtools
bgzip -c acd6.recode.vcf > acd6.recode.vcf.gz
bcftools index acd6.recode.vcf.gz 

#run a for loop outputting the desired fasta
for samp in `cat  acd6.tsv | cut -f 28 | sort | uniq`;
do
	cat TAIR10_chr_all_recode.fas | sed 's/>chr4.*/>chr4/g' | vcf-consensus -s $samp acd6.recode.vcf.gz  > temp.$samp
	bedtools getfasta -fi temp.$samp -bed acd6.bed  | sed "s/chr4/$samp/g"> $samp.acd6_copy.fasta 
	echo $samp 
	rm temp.$samp &
done

cat *.acd6_copy.fasta  | sed 's/:/\t/g'> acd6_all_samples.fasta



############################################################
# Look at 1001. The resulting vcf works to be read in by IBS.R
############################################################
#I downloaded the acd6 locus from the lab site
#I couldn't get gatk to work for this. Was outputting NA in files
# gatk VariantsToTable \
# 	-V acd6_74936542461bf48799441b87dd6ebece_fullgenome.vcf \
# 	-O output.table -disable-sequence-dictionary-validation \
# 	-GF GT -F POS=8295146

vcftools --vcf acd6_74936542461bf48799441b87dd6ebece_fullgenome.vcf --max-alleles 2 --remove-indels --recode \
 --out biallelic_acd6_74936542461bf48799441b87dd6ebece_fullgenome.vcf 

vcftools --plink-tped --out acd6_1001 --vcf biallelic_acd6_74936542461bf48799441b87dd6ebece_fullgenome.vcf.recode.vcf

############################################################
# Now let's do a PCA of the acd6 alleles with Todesco's
############################################################
#let's get a list of todesco alleles
cat acd6_todesco.fasta | grep "^>" | \
awk -F "ecotype" '{print $2}' | awk -F "At" '{print $1}' \
| awk -F "accelerated" '{print $1}' | awk -F "ankyrin" '{print $1}' > todesco_list.txt

#Now get the id numbers for these strains
cat todesco_list_id.txt | grep -v NA | grep -v "# "| cut -f 2 > samp_list.txt

#Now reduce the 1001 vcf
bcftools view -s 6897,6958,8337,6932,7169,7125,6961,7328,8264,6188,6911,6909,7127,9941,6940,6043,6929,6919 \
acd6_74936542461bf48799441b87dd6ebece_fullgenome.vcf > acd6_1001_subset.vcf

#prep the indeces and such for the merged vcf
sed -i 's/chr4/4/g' acd6.vcf
bgzip -c acd6.vcf > acd6.vcf.gz
bgzip -c acd6_1001_subset.vcf > acd6_1001_subset.vcf.gz
bedtools intersect -a acd6.vcf.gz -b acd6_1001_subset.vcf | cut -f2 > keep_sites.txt
bcftools index acd6.vcf.gz
gatk IndexFeatureFile -F acd6_1001_subset.vcf.gz

bcftools merge --merge all acd6.vcf.gz \
acd6_1001_subset.vcf --output merged_acd6.vcf

gatk VariantsToTable \
	-V acd6.vcf -GF GT --disable-sequence-dictionary-validation \
	--output my_acd6_variants_table.txt -F POS

gatk VariantsToTable \
	-V acd6_1001_subset.vcf.gz -GF GT --disable-sequence-dictionary-validation \
	--output thousand_acd6_variants_table.txt -F POS

head -1 my_acd6_variants_table.txt > my_keep_table.txt
grep -F -f keep_sites.txt my_acd6_variants_table.txt | sed 's*\./.*NA*g' >> my_keep_table.txt
head -1 thousand_acd6_variants_table.txt > thous_keep_table.txt
grep -F -f keep_sites.txt thousand_acd6_variants_table.txt| sed 's*\./.*NA*g' >> thous_keep_table.txt

cut -f2- my_keep_table.txt | sed "s#/#|#g" > temp
paste thous_keep_table.txt temp > merged_tables.txt
rm temp

############################################################
# Now it's time to subset the acd6 whole genome mappings just to the acd6 region
############################################################
#https://wurmlab.github.io/genomicscourse/2016-SIB/practicals/population_genetics/map_call
#samtools merge merged.bam -R chr4:8293409-8299949 *.bam
# Step 1: samtools mpileup
## Create index of the reference (different from that used by bowtie2)
reference="/ebio/abt6_projects9/pathodopsis_microbiomes/data/plant_genotype/TAIR10_chr_all_recode.fas"
ln -rs $reference .
samtools faidx TAIR10_chr_all_recode.fas

# Run samtools mpileup
bcftools mpileup -r chr4:8293409-8299949 -Q20 --fasta-ref TAIR10_chr_all_recode.fas   *.bam > raw_calls_acd6.bcf

# Run bcftools call
bcftools call  -v -m raw_calls_acd6.bcf > my_acd6_calls.vcf

sed -i 's/.bam//g' my_acd6_calls.vcf

#now calculate fst
vcftools --vcf my_acd6_calls.vcf \
	--weir-fst-pop ../cluster1x.pop --weir-fst-pop ../cluster3x.pop --out fst_my_acd6

#vcf2bed < my_acd6_calls.vcf > out.bed

plink --vcf my_acd6_calls.vcf --out out
plink --bfile out --r2 --ld-window-kb 100 -r2 bin > my_acd6_calls.ld
plink --bfile out --blocks 'no-pheno-req'
plink --bfile out --hap plink.blocks --hap-freq


     CHROM     POS WEIR_AND_COCKERHAM_FST
270   chr4 8294681              0.0859257
345   chr4 8295157              0.1041020
352   chr4 8295218              0.0823132
356   chr4 8295231              0.1567320
370   chr4 8295289              0.1286460
373   chr4 8295305              0.1208780
380   chr4 8295338              0.1619880
383   chr4 8295350              0.1117610
386   chr4 8295362              0.1613770
403   chr4 8295477              0.0834420
462   chr4 8295844              0.1697820***high frequency. Fst=0.544 globally. THIS IS THE BETTER ONE!
463   chr4 8295845              0.0961080
772   chr4 8297296              0.0847506
953   chr4 8297812              0.0825090
954   chr4 8297824              0.1077450***missense, 336 individuals
955   chr4 8297825              0.1335070***missense, 326 individuals, Fst=0.25 globally
1368  chr4 8299576              0.0836269
1390  chr4 8299690              0.0897623



############################################################
# Now it's time to align the sequences with the Todesco sequences
############################################################
cat acd6_all_samples.fasta acd6_todesco.fasta > all_acd6_mine_todesco.fasta
clustalw2 all_acd6_mine_todesco.fasta > all_acd6_mine_todesco_aligned.fasta
raxmlHPC-PTHREADS -T 20 -m PROTGAMMALGF -n raxml_all_acd6_mine_todesco_aligned -s all_acd6_mine_todesco_aligned.fasta


