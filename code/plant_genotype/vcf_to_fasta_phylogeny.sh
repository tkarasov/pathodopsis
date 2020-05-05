# The goal of this script is to take the vcf files and get fasta sequences from them. It's laborious and weirdly hard but the key is using a mix of bcf and vcftools
#AT4G14400 is ACD6
# position of most significant SNPs are 8295146, then also 8295218 and 8295231

############################################################
# Sorting references
############################################################
#rename TAIR10 chromosomes 
cat TAIR10_chr_all.fas | sed 's/>/>chr/g' > TAIR10_chr_all_recode.fas

# run snpEff in the snpEff.jar directory
java -jar /ebio/abt6_projects8/Pseudomonas_mapping/Programs/snpEff_latest_core/snpEff/snpEff.jar testAthalianaTair10 $output_dir/poolsGVCF.filtered_snps_final.PASS.bi.vcf > $output_dir/poolsGVCF_SNPeff_annotation.vcf

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

intersectBed -a poolsGVCF.filtered_snps_final.PASS.bi.vcf -b TAIR10_gander_20kb.bed > poolGVCF_gander.vcf

echo -e "chr4\t8294164\t8299195" >acd6.bed
#echo -e "chr4\t1\t8315146" >acd6.bed

# Take only the region around acd6
vcftools --vcf poolsGVCF.filtered_snps_final.PASS.bi.vcf \
	--chr chr4 --from-bp 8294164 --to-bp 8299195 --recode \
       	--out acd6

############################################################
# Perform fst analysis
############################################################
vcftools --vcf poolsGVCF.filtered_snps_final.PASS.bi.vcf \
	--weir-fst-pop cluster1.pop --weir-fst-pop cluster3.pop --out fst_1_3

vcftools --vcf poolGVCF_gander.vcf \
	--weir-fst-pop cluster1.pop --weir-fst-pop cluster3.pop --out poolGVCF_gander_fst_1_3

############################################################
# Convert acd6 to tped and tfam
############################################################

vcftools --plink-tped --out acd6_patho.vcf --vcf acd6.recode.vcf


vcftools --gzvcf acd6.recode.vcf.gz --freq --out allel_freq
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
# Now it's time to align the sequences with the Todesco sequences
############################################################
cat acd6_all_samples.fasta acd6_todesco.fasta > all_acd6_mine_todesco.fasta
clustalw2 all_acd6_mine_todesco.fasta > all_acd6_mine_todesco_aligned.fasta
raxmlHPC-PTHREADS -T 20 -m PROTGAMMALGF -n raxml_all_acd6_mine_todesco_aligned -s all_acd6_mine_todesco_aligned.fasta


