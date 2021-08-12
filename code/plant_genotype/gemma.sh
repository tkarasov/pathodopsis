###### THis script is meant to calculate the 





# remove mitochondra. This works
vcftools --gzvcf /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.gz \
	 --not-chr mit --not-chr chlo --out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf --recode 

#generate plink king-table
# plink2 --vcf   poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --recode \
#  --make-king square --pheno phenotype3.txt

######THIS IS THE ONE THAT WORKS 07/19/2021

# Generate bedfiles for lmekin

#This one filters at least 70% of genotyping rate per SNPS (~3000 SNPs)
plink --vcf   poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --geno 0.3 \
--out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3 --recode --make-bed #pve = 0.544807 in gemma


#This one filters at least 80% genotyping rate per SNP in gemma (~1000 SNPs)
plink --vcf   poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --geno 0.2 \
--out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.2 --recode --make-bed #pve = 1966 x 10^-6 in gemma


# Generate GRM from plink with a genotyping rate of 70%
plink2 --vcf   poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --recode \
 --make-rel square --pheno phenotype5.txt --geno 0.3

# run gemma (with genotyping rate of 70% get a h2 of 39%) 
gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  plink2.rel -o plink_rel_output_var 



gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  plink2.king -o plink_king_output_var 

gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k \
 /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.cXX.txt \
 -o plink_king_output_var 


######Work on 08/05/2021 with Gautam's new matrix
# run gemma (with genotyping rate of 70% get a h2 of 39%) 
poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.reheader.allChrs.oEAs_bgl.vcf.gz
gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k \
 /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.reheader.allChrs.oEAs_bgl.vcf.gz.cXX.txt \
-o plink_rel_output_var_impute





# Generate GRM from plink while filtering on a certain genotyping rate
plink --vcf   poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --geno 0.3 --mind 0.7 \
--out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3_mind0.7 --recode --make-bed





gemma -gk 1 -maf 0.01 -bfile /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK \
 	-a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
 	-p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex 






# make necessary files for gemma
# plink --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --make-bed --out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
#  --pheno phenotype2.txt

#plink --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --recode --pheno phenotype3.txt
plink --vcf /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.gz \
 --recode --pheno phenotype3.txt

plink --file plink --make-bed

#this one is running. Started at 10:55AM
# gemma -gk -maf 0.01 -bfile plink
# vcftools  --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --recode --out plink
# gemma -gk -maf 0.1 -bfile /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt
 # maek kinship matrix in gemma. This hangs forever
 
 #cannot find bim/bed/fam for this
 gemma -gk -maf 0.01 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
 	-a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
 	-p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex 
# Getting ready to run Gemma


 #cannot find bim/bed/fam for this
 gemma -gk 1 -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf \
 -o poo
 

 #cannot find bim/bed/fam for this
 gemma -gk 1 -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK \
 -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex
 

gemma -vc 2 -n 2 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  plink2.king -o plink_king_output_var 

# https://www.cog-genomics.org/plink/1.9/distance
 
gemma -vc 2 -n 2 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  plink2.rel -o plink_rel_output_var 





# Gautam's stwuff
g_file=/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/
gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  $g_file/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.rel -o rel_id_output_var 



gemma -vc 2 -n 1 -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -k  $g_file/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.TK.rel -o rel_id_output_var 













# #################################
# # Reorder according to genotype file
# #################################
# pref = "/ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf"
# #system("zcat /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf | grep ^#CHROM > genot.txt")
# #plink --file /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf --recode bimbam --out mybimbam
# 

# # remove mitochondra
# vcftools --gzvcf /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.gz --not-chr mit --not-chr chlo --out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf --recode 
# 
# # make necessary files for gemma
# plink --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --make-bed --out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
# --pheno phenotype2.txt
# 
# # vcftools --vcf plink --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --make-bed --out plink --vcf poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode.vcf --make-bed --out poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf --plink
# # gemma -gk -maf 0.1 -bfile /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt
#  gemma -gk -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex
# # Getting ready to run Gemma
#  
#  gemma -gk 1 -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
#  -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
#  -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex
#  
