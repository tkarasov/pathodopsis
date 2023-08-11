###### THis script is meant to calculate the kinship matrix and variance explained by genotype
cd /ebio/abt6_projects9/pathodopsis_microbiomes/data/
gcta=/ebio/abt6_projects9/pathodopsis_microbiomes/Programs/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
#make a grm
# https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
# poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3
$gcta --bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02_no_mit_no_chlo.vcf.recode_geno0.3  --thread-num 5 --out test --make-grm

# Input the GRM file and output the first 20 eigenvectors for a subset of individuals
$gcta  --grm test --pca 20  --out test 
#--keep test.indi.list 

# Estimate variance explained by genotype

# Take the top 20 eigen vectors and include as covariates in the model for anova in vegan.



# Now build pheno file to then feed to vegan for covariates. 


#  
#  gemma -gk 1 -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
#  -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
#  -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex
#  
