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
$gcta --reml --grm test --pheno phenotype3.txt --out test
# Summary result of REML analysis (phenotype is loading on PC1):
# Source  Variance        SE
# V(G)    0.010033        0.001269
# V(e)    0.009107        0.001430
# Vp      0.019140        0.001290
# V(G)/Vp 0.524175        0.062331


# Take the top 20 eigen vectors and include as covariates in the model for anova in vegan.
# eigenvectors are in test.eigenvec. First two columns are individual and family ID. These can be fed into R with the other microbe data
# now run vegan
run_vegan_reml.R


#  
#  gemma -gk 1 -maf 0.1 -bfile poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf \
#  -a /ebio/abt6_projects9/pathodopsis_HpA/host/pathodopsis_Microbiome/postVcf/poolsGVCF.filtered_snps_final.PASS.bi.cleaned.maf0_02.vcf.map \
#  -p /ebio/abt6_projects9/pathodopsis_microbiomes/data/phenotype2.txt -o ex
#  
