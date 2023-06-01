# the original name of this script is "run_fragment_placement.sh" which resides in qiime folder. This is how I overlapped Ben Brachi reads of v2/v3 with my own. 
#script to run sepp with Brachi data
# https://library.qiime2.org/plugins/q2-fragment-insertion/16/

cd /Users/talia/Dropbox/pathodopsis/brachi_analysis/qiime

#
# wget -O "sepp-refs-gg-13-8.qza"   "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"

# R write refernece file to fasta
#  writeXStringSet(x =list_of_DNAStrings, filepath = "/Users/talia/Dropbox/pathodopsis/brachi_analysis/qiime/patho_16S_repseq.fasta", 
#    format = "fasta")
# 


conda activate q2-2021.8

#make qza for ben input fasta and my top ASVs
qiime tools import --input-path /Users/talia/Dropbox/pathodopsis/brachi_analysis/brachi_data_ind_download/16S_repseq.fasta --output-path brachi_16S_repseq.qza --type 'FeatureData[Sequence]'
qiime tools import --input-path /Users/talia/Dropbox/pathodopsis/brachi_analysis/qiime/patho_16S_repseq.fasta --output-path patho_16S_repseq.qza --type 'FeatureData[Sequence]'


#now run fragment-insertion sepp on brachi samples
qiime fragment-insertion sepp \
  --i-representative-sequences brachi_16S_repseq.qza \
  --i-reference-database sepp-refs-gg-13-8.qza \
  --o-tree brachi-insertion-tree.qza \
  --o-placements brachi-insertion-placements.qza


#now run fragment-insertion sepp on patho samples
 qiime fragment-insertion sepp \
  --i-representative-sequences patho_16S_repseq.qza \
  --i-reference-database sepp-refs-gg-13-8.qza \
  --o-tree patho-insertion-tree.qza \
  --o-placements patho-insertion-placements.qza 

# now filter (why?)
qiime fragment-insertion filter-features \
  --i-table table.qza \
  --i-tree insertion-tree.qza \
  --o-filtered-table filtered_table.qza \
  --o-removed-table removed_table.qza



  #or sortmerna
#  https://github.com/biocore/sortmerna(base) tkarasov@taco:/ebio/abt6_projects9/pathodopsis_microbiomes/brachi_analysis/qiime$ 

