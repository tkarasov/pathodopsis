## Max Planck Institute for Developmental Biology
#### Talia Karasov
This repository contains scripts related to the pathodopsis microbiome study from 2018-2020.


### DEPENDENCIES
#### Read Mapping of Plant and Metagenome
* [python3](https://www.python.org/download/releases/3.0/)
* [bwa (version 0.7.17-r1188)](https://github.com/lh3/bwa)
* [picard (version 2.17.3)](https://broadinstitute.github.io/picard/)
* [samtools (version 1.6-19-g1c03df6(using htslib 1.6-53-ge539b32)](http://www.htslib.org/)
* [centrifuge (version 1.0.4)](https://ccb.jhu.edu/software/centrifuge/manual.shtml)
* [dada2 (version 1.14)](https://benjjneb.github.io/dada2/)


### 16S/ITS Analysis
#### Demultiplexing
* [Bash script to demuliplexing pipeline](https://github.com/tkarasov/pathodopsis/blob/master/code/16S/process_step1_demultiplex_16S_phyllosphere_all.sh)
* [Demultiplexing python script (called by two sequential bash scripts)](https://github.com/tkarasov/pathodopsis/blob/master/code/16S/demultiplex_16S_non_flashed.py)
* [Bash script to clip barcodes](https://github.com/tkarasov/pathodopsis/blob/master/code/16S/clip_barcodes_after_python.sh)

#### Running dada2 to filter data and call ASVs
* [Bash script to run dada2 pipeline](https://github.com/tkarasov/pathodopsis/blob/master/code/16S/process_step2_rundada2_16S.sh)
* [dada2 filter_tlk](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/dada2_filter.R)
* [dada2_inference of ASVs](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/dada2_inference_tlk.R)
* [taxonomy assignment and chimera removal](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/dada2_chimera_taxa.R)
* [generate OTU table](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/convert_dada2_out.R)
* [make filtered OTU table](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/after_dada2_make_otutable_generic.R)
* [phylogeny](https://github.com/tkarasov/pathodopsis/blob/master/code/amplicon_general/after_dada2_make_tree.R)

### Metagenomics

### Climate and Species associations
* [Koppen Geiger climate classification](https://github.com/tkarasov/pathodopsis/blob/master/code/climate/climate_classification.R)
* [General association between whole dataset and climate/species variables](https://github.com/tkarasov/pathodopsis/blob/master/code/modeling_general/permanova_OTU.R)
* [Random forest regression on various climatic variables](https://github.com/tkarasov/pathodopsis/blob/master/code/climate/random_forest.R)

### Genus and ASV-specific associations
* [OTU5 Kriging distribution](https://github.com/tkarasov/pathodopsis/blob/master/code/climate/kriging_OTU5_map.R)
