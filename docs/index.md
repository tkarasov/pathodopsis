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
* [dada2 filter]
* [dada2_inference of ASVs]
* [taxonomy assignment]
* [phylogeny]

### Metagenomics

### Spatial Analysis
