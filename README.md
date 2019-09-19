# pathodopsis
Metagenomic and 16S analysis of Eurasian collections in 2018

## 16S and ITS analysis
16S and ITS were analysed using DADA2
A modified version of the pipeline and scripts detailed here: https://github.com/LangilleLab/microbiome_helper/wiki/DADA2-16S-Chemerin-Tutorial were used for the analysis.

### Amplicon analysis with DADA2

We created a conda environment titled pathodopsis

```console
conda activate pathodopsis
```
And we demultiplexed our multi-multiplexed samples first using the in-house sample-sheet demultiplexing (for the 96-well demultiplexing) and next using barcodes.

The demultiplexing and pooling of samples across runs is done in the pipeline detailed here:
```
/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/process_step1_demultiplex_16S.sh
```
And the dada2 pipeline is here:
```
/ebio/abt6_projects9/pathodopsis_microbiomes/pathodopsis_git/code/16S/process_step2_rundada2_16S.sh
```
But Ill describe step-by-step below. The first issue is that our samples were demultiplexed into the 96-well format, but we need to demultiplex the nested barcodes within this format. This is accomplished by the following script.
```console
demultiplex_16S_non_flashed.sh <input_raw_read_direc> <output_direc>
```
Now we have a nicely demultiplexed directory and we want to remove the barcodes from the fastqs.
```console
clip_barcodes_after_python.sh <output_direc>
```
Depending on whether the sample was sequenced across multiple lanes or not, we will either pool the clipped reads or proceed directly to dada2. Once the samples are ready in fastq form (pooled or not) we will transition to the process_step2 script in which we run dada2.

At this point we have demultiplexed samples with barcodes removed ready for filtering and dada2 sequence analysis. We will use primarily scripts (and slightly modified versions) from the Langille tutorial. First filter reads at a certain truncation length

```console
dada2_filter_tlk.R \
        -f <input folder> --truncLen <trunc_from, trunc_to> --maxN <# N's allowed which should be 0> \
        --maxEE 3,7 --truncQ 2 --threads 9 --f_match _R1_.*fastq.gz \
        --r_match _R2_.*fastq.gz
```

Now we are ready to infer the error models
```console
 dada2_inference_tlk.R \
        -f filtered_fastqs/ --seed 1659 -t 9 --verbose --plot_errors
 ```
 
 Lastly we want to remove chimeras and assign taxonomy
 ```console
 dada2_chimera_taxa.R -i seqtab.rds \
  -r <path to reference sequences for taxonomy assignment> \
  -t 9 --skip_species
  ```
  
  Then to reduce clutter we merge logfiles
  ```console
  merge_logfiles.R \
        -i cutadapt_log.txt,dada2_filter_read_counts.txt,dada2_inferred_read_counts.txt,dada2_nonchimera_counts.txt \
        -n cutadapt,filter,infer,chimera -o combined_log.txt
  ```
  
  Lastly convert ot biom format
  ```console
  convert_dada2_out.R -i \
         seqtab_final.rds -b seqtab.biom.tsv -f seqtab.fasta --taxa_in tax_final.rds
```
### Analysis of OTU table
Once we have the OTU table generated and filtered using dada2, we want to analyse our microbiomes. I've written a series of Rscripts that will analyse the data. Again these scripts can be submitted to the cluster with a shell script

The following script takes the dada2 output and will filter the table to desired depths and normalize for relative abundance across the table. This script also makes comparisons between the metagenome data and the 16S data

```console
after_dada2_make_otu_table.Rmd
```
Now we want to start look at differences in the microbiomes

```console
after_dada2_do_pca.Rmd
```
### Checking for oddities of metagenomic table
In general I have found more Burkholderia in the metagenomes than in the 16S data. I would like to be able to pull the reads assigned to a given category to assess who they are (and whether the classification makes sense)

```console
./code/metagenome/pull_taxonomic_reads.py <input_R1.fastq> <taxonomic_family>
```
This output a file (This output a file (/quality_control_metagenome/centrifuge_inspect/single_family_seq.fa) with the reads 
