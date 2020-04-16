#!/bin/sh
export diamond=/ebio/abt6_projects9/pathodopsis_microbiomes/data/metagenome_modelling/ref_database/diamond
SAMPLE=/ebio/abt6_projects9/pathodopsis_microbiomes/data/metagenome_modelling/here.fq
#/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run104_2018_7_metagenome_reads/PA0882_RunId0104_LaneId5MetagenomicR1R2.fq
OUTPUT_DIR=/ebio/abt6_projects9/pathodopsis_microbiomes/data/metagenome_modelling
humann2 --input $SAMPLE --output $OUTPUT_DIR
