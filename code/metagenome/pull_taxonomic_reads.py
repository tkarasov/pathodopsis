#!/usr/bin/env python

#this script takes the output from centrifuge and pull the reads belonging to a taxonomic category, then blasts them against nr

from ete3 import NCBITaxa
import os
ncbi = NCBITaxa()

def get_family(taxid, desired_ranks):
    try:
        lineage = ncbi.get_lineage(int(taxid))
        names = ncbi.get_taxid_translator(lineage)
        lineage2ranks = ncbi.get_rank(names)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        family = names[ranks2lineage["family"]]
        return(family)

    except ValueError:
        return("none")
    except KeyError:
        return("none")

def get_readnames(classification, cent_out):
    class_dict = {}
    for line in cent_out[1:]:
        if line[2]=='0':
            pass
        else:
            tax = get_family(line[2], "family")
            print(tax)
            if tax == classification:
                try:
                    class_dict[line[2]].append(line[0])
                except KeyError:
                    class_dict[line[2]] = [line[0]]
            else:
                pass
    return class_dict


def pull_reads(read1_file):
    #use seqtk to pull the desired reads
    bashCommand = "seqtk subseq "+read1_file+" /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/name.1st > /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/single_family_seq.fq"
    os.system(bashCommand)

    os.system("seqtk seq -a /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/single_family_seq.fq > /ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/single_family_seq.fa")


def blast_nr(input_fasta, output_name):
    bashCommand = 'blastn -outfmt "6 qseqid sseqid bitscore qstart sstart" -db /ebio/abt6_projects9/abt6_databases/db/blast/prebuilt/nt -query '+input_fasta+' -out '+output_name
    os.system(bashCommand)


#read in the centrifuge output file
cent_out=[line.strip().split() for line in open("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/run116_2018_9_metagenome_reads/centrifuge_output/PA0913.R1.fq.out").readlines()]

#the first column is the read name and the second column is the taxonomic classificiation
classification="Burkholderiaceae"
class_dict = get_readnames(classification, cent_out)

#write file with names of reads that I want to pull
all_reads = []

for clas in class_dict.values():
    for read in clas:
        all_reads.append(read)

all1=list(set(all_reads))

with open("/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/name.1st", "w") as f:
    for rec in all_reads:
        f.write("%s\n" %rec)

read1_file="/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/metagenome/2018_9_metagenome_reads/illumina_ST-J00101_flowcellA_SampleIdPA0913_RunId0116_LaneId3/Pathodopsis_Metagenome_S185_L003_R1_001.fastq.gz"

pull_reads(read1_file)

output_name="/ebio/abt6_projects9/pathodopsis_microbiomes/data/processed_reads/metagenome/quality_control_metagenome/centrifuge_inspect/single_family_seq_blast_out.txt"
blast_nr(input_fasta, output_name)
#Now we have a file called single_family_seq.fa which we would like to blast against nr for classification

