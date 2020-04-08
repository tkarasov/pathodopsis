# The goal of this script is to pull sequences of ITS1 of Oomycetes to build a database of Oomycete sequences for dada2
# /ebio/abt6_projects9/metagenomic_controlled/Programs/anaconda3/envs/pathodopsis/bin/ipython
# http://biopython.org/DIST/docs/tutorial/Tutorial.pdf
import Bio
from Bio import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from Bio import Entrez
import os
import csv
from ete3 import NCBITaxa
import pytaxize

# taxid for Oomycete
ncbi = NCBITaxa()
TAXLIST = ("4762")


def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    # lineage2ranks = dict((taxid, rank) for taxid, rank in lineage2ranks.items() if rank != 'no rank')
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items() if rank != "no rank")
    ranks2lineage2 = dict((rank, ncbi.get_taxid_translator([str(taxid)]).values()[0]) for (rank, taxid) in ranks2lineage.items())
    # ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage2.get(rank, 'NA') for rank in desired_ranks}


##################################################################
# Now we pull all sequences from ncbi that are Oomycete and have description with ITS1
##################################################################


Entrez.email = "tkarasov@gmail.com"
search_term = "txid4762[Orgn] AND ITS1"
handle = Entrez.esearch(db="nucleotide", term=search_term, idtype="acc",
                        rettype="gb",
                        retmode="text",
                        retmax=10000000,
                        usehistory="y")
# didn't work: term="txid4762[Orgn] AND ITS1[Gene]"

# handle = Entrez.efetch(db="Taxonomy", id="4672", retmode="xml")
records = Entrez.read(handle)


# create elink:


gi_list = records["IdList"]
gi_str = ",".join(gi_list)
record = Entrez.read(Entrez.elink(id=gi_list, dbfrom="nuccore", db="nuccore"))
# Entrez.read(Entrez.elink(id=gi_list, dbfrom="genome", db="nuccore"))
# gi_conv = ','.join([link["LinkSetDb"][0]["Link"][0]["Id"] for link in record])


# The efeth algorithm needs a concatenated string
gi_conv = gi_str

# handle2 = Entrez.efetch(db="nuccore", id=gi_conv, rettype="gb",
#                         retmode="text",
#                         usehistory="y")

# *****It is important that the files be in xml format to allow for subsetting
handle3 = Entrez.efetch(db="nuccore", id=gi_conv, rettype="gb",
                        retmode="xml",
                        usehistory="y")


records3 = Entrez.read(handle3)


# Now rename records with taxonomical information
# converted = Entrez.read(Entrez.elink(id=gi_conv, dbfrom="nuccore", db="nuccore"))
# converted2 = Entrez.read(Entrez.elink(id=gi_str, dbfrom="nucleotide", db="taxonomy"))
# fin_handle = open("/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/oomycete_ITS1.fasta", "wb")
# records2 = list(SeqIO.parse(handle2, "genbank"))
# SeqIO.write(records2, fin_handle, "fasta")
# fin_handle.close()


desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

all_seqs = []
seqs_dict = {}
i = 0
for record in records3:
    seq = record['GBSeq_sequence']
    source = record['GBSeq_source']
    accession = record['GBSeq_primary-accession']
    tax = list(reversed(record['GBSeq_taxonomy'].split("; ")))
    print tax
    works = "NO"
    z = 0
    while works == "NO":
        curr_tax = ncbi.get_name_translator([tax[z]])
        try:
            if len(curr_tax.values()[0]) == 1:
                print curr_tax

                works = "YES"
                taxid = str(curr_tax.values()[0][0])
                ranks = get_desired_ranks(taxid, desired_ranks)
                kingdom = "k__Oomycete"
                phylum = "p__" + ranks['phylum_id']
                classp = "c__" + ranks['class_id']
                order = "o__" + ranks['order_id']
                family = "f__" + ranks['family_id']
                genus = "g__" + ranks['genus_id']
                species = "s__" + ranks['species_id']
                full_lineage = (";").join([kingdom, phylum, classp, order, family, genus, species])
                build_str = source + "|" + accession + "|" + full_lineage
                my_seqs = SeqRecord(Seq(seq, generic_dna), id=("strain_" + str(i)))
                seqs_dict["strain_" + str(i)] = build_str
                my_seqs.description = ''
                all_seqs.append(my_seqs)
                print curr_tax
        except IndexError:
            pass

        z = z + 1
    i = i + 1

fin_handle = open("/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/oomycete_all_S.fasta", "wb")

SeqIO.write(all_seqs, fin_handle, "fasta")

fin_handle.close()

# NOW run ITSx

os.system("/ebio/abt6_projects9/pathodopsis_microbiomes/Programs/ITSx_1.1.2/ITSx -i \
	/ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/oomycete_all_S.fasta \
	-o /ebio/abt6_projects9/pathodopsis_microbiomes/data/taxonomical_database/oomycete_real_ITS1.fasta --save_regions ITS1")

##################################################################
# The output from ITSx isn't yet compatible. Rename with the seqs_dict
##################################################################



# for linksetdb in record[0]["LinkSetDb"]:
# 	print(linksetdb["DbTo"], linksetdb["LinkName"], len(linksetdb["Link"]))

# >Clavaria_zollingeri|KP257141|SH1547001.08FU|reps|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Clavariaceae;g__Clavaria;s__Clavaria_zollingeri
#


# out_handle = open(filename, "w")
# out_handle.write(net_handle.read())
# out_handle.close()
# search_handle = Entrez.esearch(db="nucleotide", term="Opuntia[orgn] and rpl16",
#                                ... usehistory="y", idtype="acc")


# GENOME =$(esearch - db genome - query txid${TAX}[Organism:exp] |
#           efetch - format docsum | tee "${TAX}.genome.esearch.docsum")
#     ACC =`echo $GENOME | xtract - pattern DocumentSummary - element Assembly_Accession`
#     NAME =`echo $GENOME | xtract - pattern DocumentSummary - element Assembly_Name`
#     echo authoritative genome: $ACC $NAME
# records = open('\gene_id_list.txt').readlines()  # Rediret to the .txt directory
