from Bio import SeqIO
import sys
import gzip
from re import match

# usage python demultiplex_16S_non_flashed.sh output_filename R1.fq R2.fq

file_name = sys.argv[1]
outfile_path = sys.argv[2]
fastq_F_filename = sys.argv[3]
fastq_R_filename = sys.argv[4]


def my_filter(records_forward, records_reverse, file_name):
    '''this assumes that the forward and reverse reads are in the same order'''
    all_files_R1 = {}
    all_files_R2 = {}
    for rec1, rec2 in zip(records_forward, records_reverse):
        match_forward = [rec for rec in all_forward if str(rec1.seq).startswith(rec)]
        match_reverse = [rec for rec in all_reverse if str(rec2.seq).startswith(rec)]
        if len(match_forward) != 1 or len(match_reverse) != 1:
            continue
        forward_match = forward_dict[match_forward[0]]
        reverse_match = reverse_dict[match_reverse[0]]
        try:
            #tempF= all_files_R1[forward_match+"_"+reverse_match]
            #print("tempF is:"+str(len(tempF)))
            # tempF2=tempF+[rec1]
            all_files_R1[forward_match + "_" + reverse_match].append(rec1)
            #tempR= all_files_R2[forward_match+"_"+reverse_match]
            # tempR2=[tempR]+[rec2]
            all_files_R2[forward_match + "_" + reverse_match].append(rec2)
            print(len(all_files_R1[forward_match + "_" + reverse_match]))
        except KeyError:
            all_files_R1[forward_match + "_" + reverse_match] = []
            all_files_R1[forward_match + "_" + reverse_match].append(rec1)
            all_files_R2[forward_match + "_" + reverse_match] = []
            all_files_R2[forward_match + "_" + reverse_match].append(rec2)

    # now write out the parsed files
    for file in list(all_files_R1.keys()):
        SeqIO.write(all_files_R1[file], outfile_path + "/demult_python/" + file + "_R1_" + file_name + ".fastq", "fastq")
    for file in list(all_files_R2.keys()):
        SeqIO.write(all_files_R2[file], outfile_path + "/demult_python/" + file + "_R2_" + file_name + ".fastq", "fastq")


# path="/ebio/abt6_projects9/pathodopsis_microbiomes/data/raw_reads/16S/16S_12_2018/"

forward_F1 = "GACCTACGGGAGGCAGCAG"
forward_F2 = "TGACCTACGGGAGGCAGCAG"
forward_F3 = "CTGACCTACGGGAGGCAGCAG"
forward_F4 = "ACTGACCTACGGGAGGCAGCAG"
forward_F5 = "GACTGACCTACGGGAGGCAGCAG"
forward_F6 = "TGACTGACCTACGGGAGGCAGCAG"
reverse_F1 = "ACGGACTAC"  # ??GGGT?TCTAAT
#"ATTAGA?ACCC??GTAGTCCGT"
reverse_F2 = "TACGGACTAC"
#"ATTAGA?ACCC??GTAGTCCGTA"
reverse_F3 = "CTACGGACTAC"
#"ATTAGA?ACCC??GTAGTCCGTAG"
reverse_F4 = "ACTACGGACTAC"
#"ATTAGA?ACCC??GTAGTCCGTAGT"
reverse_F5 = "GACTACGGACTAC"
#"ATTAGA?ACCC??GTAGTCCGTAGTC"
reverse_F6 = "TGACTACGGACTAC"
#"ATTAGA?ACCC??GTAGTCCGTAGTCA"

all_forward = [forward_F1, forward_F2, forward_F3, forward_F4, forward_F5, forward_F6]
all_reverse = [reverse_F1, reverse_F2, reverse_F3, reverse_F4, reverse_F5, reverse_F6]

forward_dict = {}
reverse_dict = {}
for i in range(1, 7):
    forward_dict[all_forward[i - 1]] = "341F" + str(i)
    reverse_dict[all_reverse[i - 1]] = "806F" + str(i)


records_forward = list(SeqIO.parse(gzip.open(fastq_F_filename, 'rt'), "fastq"))
records_reverse = list(SeqIO.parse(gzip.open(fastq_R_filename, 'rt'), "fastq"))

# now do demultiplexing
my_filter(records_forward, records_reverse, file_name)
