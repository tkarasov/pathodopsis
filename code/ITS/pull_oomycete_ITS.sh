#!/bin/sh
# https://www.biostars.org/p/244479/

set -u
#TAXLIST=("Daphnia pulex" "Drosophila melanogaster" "Anopheles gambiae" "Pediculus humanus"
#"Ixodes scapularis" "Apis mellifera"  "Bombyx mori")
#TAXLIST=("Strigamia maritima")
TAXLIST=("4762")
for TAX in "${TAXLIST[@]}" ; do
echo getting genome for: $TAX
   GENOME=$(esearch -db genome -query txid${TAX}[Organism:exp] |
       efetch -format docsum | tee "${TAX}.genome.esearch.docsum")
   ACC=`echo $GENOME | xtract -pattern DocumentSummary  -element Assembly_Accession`
   NAME=`echo $GENOME |  xtract -pattern DocumentSummary -element Assembly_Name`
   echo authoritative genome: $ACC $NAME
   for RES in "$ACC[@]" ; do
      RESULT=$(esearch -db assembly -query "$RES" |
      efetch -format docsum | tee "${TAX}.assembly.esearch.docsum")
   FTPP=`echo $RESULT | xtract -pattern DocumentSummary  -element FtpPath_GenBank`
   TAXID=`echo $RESULT | xtract -pattern DocumentSummary  -element Taxid`
   echo FtpPath: $FTPP
   BASENAME=`basename $FTPP`
   FTPPATHG=$FTPP/$BASENAME'_genomic.fna.gz'
   #FTPPATHP=$FTPP/$BASENAME'_protein.faa.gz'
   echo Downloading $FTPPATHG ...

   ## get genome data
   wget $FTPPATHG
   # BASENAME=`basename $FTPPATHG`
   # gunzip -f $BASENAME
   # BASENAME=`echo $BASENAME | sed s/.gz//`
   # makeblastdb -in $BASENAME -dbtype nucl -parse_seqids -taxid $TAXID -title "$TAX $NAME genomic"
   # echo Downloading $FTPPATHP ...
   # ## get protein data
   # wget $FTPPATHP
   # BASENAME=`basename $FTPPATHP`
   # gunzip -f $BASENAME
   # BASENAME=`echo $BASENAME | sed s/.gz//`
   # makeblastdb -in $BASENAME -dbtype prot -parse_seqids -taxid $TAXID -title "$TAX $NAME proteins"


done

### Make a blast db of all fasta files in directory
ls  *.fna > dblist_nuc.txt
blastdb_aliastool -dblist_file dblist_nuc.txt -out SelectedArthropds -title 'selection of arth genomes' -dbtype nucl
ls *.faa > dblist_prot.txt
blastdb_aliastool -dblist_file dblist_prot.txt -out SelectedArthropds -title 'selection of arth proteins' -dbtype prot