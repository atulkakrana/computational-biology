#!/bin/bash

## Shell scrip to run BLAST against database of annotation from miRBase, PArse the results and annotate based on parameeters used in Amborella paper
## Required components - 1) input query in FASTA 2) Output file names 3) BLAST DB of miRNAS 4) companion python script to process annotation - miRanno.py

#makeblastdb -in Plant_miRNAs_v20.tsv.fa -out mirbase_plantall_miRs -dbtype nucl -title mirBASEplants

## Input file
inQuery="AsparagusStrigent_23_24_miRNA.fa"
## Output table
outTable="BLASTAsparagusStrigent_23_24_miRNA.txt"
## Output Alignment
outAlign="BLASTAsparagusStrigent_23_24_miRNA.align"
## miRBASE/personal sequences to be used as reference - format is BLAST DB
db="./miRDB/mirbase_plantall_miRs"
## Python3 path - May differ on different systems - use command "which python3" to get path
python3path="/usr/local/bin/python3"

## Tabular view
blastn -query $inQuery -strand plus -task blastn-short -db $db -out $outTable -perc_identity 75 -num_alignments 5 -no_greedy -ungapped -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
# ## Alignment view
blastn -query $inQuery -strand plus -task blastn-short -db $db -out $outAlign -perc_identity 75 -num_alignments 5 -no_greedy -ungapped -outfmt 0
echo "BLAST finished"

# Process annotations and prepare final results
$python3path miRanno_v2.py $outTable
echo "Annotation finished - Result file 'Final_$outTable' generated"

## Version2