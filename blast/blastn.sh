#!/bin/bash

## Shell scrip to run BLAST against database of annotation from miRBase, PArse the results and annotate based on parameeters used in Amborella paper
## Required components - 1) input query in FASTA 2) Output file names 3) BLAST DB of miRNAS 4) companion python script to process annotation - miRanno.py

#makeblastdb -in Plant_miRNAs_v20.tsv.fa -out mirbase_plantall_miRs -dbtype nucl -title mirBASEplants

## Input file
inQuery="DCL_mRNA_query.fa"
## Output table
outTable="DCL_mRNA_query.txt"
## Output Alignment
outAlign="DCL_mRNA_query.align"
## miRBASE/personal sequences to be used as reference - format is BLAST DB
db="/home/kakrana/98.BLASTdb/aspa.v2/aspa.v2.genome"
## Python3 path - May differ on different systems - use command "which python3" to get path
python3path="/usr/local/bin/python3"

## Tabular view
blastn -query $inQuery -db $db -out $outTable -perc_identity 65 -num_alignments 50 -best_hit_overhang 0.15 -best_hit_score_edge 0.15 -num_threads 10 -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
# ## Alignment view
blastn -query $inQuery -db $db -out $outAlign -perc_identity 65 -num_alignments 50 -best_hit_overhang 0.15 -best_hit_score_edge 0.15 -num_threads 10 -outfmt 0
echo "BLAST finished"

# Process annotations and prepare final results
$python3path miRanno_v2.py $outTable
echo "Annotation finished"

## Version2