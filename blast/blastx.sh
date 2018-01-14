#!/bin/bash

## Shell scrip to run BLAST against database of annotation from miRBase, PArse the results and annotate based on parameeters used in Amborella paper
## Required components - 1) input query in FASTA 2) Output file names 3) BLAST DB of miRNAS 4) companion python script to process annotation - miRanno.py

# makeblastdb -in /data2/homes/kakrana/0.Atul/5.Maize_ncRNA/index/regger_ms8_peptides.fa -out /data2/homes/kakrana/0.Atul/5.Maize_ncRNA/index/regger_ms8_peptides -dbtype prot -title peptidesStanford

## Input file
inQuery="Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.fasta"
## Output table
outTable="Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.txt"
## Output Alignment
outAlign="Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.align"
## miRBASE/personal sequences to be used as reference - format is BLAST DB
db="/home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_s"
## Python3 path - May differ on different systems - use command "which python3" to get path
python3path="/usr/local/bin/python3"

## Tabular view
blastx -query $inQuery -db $db -out $outTable -best_hit_overhang 0.15 -best_hit_score_edge 0.15 -num_alignments 3 -evalue 1e-20 -num_threads 48 -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
# ## Alignment view
#blastp -query $inQuery -db $db -out $outAlign -best_hit_overhang 0.15 -best_hit_score_edge 0.15 -num_alignments 100 -outfmt 0 -html
echo "BLAST finished"

# Process annotations and prepare final results
$python3path miRanno_v2.py $outTable
echo "Annotation finished - Result file 'Final_$outTable' generated"

## Version2
