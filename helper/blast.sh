#!/bin/bash

##MakeDB
# makeblastdb -in Asparagus.meta.unigenes.fasta -out AspaMeta -dbtype nucl -title AspaMeta

## Short BLAST
#blastn -query RealaxNotStringent_20_22.csv.fa -strand plus -task blastn-short -db ./indexes/mirbase_plantall_miRs -out BLAST_RealaxNotStringent_20_22.txt -perc_identity 75 -num_alignments 5 -no_greedy -ungapped -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore"

## EPRV Example
#blastn -db /data2/homes/kakrana/database/nt/at_genomic_db -query RTV_BSV.fa -perc_identity 80 -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -num_threads 6 -out RTV_BSV_BLAST_best_hit_pid80

## tBLASTn
# tblastn -query AT_OS_ZMA_AGO_Clean.fa -db ./indexes/AspaMeta -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -num_alignments 5 -num_threads 24 -outfmt "6 qseqid sseqid pident ppos gaps nident positive length qlen qstart qend sstart send sstrand evalue bitscore" -out BLAST_AspaMeta_AT_OS_ZMA_AGO_Clean.tsv

##blastn
blastn -query BLAST_Meta_IDS.fa -db ./indexes/AspaTrinity -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -perc_identity 70 -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" -num_threads 36 -out BLAST_BLAST_Meta_IDS_AspaTrinity.tsv