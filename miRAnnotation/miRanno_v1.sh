#!/bin/bash

#makeblastdb -in Plant_miRNAs_v20.tsv.fa -out mirbase_plantall_miRs -dbtype nucl -title mirBASEplants

## Input file
inQuery="AsparagusNewLibsRelaxed21_22notInFInalv1.0.fa"
## Output table
outTable="BLASTAsparagusNewLibsRelaxed21_22notInFInalv1.0.txt"
## Output Alignment
outAlign="BLASTAsparagusNewLibsRelaxed21_22notInFInalv1.0.align"


## Tabular view
blastn -query $inQuery -strand plus -task blastn-short -db ../indexes/mirbase_plantall_miRs -out $outTable -perc_identity 75 -num_alignments 5 -no_greedy -ungapped -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
# ## Alignment view
blastn -query $inQuery -strand plus -task blastn-short -db ../indexes/mirbase_plantall_miRs -out $outAlign -perc_identity 75 -num_alignments 5 -no_greedy -ungapped -outfmt 0
echo "BLAST finished"
## Process annotations and prepare final results
# python3 $outTable
# echo "Annotation finished"