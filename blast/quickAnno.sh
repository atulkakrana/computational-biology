#!/bin/bash


## FOR CDS ###
# blastx -query 24miR.targets.genic.cds.fa -db /home/kakrana/98.BLASTdb/uniprot_oct8_trinotate/uniprot_sprot.trinotate.pep -num_threads 40 -max_target_seqs 1 -evalue 0.0001 -outfmt "6 qseqid sseqid stitle pident length mismatch qstart qend sstart send evalue bitscore" > blastx.outfmt6

### FOR PROTEINS ###
# blastp -query Asparagus.V2.0.genome.stable.protein.headfix.fa -db /home/kakrana/98.BLASTdb/uniprot_oct8_trinotate/uniprot_sprot.trinotate.pep -num_threads $cores -max_target_seqs 1 -evalue 0.0001 -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore'  > blastp.outfmt6
