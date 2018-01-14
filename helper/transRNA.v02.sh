#!/bin/bash

## Python3 path - May differ on different systems - use command "which python3" to get path
python3path="/usr/local/bin/python3"

## Required components - 1) input query in FASTA 2) companion python script to process annotation - miRanno.py 3) cleanFASTA.v.1.3 or above

# ## Input file


### MANNUAL STEPS #####################################

## Clean FASTA file - Use cleanFasta.v02.py or above
## Complement your transcripts FASTA file - Use cleanFasta.v1.03.py or above [mode 3]
## Reverse complement transcripts file - use cleanFasta [mode 2]

# ## Make a BLAST DB with normal (clean) FASTA
# makeblastdb -in 24phas.fa -out BLAST_NOR -dbtype nucl -title BLAST_NOR

# ## Make a BLAST DB with revcomp FASTA
# $python3 cleanFasta.v1.3.py $infile 2
# makeblastdb -in 24phas.clean.revcomp.fa -out BLAST_RC -dbtype nucl -title BLAST_RC 
########################################################


#### AUTOMATIC STEPS ###################################
## Main query
inQuery="24phas.fa" ## Main transcripts file to identify IRs (with BLAST_RC) and Isoforms (with BLAST_Norm)
compQuery="24phas.clean.comp.fa" ## Complemented transcripts to identify transcripts from other stand

## Required BLAST DBs
# cd /home/kakrana/0.Atul/1.PHAS/25.Lilium/2.24PHAS/analysis/IRphasi/analysis
$BLAST_RC="./BLAST_RC" ### There must be a BLAST DB made from reverse complement of main transcripts, this will be used to find IRs
$BLAST_NOR="./BLAST_NOR" ### There must be a DB of this name made from normal transcipts, this will be used to identify overlappings isofolms or transcripts in other strands

## BLAST to check for trans NAT RNAs
blastn -query $inQuery -strand plus -db BLAST_RC -out BLAST_RC.txt -perc_identity 75 -num_alignments 5 -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
## BLAST to check isoforms
blastn -query $inQuery -strand plus -db BLAST_NOR -out BLAST_NOR.txt -perc_identity 75 -num_alignments 5 -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
# ## BLAST to check overlapping transcripts from other strands
blastn -query $compQuery -strand plus -db BLAST_NOR -out BLAST_COMP.txt -perc_identity 75 -num_alignments 5 -outfmt "6 qseqid query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore"
echo "BLAST finished"

# # Process annotations and prepare final results
$python3path miRanno_v2.py BLAST_RC.txt
$python3path miRanno_v2.py BLAST_NOR.txt
$python3path miRanno_v2.py BLAST_COMP.txt
echo "Annotation finished"

### RUN FINAL get IR python script 

## Version-1