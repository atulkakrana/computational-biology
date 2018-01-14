#!/bin/bash

## Make BLAST DB from 'nr' file - See below how to get it
# ./makeblastdb -dbtype prot -in /path/to/yourfilewithproteinsequence.fasta -parse_seqids -out myformattedDBname


## BLAST your seqeunces for input results to B2GO
./blastx -db ~/home/kakrana/tools/cpc-0.9-r2/data/prot_db -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads 5 -out local_blast.xml -query myDNAsequenceTOblast.fasta


## Optional - Get BLAST files from NCBI
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz 

## Optional - Extract using gunzip
gunzip nr.gz

#### Tutorial from BLAST2GO

# Example:
# Please have a look on the following example as there is the need to use the “sed” command line in linux in order to have the fasta file in the desired format.
# Download this fasta Viridiplantae from Uniprot.

# Open a terminal window and see how the fasta file looks like.
#  head uniprotkb_viridiplantae.fasta
#  >TR:A0A022_9ACTO A0A022 Putative dehydrogenase OS=Streptomyces ghanaensis PE=4 SV=1
#  MPSMLDAVVVGAGPNGLTAAVELARRGFSVALFEARDTVGGGARTEELTLPGFRHDPCSA
# Note: Have a look at the first line of the fasta file. The accession ID A0A022 is not in between “|”. There is the need to reformat the whole fasta for the accessions IDs.

# Rectify the fasta file in order to have the correct format.
# sed -E 's/(>[A-Z0-9:_]+) ([A-Z0-9]+) (.*)/\1 gnl|\2| \3/g' uniprotkb_viridiplantae.fasta > uniprotkb_viridiplantae_mod.fasta 
# Note: This is an example and not an universal command. The user will need to understand from their sequences how to change them in order to obtain the correct format.

# Lets have a look at the modified file.
#  head uniprotkb_viridiplantae_mod.fasta
#  >TR:A0A022_9ACTO|A0A022| Putative dehydrogenase OS=Streptomyces ghanaensis PE=4 SV=1
#  MPSMLDAVVVGAGPNGLTAAVELARRGFSVALFEARDTVGGGARTEELTLPGFRHDPCSA
# Now it looks very similar to what we wanted.

# It is now safe to create the database using the Blast+ executables.
# ./makeblastdb -dbtype prot -in /path/to/uniprotkb_viridiplantae_mod.fasta -parse_seqids -out uniprotkb_viridiplantae_mod_db

# Run blast.
# ./blastx -db ~/path/to/your/uniprotkb_viridiplantae_mod_db/ -outfmt 5 -evalue 1e-3 -word_size 3 -show_gis -num_alignments 20 -max_hsps 20 -num_threads 5 -out local_blast.xml -query 10_seq.fasta

# Load your local_blast.xml file into Blast2GO (File -> Load -> Load Blast Results -> XML files) and visualize several Blast results (Show Blast Results) to see if the accession appears in the right place (ACC).

# You can proceed with the mapping step as usual.




