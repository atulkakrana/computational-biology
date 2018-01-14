#!/usr/bin/bash

##################################################################################
### PART-A: Generate BLAST. and other files ######################################

## Set Cores
cores=32
cd /home/kakrana/0.Atul/7.Assembl/3.aspa_trinity

### [mandatory] Identify long ORFS (>=100 AA)
echo "###############################"
echo TRANSDECODER Running
echo "###############################"
TransDecoder -S -t ./trinity_out_dir/Trinity.fasta

## [mandatory] Search Trinity transcripts
echo "###############################"
echo BLASTX Running
echo "###############################"
blastx -query ./trinity_out_dir/Trinity.fasta -db /home/kakrana/98.BLASTdb/uniprot_oct8_trinotate/uniprot_sprot.trinotate.pep -num_threads $cores -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

## [mandatory] Search Transdecoder-predicted proteins
echo "###############################"
echo BLASTP Running
echo "###############################"
blastp -query Trinity.fasta.transdecoder.pep -db /home/kakrana/98.BLASTdb/uniprot_oct8_trinotate/uniprot_sprot.trinotate.pep -num_threads $cores -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

#### [Optional] perform similar searches using uniref90 as the target database, rename output files accordingly.
## blastx -query Trinity.fasta -db /data1/homes/kakrana/98.BLASTdb/uniref_oct8_trinotate/uniprot_uniref90.trinotate.pep -num_threads 52 -max_target_seqs 1 -outfmt 6 > uniref90.blastx.outfmt6
## blastp -query Trinity.fasta.transdecoder.pep -db /data1/homes/kakrana/98.BLASTdb/uniref_oct8_trinotate/uniprot_uniref90.trinotate.pep -num_threads 52 -max_target_seqs 1 -outfmt 6 > uniref90.blastp.outfmt6

## [mandatory] HMMER to identify protein domains
echo "###############################"
echo HMMER Running
echo "###############################"
hmmscan --cpu 52 --domtblout TrinotatePFAM.out /data1/homes/kakrana/98.BLASTdb/PFAM-A_Oct_12_trinotate/Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log

## [mandatory] SIGNAL-P
echo "###############################"
echo SIGNAL-P Running
echo "###############################"
signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep

## [mandatory] TMHMM Server
echo "###############################"
echo TMHMM Running
echo "###############################"
tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out

## [mandatory] Remove rRNAs
echo "###############################"
echo RNAMMER Running
echo "###############################"
/home/kakrana/tools/Trinotate-2.0.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ./trinity_out_dir/Trinity.fasta --path_to_rnammer /usr/local/rnammer-1.2/rnammer

####################################################################################
### PART-B: Prepare SQL DB #########################################################

## 1. Retrieve the Trinotate Pre-generated Resource SQLite database - Contains both SWISSPROT and UNIREF
wget "ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz

## DB to be used for next steps
sqlite_db="Trinotate.sqlite"

## 2. [mandatory] Load transcripts and coding regions ########
echo "###############################"
echo Loading transcripts and peptides
echo "###############################"
/usr/local/trinityrnaseq-2.0.6/util/support_scripts/get_Trinity_gene_to_trans_map.pl ./trinity_out_dir/Trinity.fasta >  Trinity.fasta.gene_trans_map

## 2.B [mandatory] Initiate SQLite DB
Trinotate $sqlite_db init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta ./trinity_out_dir/Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep

## 3. [mandatory] Loading BLAST homologies ####################
## [mandatory] load protein hits
echo "###############################"
echo Loading protein set
echo "###############################"
Trinotate $sqlite_db LOAD_swissprot_blastp blastp.outfmt6

## [mandatory] load transcript hits
echo "##############################"
echo Loading blast results
echo "##############################"
Trinotate $sqlite_db LOAD_swissprot_blastx blastx.outfmt6
## load protein hits
Trinotate $sqlite_db LOAD_trembl_blastp uniref90.blastp.outfmt6
## load transcript hits
Trinotate $sqlite_db LOAD_trembl_blastx uniref90.blastx.outfmt6

## 4. [mandatory] Load Pfam domain entries ##################
echo "#############################"
echo Loading PFAM results
echo "#############################"
Trinotate $sqlite_db LOAD_pfam TrinotatePFAM.out

## 5. [mandatory] Load transmembrane domains ################
echo "############################"
echo Loading TMHMM results
echo "############################"
Trinotate $sqlite_db LOAD_tmhmm tmhmm.out

## 6. [mandatory] Load signal peptide predictions ##########
echo "###########################"
echo Loading RNAMMER results
echo "###########################"
Trinotate $sqlite_db LOAD_signalp signalp.out

###################################################################################
### [mandatory] PART-3 Trinotate: Output an Annotation Report #################################
Trinotate $sqlite_db report > trinotate_annotation_report.xls