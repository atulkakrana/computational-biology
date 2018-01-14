#!/usr/bin/bash


### TRINITY BASED ABUNDANCE ESTIMATIOM ##############################
#####################################################################

## This method to be used whn genome for species in not available
## You must have a Trinity assembly done before using this
## Link: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification

## Just prepare the reference for alignment and abundance estimation
echo "###################################################"
echo "Preparing trinity assembly for abundance estimation"
echo "###################################################"
/usr/local/trinityrnaseq-2.0.6/util/align_and_estimate_abundance.pl --transcripts /home/kakrana/0.Atul/7.Assembl/4.strep_trinity/trinity_out_dir/Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --thread_count 32
echo "Prepration done"

## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.) - 
echo "###################################################"
echo "Performing abundance estimation"
echo "###################################################"
/usr/local/trinityrnaseq-2.0.6/util/align_and_estimate_abundance.pl --transcripts /home/kakrana/0.Atul/7.Assembl/4.strep_trinity/trinity_out_dir/Trinity.fasta --seqType fq --left Strep_an_30.chopped.pair_1.trimmed.fastq --right Strep_an_30.chopped.pair_2.trimmed.fastq --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --thread_count 32
echo "Estimation done"
