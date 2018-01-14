#!/usr/bin/bash


### TRINITY BASED ABUNDANCE ESTIMATIOM ##############################
#####################################################################

## This method to be used whn genome for species in not available
## You must have a Trinity assembly done before using this
## Link: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification

## Just prepare the reference for alignment and abundance estimation
## For DEGs this step needs to be done on every library and compiled as matrix later (below)
# echo "###################################################"
# echo "Preparing trinity assembly for abundance estimation"
# echo "###################################################"
# /usr/local/trinityrnaseq-2.0.6/util/align_and_estimate_abundance.pl --transcripts ./trinity_out_dir/Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --thread_count 64
# echo "Prepration done"

# ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.) - 
# echo "###################################################"
# echo "Performing abundance estimation"
# echo "###################################################"
# /usr/local/trinityrnaseq-2.0.6/util/align_and_estimate_abundance.pl --transcripts ./trinity_out_dir/Trinity.fasta --seqType fq --left left_SE.fq --right right.fq --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --thread_count 64
# echo "Estimation done"


### Make matrix of all samples for DEG analysis in edgeR or so
### Needs computation of abundances for each sample in a seprate folder see above command or dashi.dbi.udel.edu/data1/homes/kakrana/0.Atul/7.Assembl/4.strep_trinity/raw/Strep_an_30
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#build-transcript-and-gene-expression-matrices
cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
/usr/local/trinityrnaseq-2.0.6/util/abundance_estimates_to_matrix.pl --est_method RSEM --name_sample_by_basedir ./raw/Strep_an_1_5r/RSEM.genes.results ./raw/Strep_an_3_0r/RSEM.genes.results ./raw/Strep_an_15/RSEM.genes.results ./raw/Strep_an_30/RSEM.genes.results