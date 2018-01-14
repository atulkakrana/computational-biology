#!/usr/bin/bash


### TRINITY BASED ABUNDANCE ESTIMATIOM ##############################
#####################################################################

## Declare working directory
cd /data1/homes/kakrana/0.Atul/7.Assembl/3.aspa_trinity/RSEM

### 1. Compute Lib-wise abundances
### See abundances.v2.sh - This must be run on each library under their individual folder, these folder names are used as input in step-2

### 2. Make matrix of all samples for DEG analysis in edgeR or so
### Needs computation of abundances for each sample in a seprate folder see above command or dashi.dbi.udel.edu/data1/homes/kakrana/0.Atul/7.Assembl/4.strep_trinity/raw/Strep_an_30
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#build-transcript-and-gene-expression-matrices

# /usr/local/trinityrnaseq-2.0.6/util/abundance_estimates_to_matrix.pl --est_method RSEM --name_sample_by_basedir Asp_0_5_ant_b/RSEM.genes.results Asp_0_5_ant_budr/RSEM.genes.results Asp_1_ant_b/RSEM.genes.results Asp_1_ant_budr/RSEM.genes.results Asp_le/RSEM.genes.results  Asp_leaf/RSEM.genes.results

### 3. Counting Numbers of Expressed Genes (for transcripts see above link) ###
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#counting-numbers-of-expressed-transcripts-or-genes
## Script name is different then on webpage -tested OK
# /usr/local/trinityrnaseq-2.0.6/util/misc/count_matrix_features_given_MIN_FPKM_threshold.pl matrix.not_cross_norm.fpkm.tmp | tee trans_matrix.fpkm.not_cross_norm.counts_by_min_fpkm

### 4. Compute ExN50 ###################################
# /usr/local/trinityrnaseq-2.0.6/util/misc/contig_E_statistic.pl matrix.TMM.fpkm.matrix ../trinity_out_dir/Trinity.fasta | tee ExN50.stats