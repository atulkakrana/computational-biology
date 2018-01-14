#!/usr/bin/bash


####### Computing N50 ########################################################################################
####### https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats ##########
# /usr/local/trinityrnaseq-2.0.6/util/TrinityStats.pl Trinity.fasta > assembly_summary.txt


####### Full-length transcript analysis for model and non-model organisms using BLAST+ ########################
####### https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts ##########

##  Build a blastable database
# makeblastdb -in uniprot_sprot.fasta -dbtype prot

# ## Perform the blast search, reporting only the top alignment:
# cd /home/kakrana/0.Atul/7.Assembl/2.daylily_trinity
# blastx -query ./trinity_out_dir/Trinity.fasta -db /home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_s -out ./blastx.outfmt6 -evalue 1e-20 -num_threads 48 -max_target_seqs 1 -outfmt 6

## Examine the percent of the target being aligned to by the best matching Trinity transcript
# cd /home/kakrana/0.Atul/7.Assembl/2.daylily_trinity
# /usr/local/trinityrnaseq-2.0.6/util/analyze_blastPlus_topHit_coverage.pl ./blastx.outfmt6 ./trinity_out_dir/Trinity.fasta /home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_sprot.fasta > blast_summary.txt

# ## Group the multiple HSPs per transcript/database_match pairing like so - Scripts missing in trinity package
cd /home/kakrana/0.Atul/7.Assembl/2.daylily_trinity
/usr/local/trinityrnaseq-2.0.6/util/misc/blast_outfmt6_group_segments.pl ./blastx.outfmt6  ./trinity_out_dir/Trinity.fasta  /home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_sprot.fasta > blast.outfmt6.grouped

# ## Then compute the percent coverage by length histogram like so.
cd /home/kakrana/0.Atul/7.Assembl/2.daylily_trinity
/usr/local/trinityrnaseq-2.0.6/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl ./blast.outfmt6.grouped > blast_summary2.txt


######### Computing N50X ######################################################################################
######### https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats - Requires quantification step first ########
# /usr/local/trinityrnaseq-2.0.6/util/misc/contig_E_statistic.pl transcripts.TMM.EXPR.matrix ./trinity_out_dir/Trinity.fasta | tee ExN50.stats