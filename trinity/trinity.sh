#!/usr/bin/bash


############# 2. TRINITY ASSEMBLY#########################

#### Assembly w/o genome #################################
# Trinity --seqType fq --left left_SE.fq --right right.fq --SS_lib_type RF --min_contig_length 150 --CPU 10 --max_memory 500G > run.log 2>&1


#### Assembly with genome ################################
#### Map reads to genome using GSNAP and provide a sorted BAM file

# Make GMAP DB
# cd /home/kakrana/gmap_db/
# gmap_build -D ~/gmap_db/ -d Streptochaeta.v1.0 Streptochaeta_MaSuRCA_scaffolds_gt871_edit3.fasta

# # Run GMAP to map polished HQ transcripts to genome
# cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
# gmap -D ~/gmap_db -d Asparagus.v2.0 -f samse -t 24 -n 0 ./Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.fasta > ./Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sam 2> ./Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sam.log

# # # Sort\
# cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
# sort -k 3,3 -k 4,4n Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sam  > Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sorted.sami

# Trinity --seqType fq --left left.fq --right right.fq --SS_lib_type RF --min_contig_length 150 --genome_guided_bam rnaseq.coordSorted.bam --genome_guided_max_intron 10000 --CPU 10 --max_memory 500G > run.log 2>&1

#### END OF ASSEMBLY ######################################


########### 2. QUALITY ASSESMENT ##########################

#### Computing N50 ######################################  
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats
cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
/usr/local/trinityrnaseq-2.0.6/util/TrinityStats.pl ./trinity_out_dir/Trinity.fasta > assembly_summary.txt

#### Assessing the Read Content #########################
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly
## For strand-specific RNA-Seq data use '--SS_lib_type' parameter, and put this parameter before the '--' above, since all the parameters after '--' are applied to the bowtie aligner.
# /usr/local/trinityrnaseq-2.0.6/util/bowtie_PE_separate_then_join.pl --seqType fq --left left.fq --right right.fq --SS_lib_type RF --target ./trinity_out_dir/Trinity.fasta --aligner bowtie -- -p 52 --all --best --strata -m 300 

#### Full-length transcript analysis for model and non-model organisms using BLAST+ ########################
#### https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts ##########

##  Build a blastable database
# makeblastdb -in uniprot_sprot.fasta -dbtype prot

## Perform the blast search, reporting only the top alignment:
#cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
#blastx -query ./trinity_out_dir/Trinity.fasta -db /home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_s -out ./blastx.outfmt6 -evalue 1e-20 -num_threads 48 -max_target_seqs 1 -outfmt 6

## Examine the percent of the target being aligned to by the best matching Trinity transcript
#cd /home/kakrana/0.Atul/7.Assembl/4.strep_trinity
#/usr/local/trinityrnaseq-2.0.6/util/analyze_blastPlus_topHit_coverage.pl ./blastx.outfmt6 ./trinity_out_dir/Trinity.fasta /home/kakrana/98.BLASTdb/swissprot_Sep30/uniprot_sprot.fasta > blast_summary.txt

#########################################################
