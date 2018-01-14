#!/usr/bin/bash

# Make GMAP DB
# cd /home/kakrana/gmap_db/
# gmap_build -D ~/gmap_db/ -d Asparagus.v2.0 Asparagus.V2.0.genome.stable_modifed.fa

# # Run GMAP to map polished HQ transcripts to genome
cd /home/kakrana/0.Atul/7.Assembl/3.aspa_trinity
## Get SAM file
gmap -D ~/gmap_db -d Asparagus.v2.0 -f samse -t 52 -n 0 ./trinity_out_dir/Trinity.fasta > ./Trinity.sam 2> ./Trinity.sam.log
# ## Get GFF format
# gmap -D ~/gmap_db -d Asparagus.v2.0 -f 2 -t 52 -n 0 ./trinity_out_dir/Trinity.fasta > ./Trinity.gff3 2> ./Trinity.sam.log

## Convert SAM to sorted BAM and make index
cd /home/kakrana/0.Atul/7.Assembl/3.aspa_trinity
samtools view -bS Trinity.sam | samtools sort - Trinity.sorted
samtools index Trinity.sorted.bam Trinity.sorted.bai

### ALternative sorting options
# sort -k 3,3 -k 4,4n Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sam  > Aspar_1-3all_ref_polished_high_qv_consensus_isoforms.mod.sorted.sam