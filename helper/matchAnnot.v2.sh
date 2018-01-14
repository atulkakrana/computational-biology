#!/usr/bin/bash

#PBS -N PacMap -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l walltime=96:00:00
#PBS -m ae  Mail status when the job completes
#PBS -M kakrana@udel.edu  Mail to this address

# Enter SMRT tools [As e-mailed by karol] - Do these commands manuualy
# SMRT_ROOT=/opt/smrtanalysis
# . $SMRT_ROOT/current/etc/setup.sh

# #### Make GMAP DB ######################
# ########################################
# cd /home/kakrana/gmap_db/
# gmap_build -D ~/gmap_db/ -d Asparagus.v2.0 Asparagus.V2.0.genome.stable_modifed.fa

# #### SET PATH #########################
# #######################################
cd /home/kakrana/0.Atul/5.Maize_ncRNA/4.matannot/v2.hybrid.agp.v2

# #### Get SAM file #####################
# #######################################
gmap -D ~/gmap_db -d Zea_mays.AGP.mod.v2.17 -f samse -t 48 -n 0 Trinity.collapsed.rep.fa > ./Trinity.collapsed.rep.sam 2> ./Trinity.collapsed.rep.sam.log

# #### SORT SAM FILE ####################
# #######################################
sort -k 3,3 -k 4,4n Trinity.collapsed.rep.sam  > Trinity.collapsed.rep.sorted.sam 

## Sort using Samtools - Convert to BAM, sort, reconvert to SAM
# samtools view -S -b -o Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.bam Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.sam
# samtools sort Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.bam Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.sorted
# samtools view -h Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.sorted.bam > Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv_fusion.rep.sorted.sam

## Get GTF from Plant ensembl
# wget ftp://ftp.ensemblgenomes.org/pub/release-27/plants/gtf/zea_mays/Zea_mays.AGPv3.27.gtf.gz
# gunzip Zea_mays.AGPv3.27.gtf.gz

# #### RUN MATCHANNOT ##################
# #### MatchAnnot commands - https://gist.github.com/rlleras/f31cfec1e36efd289ae0
# ######################################
# matchAnnot.py --gtf /home/kakrana/gmap_db/AGPv2/Zea_mays.AGPv2.17.mod.gtf --format alt Trinity.collapsed.rep.sorted.sam > trinity.hybird.collapsed.matchedAnnot.agpv2.out
