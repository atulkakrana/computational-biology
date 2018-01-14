#PBS -N GMAP -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l walltime=96:00:00
#PBS -m ae  Mail status when the job completes
#PBS -M kakrana@udel.edu  Mail to this address

# Enter SMRT tools [As e-mailed by karol]
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
sort -k 3,3 -k 4,4n Trinity.collapsed.rep.sam  > Trinity.collapsed.rep.sam 

# #### GET GFF FORMAT ################### 
# #######################################
# gmap -D ~/gmap_db -d Zea_mays.AGP.mod.v2.17 -f 2 -t 8 -n 0 Trinity.collapsed > ./Trinity.gff3 2> ./Trinity.gff.log

# #### SAM TO SORTED BAM ################
# #######################################
# samtools view -bS Trinity.collapsed.rep.sam | samtools sort - Trinity.collapsed.rep.sorted
# samtools index Trinity.collapsed.rep.sorted.bam Trinity.collapsed.rep.sorted.bai



## V1- > v2
## Cleaner code


## USE
## Use for gmap operarions - creating DB, Mmapping transcripts to get SAM, sorted SAM, BAM and BAM index

