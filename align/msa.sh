#!/usr/bin/bash

## This script is used to geenrate MSA and perform phylogenetic analysis

input="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalcands.fa"
output="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalcands.aln"
output2="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalcands.phyi"


##### MSA BY CLUSTAL OMEGA #######
##################################
#### The output mighy be incompatibel with PhyML but should work with PhyLip

## Clustal output for other downstream analysis
# clustalo -i $input --infmt=fa -o $output --outfmt=clu -iter=1000 -v --auto --threads=6
## Phylip output for PhyML
# clustalo -i $input --infmt=fa -o $output2 --outfmt=clu -iter=1000 -v --auto --threads=6

#### MSA BY MUSCLE3 #############
#################################
# Output in clustal for basic tools and phylip (interleaved format) for PhyMl
# muscle3 -in $input -clwout $output -phyiout $output2

##### CONVERT FORMAT #################
######################################
### Clustalw for mat needs to be changed to NEXUS format
#### Use Online converter: http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_nexus.php
#### And remove 'missing=?' and "gap=-" from header of NEXUS file as these are not accepted by PhyML
# seqmagick convert --output-format nexus --alphabet dna input.fasta output.nex

##### PHylogentic analysis by PhyML ######
##########################################
## PhyML phlogenetic analysis
# PhyML -i $output2 --datatype aa --search BEST --bootstrap 10


mpirun -n 48 phyml-mpi -i $output2 --datatype aa --search BEST --bootstrap 10