#!/usr/bin/bash

input="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.fa"
output="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.align"
output2="Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.phyx"

mafft --genafpair --maxiterate 1000 --thread -48 $input > $output

#### Sequence conversion #################
##########################################
## Manual https://media.readthedocs.org/pdf/seqmagick/latest/seqmagick.pdf
seqmagick convert --output-format phylip-relaxed --input-format fasta --alphabet protein $output $output2

##### PHylogentic analysis by PhyML ######
##########################################
PhyML -i $output2 --datatype aa --search BEST


### Seqmagik formats
# Extension Format
# .afa fasta
# .aln clustal
# .fa fasta
# .faa fasta
# .fas fasta
# .fasta fasta
# .fastq fastq
# .ffn fasta
# .fna fasta
# .fq fastq
# .frn fasta
# .gb genbank
# .gbk genbank
# .needle emboss
# .nex nexus
# .phy phylip
# .phylip phylip
# .phyx phylip-relaxed
# .qual qual
# .sff sff-trim
# .sth stockholm
# .sto stockholm
# Note: NEXUS-format output requires the --alphabet flag