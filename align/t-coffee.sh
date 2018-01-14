#!/usr/bin/bash

t_coffee Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.fa -n_core 48 -mode accurate -output=clustalw,fasta_aln,msf

#### Sequence conversion #################
##########################################
## Manual https://media.readthedocs.org/pdf/seqmagick/latest/seqmagick.pdf
seqmagick convert --output-format phylip-relaxed --alphabet protein Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.aln Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.phyx

##### PHylogentic analysis by PhyML ######
##########################################
PhyML -i Os_Zm_Gm_At_Ao_DL_La_Atr_Zmr_AGOs.all.finalnames.phyx --datatype aa --search BEST


# Sequence format:
# - clustalw_aln, clustalw
# - gcg
# - msf_aln
# - pir_aln, pir_seq
# - fasta_aln, fasta_seq
# - phylip

# Output format:
# - score_ascii : causes the output of a reliability flag
# - score_html  : causes the output to be a reliability plot in HTML
# - score_pdf   : idem in PDF (if ps2pdf is installed on your system)
# - score_ps    : idem in postscript

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