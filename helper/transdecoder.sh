#!/bin/bash

### Change Directory ###
echo "Changing directory to:"
pwd
cd /home/kakrana/8.Collab/4.Kak/5.lncrna/8.wbdata/WB_RNAseq/transdecoder

# ### 1. Extract CDS from cufflinks merged assembly ####
# echo "Extracting CDS"
# /home/kakrana/tools/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl ../final_merged.gtf /home/kakrana/8.Collab/4.Kak/99.genome/GenomeM.dna_sm.toplevel.fa > final_merged.cds.fa
# echo "Done"

# ### 2. Convert GTF to GFF ####
# echo "Converting to GFF"
# /home/kakrana/tools/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl ../final_merged.gtf > ./final_merged.gff3
# echo "Done"

# #### 3. Generate best candidate ORF predictions #####
# echo "Predicting CDS"
# TransDecoder.LongOrfs -t ./final_merged.cds.fa
# ### optionally, identify peptides with homology to known proteins)
# TransDecoder.Predict -t ./final_merged.cds.fa
# echo "Done"

#### 4. Generate a genome-based coding region annotation file ####
echo "Generating final annotations"
/home/kakrana/tools/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl ./final_merged.cds.fa.transdecoder.gff3 final_merged.gff3 final_merged.cds.fa > final_merged.cds.fa.transdecoder.genome.gff3
echo "Done"
