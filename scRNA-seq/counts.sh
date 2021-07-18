#!/bin/bash

## code to compute
## counts for 10x
## FASTQ Files

## Notes
## 1. change "--id" to any string/text that you want to name folder for output
## 2. change "--sample" to exact sample prefix as FASTQ files

## Generate Counts For Samples - E16-5-lens
/home/anand/tools/cellranger-6.0.2/bin/cellranger count --id=run_count_E16_5_lens \
--fastqs=/home/anand/0.work/1.scseq/data \
--sample=E16-5-lens \
--transcriptome=/home/anand/references/refdata-gex-mm10-2020-A


## Generate Counts For Samples - P0
/home/anand/tools/cellranger-6.0.2/bin/cellranger count --id=run_count_P0 \
--fastqs=/home/anand/0.work/1.scseq/data \
--sample=P0 \
--transcriptome=/home/anand/references/refdata-gex-mm10-2020-A