#!/bin/bash

java -jar /home/kakrana/tools/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 60 Lilium_4mm_ana_L1_1.fastq Lilium_4mm_ana_L1_2.fastq Lilium_4mm_ana_L1_1.pair_1.trimmed.fastq Lilium_4mm_ana_L1_1.unpair_1.trimmed.fastq Lilium_4mm_ana_L1_2.pair_2.trimmed.fastq Lilium_4mm_ana_L1_2.unpair_2.trimmed.fastq ILLUMINACLIP:/data1/homes/kakrana/tools/Trimmomatic-0.32/adapters/TruSeq-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > allout.txt 2>&1
