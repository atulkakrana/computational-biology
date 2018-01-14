#!/usr/bin/bash

## SAM TO BAM
## bam extenstion will be added automatically
samtools view -bS FILE.sam | samtools sort - FILE.sorted

## BAM to sorted BAM
samtools index FILE.sorted.bam FILE.mod.sorted.bai