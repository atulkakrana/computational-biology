#!/bin/bash

## code to combine
## counts for 10x
## for different 
## libraries; the
## output folder 
## from "counts.sh"
## serve as inputs for this one

## Notes
## 1. change "--id" to any string/text that you want to name folder for output
## 2. change "--csv" to path for CSV file with sample_id, molecule_h5 file path ,batch , group

## Aggregate Counts For Samples - E16.5 and P0
/home/anand/tools/cellranger-6.0.2/bin/cellranger aggr --id=run_aggr_E16_5_P0 \
--csv=/home/anand/0.work/1.scseq/aggr.csv