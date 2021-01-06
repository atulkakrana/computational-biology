## Rocket v4b1
Note: the data files are deleted during the run to save disk space, so, please keep a copy of your FASTQ files elsewhere for longterm storage (on S3 or other device)


### Data Prepration

**Paired End Data**
* Rename the paired (R1 and R2) FASTQ file by adding suffixes - "_1" and "_2".
    For Example: Paired-end samples with names "AB_S7_L001_R1_001.fastq" and "AB_S7_L001_R2_001.fastq" can be renamed to "AB_S7_L001_1.fastq" and "AB_S7_L001_2.fastq"

* When providing names in the `sampleinfo` file, exclude the "_1" and "_2" suffixes.
    For Example: Paired end samples renamed above - "AB_S7_L001_1.fastq" and "AB_S7_L001_2.fastq" can be added to `sampleInfo` file as "AB_S7_L001" and ""AB_S7_L001"

**Data Preprocessing**
* NOTE: From v4, to save disk space, we delete input FASTQ files after they are preprocessed. Make sure to store your libraries elsewhere for future use. 

### Future Updates
* GTF files processing  - the current modified GTFs may not be best options (especially for lncRNAs)
* Genome assembly GRC39 - Update genome, corresponding indexes to grc39 (for mouse and human)

### AWS Instances
`t2.micro` for basic operations      
`m5a.2xlarge` for **testing** on 6 paired end libraries
