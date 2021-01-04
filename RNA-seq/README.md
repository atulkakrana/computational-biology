## Rocket v4b1

### Data Prepration
**Paired end data**
* rename the paired (R1 and R2) files to have "_1" and "_2"
    For Example: paired samples "68A_S7_L001_R1_001.fastq"  and "68A_S7_L001_R2_001.fastq" are renamed to "68A_S7_L001_1.fastq" and "68A_S7_L001_2.fastq"

* in sampleinfo file when providing names, exclude the "_1" and "_2" suffixes

### To Test
* GTF files, and their processing - the current modified GTFs may not be best options (especially for lncRNAs)
* Update genome, corresponding indexes to grc39 (for mouse and human)
* provide average read-length to prepDE (stringQuant function) computed from statsWriter

### AWS
m5a.2xlarge for **testing** on 2 paired end libraries
