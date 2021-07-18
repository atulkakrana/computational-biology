## scRNA-Seq 

### Preprocessing
1. Download sample speciific FASTQ files to local machine or AWS instances
2. Use `count.sh` to generate gene-level counts and Loupe file for each sample
3. Use `aggregate.sh` to combine counts from all sample-speciific counts (and a new Loupe file)
4. 


### Identify Cell Types
Loupe file from the "aggregation" step can be used to check for key markers for each cluster and assign a cell-type based on prior knowledge. Third-party computational tools developed for automated identification of cell types.