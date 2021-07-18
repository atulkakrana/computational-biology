## Single-Cell RNA-Seq 

### Preprocessing
1. Download sample speciific FASTQ files to local machine or AWS instances

2. Use `count.sh` to generate gene-level counts and Loupe file for each sample
    1. change "--id" to any string/text that you want to name folder for output
    2. change "--sample" to exact sample prefix as FASTQ files

3. Use `aggregate.sh` to combine counts from all sample-speciific counts (and a new Loupe file)
    1. change "--id" to any string/text that you want to name folder for output
    2. change "--csv" to path for CSV file with sample_id, molecule_h5 file path ,batch , group

4. Explore expression and clusters using:
    1. Loupe browser (load `cloupe` file from aggregation step), or
    2. Use `scRNA.R` to explore the data

### Identify Cell-Types
Loupe file from the "aggregation" step can be used to check for key markers for each cluster and assign a cell-type based on prior knowledge. Third-party computational tools developed for automated identification of cell types.

Links:
1. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/tutorials/gex-analysis-nature-publication#header
2. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
