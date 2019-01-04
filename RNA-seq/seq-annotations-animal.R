## Script to add annotations to RNA-seq data
## https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

### Annotationa packages
require(AnnotationDbi)
require(Mus.musculus)
require(org.Mm.eg.db)
# require(EnsDb.Mmusculus.v79)

### Read output file
setwd('/home/kakrana/8.Collab/4.Kak/5.lncrna/8.wbdata/WB_RNAseq2'); getwd()
seq.results = read.delim("gene_count_matrix.csv",head = TRUE,sep = ',')
dim(seq.results)

### Check keys and output columns
### in the selected annotation database
### for key type to use and for 
### selecting annotations
annodb = Mus.musculus
# annodb = org.Mm.eg.db
message("Keytypes"); keytypes(annodb);
message("Columns"); columns(annodb)

### Get annotations using organism level package
### Genes in haplotypic region will have more
### then one annotation
keys    = seq.results[,1]; head(keys)
cols    = c("SYMBOL","GENENAME","ENTREZID","CHR","TXID")
annots  = select(annodb, keys=as.character(keys),columns=cols, keytype="ENSEMBL")
message("Annotations table dim:");dim(annots)

### Merge annotation to results frame and
### seprate annotated one for collapsing
merged.annots = merge(seq.results,annots,by.x="gene_id",by.y="ENSEMBL",all.x=TRUE); head(merged.annots)
merged.annots.na = merged.annots[is.na(merged.annots$ENTREZID),]
merged.annots.1  = merged.annots[!is.na(merged.annots$ENTREZID),]

message("Seq results dim:");dim(seq.results)
message("Merged dataframe dim:");dim(merged.annots)
message("Rows with no annotation:");dim(merged.annots.na)
message("Rows with an annotation:");dim(merged.annots.1)


### Collapse on EntrezID/GENEID (both seems same)
### this will also remove N/A (seprate them, collapse and join)
merged.annots.2 = merged.annots.1[!duplicated(merged.annots.1$ENTREZID),]
message("collapsed annotated rows dim:");dim(merged.annots.2)

### Merge collapsed annotated genes with unannotated ones
### Genes with different ENTREZID but same ENSEMBL may be
### present as these belong to haplotypic region
merged.annots.3 = rbind(merged.annots.2,merged.annots.na)
message("merged and collapsed annotated rows dim:");dim(merged.annots.3)

## Write output table
write.table(merged.annots.3,"anno-table.txt", sep = '\t',row.names=FALSE)

