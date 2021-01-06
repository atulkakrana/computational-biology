## RNA-Seq Ananlysis
## Use Rocket.v4 to process
## RNA-seq data to matrices

## Guide and Workflows
## https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual [Using StringTie with DESeq2]
## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

## NOTE:
## Irrespective of method to compute counts, these are reffered to as 'cpm' in this script
## The input file, if generated from Rocket/Stringtie includes TPM counts

## To Do
## [1] For input file from StringTie assign genes names to stringtie transcript (MSTRG.) using
###### 'IsoformSwitchAnalyzeR' as suggested: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
## [2] Annotate the Output from StringTie/IsoformSwitchAnalyzeR

## ENVIRONEMENT ###################################################
HOME = path.expand("~")
print(HOME)

## IMPORTS #########################################################
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

### DATA ###########################################################
## set workig directory
getwd()
setwd(paste(HOME,'0.work/1.seq-analysis/1.foxe3_rnaseq', sep="/"))

## read data and pheno file
countData   = read.table("gene_count_matrix.csv",header = T,row.names = 1, sep = ",")
colnames(countData); dim(countData)
phenoData   = read.table('phenofile.txt', sep ='\t', header = T);phenoData
colnames(phenoData)

## sanity check tp see all
## libs in phenodata are 
## present in counts data file
colmatch = all((phenoData$code) %in% colnames(countData))
ifelse(colmatch == TRUE, "All Good!", "STOP: Check colnames in pheno file and data file!!!")

## extract data for samples
## provided in the pheonofile
countData.sel = countData[, phenoData$code]
colnames(countData.sel)
dim(countData); dim(countData.sel)

## Prepare expression data
## object, with groups for 
## comparisions
grp   = phenoData$sample
myData = DGEList(counts=countData.sel, group=grp); myData
####################################################################

### Presence/Absence filter#########################################
## select genes that are robustly present
## the cutoff is arbitraly decided or 
## use the edgeR method below (filterByExpr)
# presentList     = rowSums(cpm(myData)>2) >= 2
# present.genes   = myData[presentList,]
# dim(myData);dim(present.genes)

presentList     = filterByExpr(myData, group=grp)
present.genes   = myData[presentList,, keep.lib.sizes=FALSE]
dim(myData);dim(present.genes)
present.genes[1:5,]

## convert to counts per million
## just for testing
# countDF.cpm   = cpm(present.genes$counts)
# write.table(countDF.cpm,'testCountsCPM2.tsv',sep='\t')
# countDF.cpm[1:5,]
####################################################################


### NORMALIZATION ##################################################
### Reset lib Sizes, and normalize
present.genes$samples$lib.size = colSums(present.genes$counts)
present.genes$samples
norm.genes = calcNormFactors(present.genes, method = 'TMM')
norm.genes; norm.genes$samples$norm.factors
boxplot(log2(getCounts(present.genes)+1),col = 'red')

##FOr Exact test
# norm.tags$samples$groups = c(1,1,2,2)## Group can be introduced 
# norm.tags$samples

# ##For GLM- Page 22
# group <-factor(c(1,1,2,2))
# design<- model.matrix(~0+group, data = norm.tags$samples)##  required to avoid intercept
# colnames(design) <- c('FL','NO')
# contrast.matrix <- makeContrasts(AvsB=FL-NO,BvsA=NO-FL,levels=design)
# design
# fit <- lmFit(norm.tags,design)
# ##Effective library size
# effect.lib.size.upQT <- fact.upQT$samples$lib.size * fact.upQT$samples$norm.factors ## Seems the data is pre-normalized
# effect.lib.size.TMM <- fact.TMM$samples$lib.size * fact.TMM$samples$norm.factors
####################################################################

### Plots ##########################################################
# To view the plot immediately
plotMDS.DGEList(norm.genes , main = "MDS Plot for Count Data", labels = phenoData$sequencer_code, cex=0.6)
boxplot(log2(getCounts(norm.genes)+1),col = 'green')
####################################################################

#### Visualize Batch effects ######################################
##Create Expression set
# exprs = norm.genes$counts
# class(exprs);dim(exprs);exprs[1:5,]
# phenoFile = read.table('phenoData.csv', sep =',', header = T,row.names = 1);phenoFile[1:5,]
# phenoFile.sel = phenoFile[-c(1,2,7),]
# rownames(phenoFile.sel)
# 
# # colnames(exprs) = rownames(phenoFile.sel);colnames(exprs)
# 
# head(exprs[,1:4])
# all(rownames(phenoFile.sel)==colnames(exprs))
# 
# require(Biobase)
# phenoData <- new("AnnotatedDataFrame",data=phenoFile.sel)
# phenoData
# exampleSet <- ExpressionSet(assayData=exprs,phenoData=phenoData)
# exampleSet
# 
# library(pvca)
# pct_threshold <- 0.6
# pData(exampleSet)
# batch.factors <- c("stage",)
# pvcaObj <- pvcaBatchAssess(exampleSet, batch.factors, pct_threshold)
# save.image("PVCA.RData")
# 
# bp <- barplot(pvcaObj$dat,  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"),las=2,main="Principle variance estimation bar chart")
# axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
# values = pvcaObj$dat
# new_values = round(values , 3)
# text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8) 
# 
# ##PCA
# plotPCA(exprs,groups = phenoFile$stage,addtext = phenoFile$stage,pcs = c(1,2),legend=F)
# plotPCA(exprs,groups = phenoFile$stage,addtext = phenoFile$stage,pcs = c(3,1),legend=F)
# plotPCA(exprs,groups = phenoFile$stage,addtext = phenoFile$stage,pcs = c(2,3),legend=F)

## PCA DATA ########################################################


##### Design #######################################################
# seq.group <-factor(phenoFile$stage)
# design = model.matrix(~0+seq.group);design
design = model.matrix(~0+grp);design
colnames(design) = gsub("grp", "", colnames(design))
colnames(design)
contrast.matrix = makeContrasts(foxe3-wt, levels = design)
design

# colnames(design) <- levels(norm.genes$samples$group)
####################################################################


#### Dispersion ####################################################
## Exact Test
# norm.genes.ex <- estimateCommonDisp(norm.genes, verbose=TRUE)
# norm.genes.ex <- estimateTagwiseDisp(norm.genes.ex)

## GLM
norm.genes.glm <- estimateGLMCommonDisp(norm.genes, design, verbose=TRUE) # Estimates common dispersions
norm.genes.glm <- estimateGLMTrendedDisp(norm.genes.glm, design) # Estimates trended dispersions
norm.genes.glm <- estimateGLMTagwiseDisp(norm.genes.glm, design) # Estimates tagwise dispersions
####################################################################


#### RESULTS ########################################################

### Exact Approach ###
# exact.test <- exactTest(norm.genes.ex)
# topTags(exact.test,n=20)
# all.Res <- topTags(exact.test, n=Inf)

### GLM approach
fit             = glmFit(norm.genes.glm, design) # Returns an object of class DGEGLM
logratio.test   = glmLRT(fit, contrast = contrast.matrix ) # Takes DGEGLM object and carries out the likelihood ratio test.
topTags(logratio.test,n=20)

all.Res         = topTags(logratio.test, n=Inf)
summary(decideTestsDGE(logratio.test))

### Write results
write.table(all.Res, "edgeR.tsv", sep = '\t' )

### Get counts
detags      = rownames(all.Res)
countFile   = cpm(norm.genes)[detags,]
write.table(countFile,"counts.tsv", sep= '\t')
##################################################################

############ Annotate and Combine ################################
## DG and Counts
final.res   = read.delim("edgeR.tsv")
final.counts= read.delim("counts.tsv")
anno.result = data.frame(final.res,final.counts)
write.table(anno.result,"rnaseq_results.tsv",sep='\t')

## Additional annotations
final.res[1:5,];dim(final.res)
res.IDs     = rownames(final.res)
annoFile    = read.delim("genes.attr_table",stringsAsFactors=FALSE) ## Generated by Rocket [not by v4]
anno        = annoFile[,c(1,4,5,6,7)]
anno.res    = anno[with (anno, match(res.IDs, anno$tracking_id)),]
anno.res[1:5,];dim(anno.res);dim(final.res)

anno.result = data.frame(anno.res,final.res,final.counts)
write.table(anno.result,"rnaseq_results.tsv",sep='\t')
#################################################################

##################### Plots #####################################
## Plot changed genes - Should be used when only results for single co-efficent are generated
summary(de <- decideTestsDGE(logratio.test))
DEGs <- as.logical(de)
plotSmear(logratio.test, de.tags=DEGs)
abline(h=c(-2, 2), col="red")
##################################################################

#################### ANNOTATIONS ################################
getwd()
setwd('/data2/homes/kakrana/3.ProjectR/3.Anther/2.Seq')
seq.results = read.delim("ResAnnotatedUnpairedTtest_NonRed.txt",head = TRUE)
dim(seq.results)

require(biomaRt)
listMarts() ## Check the marts - plant_mart_XX

# listDatasets(plantMart) ## Find species dataset
# plantMart = useMart("plants_mart_23",dataset="zmays_eg_gene")
attribs = listAttributes(plantMart);attribs[1:100,]
# filters = listFilters(plantMart);filters[1:100,]

plantMart = useMart("plants_mart_24")
plantMart = useDataset("zmays_eg_gene",mart = plantMart)

names(seq.results)
# test.ids = c('GRMZM2G180251','GRMZM2G019971')
geneIDs = seq.results$gene_short_name
length(geneIDs);geneIDs[1:5]

## Filters is like "where" of mySQL - In this case we specificy in query WHERE 'ensembl_gene_id' =test.ids
test.anno = getBM(attributes = c('ensembl_gene_id','chromosome_name','start_position','end_position','strand',"description"),
                  filters = 'ensembl_gene_id',values = geneIDs, mart = plantMart) ## Server is Down Sometimes
dim(test.anno);dim(seq.results);names(test.anno)
# write.table(test.anno,"seqLimmaAnnoTest.tsv",sep = '\t',row.names=F)

## Merge annotation to results frame
final.seq = merge(seq.results,test.anno,by.x="gene_short_name",by.y="ensembl_gene_id",all.x=TRUE)
dim(final.seq)
write.table(final.seq,"ResAnnotated2UnpairedTtest.txt",sep = "\t",row.names=FALSE)
