library(edgeR)

### DATA ###########################################################
getwd()
setwd('/data2/homes/kakrana/3.ProjectR/3.Anther/2.Seq')

countDF = read.table("genes.count_table",header = T,row.names = 1)
names(countDF); dim(countDF)
countDF.sel = countDF
dim(countDF);dim(countDF.sel)

phenoFile = read.table('AntherWalbotPheno.tsv', sep ='\t', header = T);phenoFile[1:5,]
colnames(countDF.sel) = phenoFile$Code
countDF.sel[1:5,];dim(countDF.sel)
names(phenoFile)
phenoFile
grp = phenoFile$Sample
myData = DGEList(counts=countDF.sel, group=grp) # Constructs DGEList object
####################################################################

### Presence/Absence filter#########################################

presentList <- rowSums(cpm(myData)>2) >= 2 ## How much CPM in how many samples
present.genes <- myData[presentList,]
dim(myData);dim(present.genes)

countDF.cpm = cpm(present.genes$counts)
write.table(countDF.cpm,'testCountsCPM2.tsv',sep='\t')
countDF.cpm[1:5,]
present.07 = countDF.cpm[countDF.cpm[[2]] > 10,];dim(present.07)

### Reset lib Sizes
present.genes$samples$lib.size <- colSums(present.genes$counts)
present.genes$samples
boxplot(log2(getCounts(present.genes)+1),col = 'red')
####################################################################


### NORMALIZATION ##################################################
norm.genes <- calcNormFactors(present.genes)

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
plotMDS.DGEList(norm.genes , main = "MDS Plot for Count Data",labels = phenoFile$Group2,cex=0.6)
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
# pvcaObj <- pvcaBatchAssess (exampleSet, batch.factors, pct_threshold)
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

##PCA DATA ####################################################


##### Design #######################################################
# seq.group <-factor(phenoFile$stage)
# design = model.matrix(~0+seq.group);design

design = model.matrix(~0+group,data=norm.genes$samples);design ## Tested OK
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)
contrast.matrix = makeContrasts((mac1_0_4mm+mac1_0_7mm+mac1_1_0mm+mac1_1_5mm+mac1_2_0mm+ms23_0_4mm+ms23_0_7mm+ms23_1_0mm+ms23_1_5mm+ms23_2_0mm+ocl4_0_4mm+ocl4_0_7mm+ocl4_1_5mm+ocl4_1mm+ocl4_2_0mm)/15 
              - (W23_0_2mm+W23_0_4mm+W23_0_7mm+W23_1_0mm+W23_1_5mm+W23_2_0mm+W23_2_5mm+W23_3mm+W23_4mm+W23_5mm)/10, levels = design)
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
save.image("AntherSeq.RData")

# load("AntherSeq.RData")
####################################################################


#### RESULTS ########################################################

### Exact Approach ###
# exact.test <- exactTest(norm.genes.ex)
# topTags(exact.test,n=20)
# all.Res <- topTags(exact.test, n=Inf)

### GLM approach
fit = glmFit(norm.genes.glm, design) # Returns an object of class DGEGLM

# contrast.matrix<-makeContrasts(groupmm0.4-groupmm1.5,levels=design)##Old WB and not new
# contrast.matrix<-makeContrasts(mm0.4-mm0.7,levels=design)##Old WB and not new

logratio.test = glmLRT(fit, contrast = contrast.matrix ) # Takes DGEGLM object and carries out the likelihood ratio test.
topTags(logratio.test,n=20)
all.Res <- topTags(logratio.test, n=Inf)
summary(dt <- decideTestsDGE(logratio.test))

### Write results
write.table(all.Res, "ResUnpairedTtest_CPM2.txt", sep = '\t' )
### Get counts
detags <- rownames(all.Res)
countFile = cpm(norm.genes)[detags,]
write.table(countFile,"TagUnpairedTtest.txt", sep= '\t')
##################################################################

############ Annotate and Combine ################################
final.res = read.delim("ResUnpairedTtest_CPM2.txt")
final.counts =  read.delim("TagUnpairedTtest.txt")

final.res[1:5,];dim(final.res)
res.IDs = rownames(final.res)

annoFile = read.delim("genes.attr_table",stringsAsFactors=FALSE) ## Generated by Rocket
anno = annoFile[,c(1,4,5,6,7)]
anno.res = anno[with (anno, match(res.IDs, anno$tracking_id)),]
anno.res[1:5,];dim(anno.res);dim(final.res)

anno.result = data.frame(anno.res,final.res,final.counts)
write.table(anno.result,"ResAnnotatedUnpairedTtest.txt",sep='\t')
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
