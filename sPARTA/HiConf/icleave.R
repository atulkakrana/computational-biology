##Scaling sRNA libraries - before degradome analysis

library(edgeR)

### DATA ###########################################################
getwd()
setwd('/data2/homes/kakrana/Rstudio/1.PAREPipe/data')
# files =c('AntW82_2.txt','AntW82s_1.txt','FlBW82_2.txt','FlBW82s_1.txt','Nod15_1.txt','Nod15_2.txt')
files.reduced = c('AntW82_2.txt','AntW82s_1.txt','Nod15_1.txt','Nod15_2.txt')
miRData <- readDGE(files.reduced) ## DGEList object
summary(miRData$counts)## Page 35 edgeR mannual- for why?
####################################################################

### Presence/Absence filter#########################################
### Filter low expression tags - Potential expression cutoff is 1-100 counts per million
presentList <- rowSums(cpm(miRData)>1) >= 2 ### How to calculate cpm cutoff? How many real tags it equals to? Use raw data not normalized
present.tags <- miRData[presentList,]
dim(present.tags)

### RESET LIB SIZES ################################################
##Page 36 edgeR
present.tags$samples$lib.size <- colSums(present.tags$counts)
present.tags$samples
####################################################################

### NORMALIZATION ##################################################
# data <- as.matrix(read.table("Norm_test.txt", sep="\t",row.names=1,as.is=TRUE))
# fact.RLE <-calcNormFactors(present.tags,method=c("RLE"))
# fact.RLE
# fact.TMM <- calcNormFactors(present.tags,method=c("TMM"))
# fact.TMM
# fact.upQT <- calcNormFactors(present.tags,method=c("upperquartile"))
# fact.upQT
norm.tags <- calcNormFactors(present.tags)

##FOr Exact test
norm.tags$samples$groups = c(1,1,2,2)## Group can be introduced 1) if data is read by sampleinfo file (page 34) 2) like this 3) in GLM - like microarray matrix
norm.tags$samples

##For GLM- Page 22
group <-factor(c(1,1,2,2))
design<- model.matrix(~0+group, data = norm.tags$samples)##  required to avoid intercept
colnames(design) <- c('FL','NO')
contrast.matrix <- makeContrasts(AvsB=FL-NO,BvsA=NO-FL,levels=design)
design
# fit <- lmFit(norm.tags,design)
# ##Effective library size
# effect.lib.size.upQT <- fact.upQT$samples$lib.size * fact.upQT$samples$norm.factors ## Seems the data is pre-normalized
# effect.lib.size.TMM <- fact.TMM$samples$lib.size * fact.TMM$samples$norm.factors
####################################################################

### DISPERSION ######################################################
plotMDS(norm.tags) ### Dimension 1 shows that difference is clear b/w conditions though replicated might vary
# fact.TMM <- estimateCommonDisp(fact.TMM, verbose=TRUE)
###Exact test
norm.tags <- estimateCommonDisp(norm.tags, verbose=TRUE)
norm.tags <- estimateTagwiseDisp(norm.tags)

###GLM
norm.tags <- estimateGLMCommonDisp(norm.tags, design, verbose=TRUE)
norm.tags <- estimateGLMTrendedDisp(norm.tags, design)
norm.tags <- estimateGLMTagwiseDisp(norm.tags, design)

fit <- glmFit(norm.tags, design)

plotBCV(norm.tags, cex=0.4) ## Plots tag wise dispersion
######################################################################

### DGE ##############################################################
### Using exact test instead of generalized linear model - because we are just interested in up/downregulated rather than stages in which up/down regulated
###EXACT######
exact.test <- exactTest(norm.tags)
topTags(exact.test,n=20)

##Write results
all.Res <- topTags(exact.test, n=Inf,sort.by='logFC')
setwd('/data2/homes/kakrana/Rstudio/1.PAREPipe')
write.table(all.Res, "ExactResCPM1.tsv", sep = '\t' )
allRownames <- rownames(topTags(exact.test, n=Inf,sort.by='logFC'))
selected <- all.Res[length(allRownames) < 23]

###Getting counts per million
detags <- rownames(topTags(exact.test, n=Inf))
countFile = cpm(norm.tags)[detags,]
write.table(countFile,"miRCount.txt", sep= '\t')

summary(de <- decideTestsDGE(exact.test))## Total genes at 5% FDR
results.default <- decideTestsDGE(exact.test, adjust.method = 'BH')##Default: method = separate

##Plot changed genes
detags <- rownames(norm.tags)[as.logical(de)]
plotSmear(exact.test, de.tags=detags)
abline(h=c(-1, 1), col="blue")

##GLM ######
lrt <- glmLRT(fit) ### Log liklihood ratio
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
summary(de <- decideTestsDGE(lrt))

######################################################################

###Normalize library
##We call the product of the original library size and the scaling factor the effective library size. The
##effective library size replaces the original library size in all downsteam analyses.
##http://grokbase.com/t/r/bioconductor/123gtpjqx4/bioc-edger-how-to-obtain-a-list-of-normalized-expression-values
# ##https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2011-February/001810.html 
# effective.lib.size <- fact.upQT$samples$lib.size * fact.upQT$samples$norm.factors
# effective.lib.size
# norm.Lib <- log2(t(t(fact.upQT$counts+0.25)/(effective.lib.size+0.25)))

##Also
###https://stat.ethz.ch/pipermail/bioconductor/2012-October/048532.html
y = cpm(data) ## Counts per million

# sampletype <- c('C1','C2')
# group <-factor(sampletype)
# design<- model.matrix(~0+group)##  required to avoid intercept
# y = estimateGLMCommonDisp(factors.TMM,design,verbose=TRUE)
# 
# y = estimateGLMTrendedDisp(y,design)
# y = estimateGLMTagwiseDisp(y,design)
# fit <- glmFit(y, design)


