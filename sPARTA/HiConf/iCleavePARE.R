### PARE analysis - No replicates

library(edgeR)
library(splines)

### DATA ###########################################################
getwd()
setwd('/data2/homes/kakrana/Rstudio/1.PAREPipe/data')##Office
# files =c('15DNodPARE.txt','UOF_dPARE.txt','Ant_dPARE.txt')
files.reduced = c('Ant_dPARE.txt','15DNodPARE.txt')
data <- readDGE(files.reduced) ## DGEList object
summary(data$counts)## Page 35 edgeR mannual- for why?
####################################################################

### Presence/Absence filter#########################################
### Filter low expression tags - Potential expression cutoff is 1-100 counts per million
presentList <- rowSums(cpm(data)>0.5) >= 1 ### How to calculate cpm cutoff? How many real tags it equals to? Use raw data not normalized
present.tags <- data[presentList,]
dim(present.tags)
####################################################################

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
####################################################################

### DE #############################################################
## Calculating BCV if no replicates
norm.tags2 <- norm.tags
norm.tags2$samples$groups = c(1,1) ## Grouping just to calculate BCV
norm.tags2 <- estimateDisp(norm.tags2,robust=TRUE,winsor.tail=c(0.05,0.2)) ## From https://stat.ethz.ch/pipermail/bioconductor/2013-August/054239.html
plotBCV(norm.tags2)
norm.tags$samples$groups = c(1,2)## Re-setting the grouping to correct order
results = exactTest(norm.tags,dispersion=norm.tags2$trended.dispersion)### Not working because group is just one - E-mail

### Arbitrarily choosen bcv - EdgeR page 18
norm.tags$samples$groups = c(1,2)
norm.tags$samples
bcv <- 0.3
results = exactTest(norm.tags,dispersion=bcv^2)

### Write DE results
all.Res <- topTags(results, n=Inf,sort.by='logFC')
setwd("/data2/homes/kakrana/Rstudio/1.PAREPipe/")
write.table(all.Res, "ExactResPARE_CPM05.tsv", sep = '\t' )
summary <- summary(de <- decideTestsDGE(results))## Total genes at 5% FDR
write(summary, "ExactResPARE_CPM05Sum.txt")

##################################################################

#### GET COUNTS ##################################################
###Getting counts per million
detags <- rownames(topTags(results, n=Inf))
countFile <- cpm(norm.tags)[detags,]
write.table(countFile,"PAREcountCPM05.tsv",sep = '\t' )
##################################################################