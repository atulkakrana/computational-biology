## Sample Graphics script
## Atul Kakrana

#### 24-nt PHAS loci v2 genome - Nov28 ####################################################################
########################################
getwd()
setwd("/new_data/data2/homes/kakrana/3.ProjectR/10.Asparagus")
exprs_data <- read.delim("24Phas.v2.txt",sep= "\t",header=TRUE,row.names=1,as.is=TRUE)
names(exprs_data);rownames(exprs_data);dim(exprs_data)
exprs_data[1:5,1:5]
data = as.matrix(exprs_data[13:18,1:65]); dim(data) ## Choose reproductive or all samples
test = data[,c(2,34,51,56)];test ## These will be removed
e = data[,-c(2,34,51,56)]; dim(e) ## Remove fasle positive

e=as.matrix(log2(data+1))

## Before clustering - Remove those with no variability i.e. all stages have 0 fold change
ind <- apply(e, 1, var) == 0
e <- e[!ind,];dim(e)

## Clustering - For 58 S/N panicle - CLUST1
hr <- hclust(as.dist(1-cor(t(e), method="pearson")), method="complete") ## Do you need to transpose??
hc <- hclust(as.dist(1-cor(e, method="spearman")), method="complete")


## Clustering - For 58 S/N Anther - CLUST2
# row_distance = dist(e, method = "euclidean") ## check available methods in 'dist' help file
# hr = hclust(row_distance, method = "complete") ## check available methods in 'hclust' help file
col_distance = dist(t(e), method = "manhattan")
hc = hclust(col_distance, method = "complete")

## Color
hmcol = colorRampPalette(c("white","white","pink","pink2","red","red3"))(1024)
# hmcol = colorRampPalette(brewer.pal(9,"Reds"))(512)
# hmcol = colorRampPalette(brewer.pal(6,"Blues"))(1024)

## Plot w/o clustering
# heatmap.2(e,col=hmcol,density.info="none",trace="none",scale='none',Colv=FALSE,Rowv=FALSE,
#           cexRow=0.8,cexCol=0.8,dendrogram =c('none'),key=TRUE,keysize=1)##lhei=c(2.5,5.0),lwid=c(2.5,5.0)

## Plot with clustering
heatmap.2(e,col=hmcol,density.info="none",trace="none",scale='column',Colv=as.dendrogram(hc),Rowv=FALSE,
          cexRow=0.6,cexCol=0.6,dendrogram =c('column'),key=TRUE,keysize=1)##lhei=c(2.5,5.0),lwid=c(2.5,5.0)

###########################################################################################################

################ Co-relation plot ####################################
######################################################################
require(psych)

getwd()

## Get data
setwd("/new_data/data2/homes/kakrana/3.ProjectR/8.Collab/5.UGA/co-relation/")
corr.data <- (read.table("lib-normalized-nt-normalized.abundances.veg.mod.w-inter.txt",sep= "\t",header=TRUE,row.names=1,as.is=TRUE))##
names(corr.data)
corr.data[1:5,1:5]

### Compute co-relations
res = corr.test(corr.data)
our.corr = round(res$r,6);our.corr
our.pvals = round(res$p,6);our.pvals

## Get colors
require(Heatplus);require(RColorBrewer);require(gplots)
hmcol = colorRampPalette(c("#40778a","grey88","white","white","grey88","#d33747"))(256)

## Make correlation plot
require(corrplot)
# corrplot(our.corr,p.mat = our.pvals,insig = "p-value", sig.level = -1,
#          method = "circle",type = "lower",col = hmcol,tl.cex=0.6,cl.cex=0.6) ## p.mat = cor.matrix$p, insig = "p-value"
corrplot(our.corr,method = "circle",type = "lower",col = hmcol,tl.cex=0.6,cl.cex=0.6) ## p.mat = cor.matrix$p, insig = "p-value"

## Make heatmaps
heatmap.2(our.corr,col=hmcol,density.info="none",trace="none",scale="none",Colv=F,
          Rowv=F,cexRow=0.6,cexCol=0.6,dendrogram =c('none'),
          key=TRUE,keysize=1)##lhei=c(2.5,5.0),lwid=c(2.5,5.0)

library(corrplot)
M <- cor(mtcars)
corrplot(M, method = "circle")
########################################################################

#### Plot-2 - direct IRs ###################################
##############################################################

require(ggplot2)
## Read ###
setwd("/new_data/data2/homes/kakrana/3.ProjectR/10.Asparagus/2.figs/IRplot")
adata <- read.delim2("direct.consensus.merged.size.mod.txt",header=T)
names(adata);dim(adata)
adata[1:15,]

adata$Abun.2 = as.numeric(as.character(adata$Abun.2)) ## Data transformed to numerical, required for floating number
adata$Abun = as.numeric(as.character(adata$Abun)) ## Data transformed to numerical, required for floating number
adata$Abun.sqrt = as.numeric(as.character(adata$Abun.sqrt)) ## Data transformed to numerical, required for floating number

adata[1:32,]
amax = max(adata$Abun.sqrt) 
bmin = min(adata$Abun.sqrt)
amax;bmin

### Prepare ###
maxpos = 554 ## Length from left to where plot is to be done
plot.data = adata[adata$Pos <= maxpos & (adata$Abun.sqrt < -1 | adata$Abun.sqrt > 1),]

abreaks = seq(from = 4, to = maxpos, by = 24)
bbreaks = seq(from = 1, to = maxpos, by = 24)
xbreaks = c(abreaks,bbreaks)


yabreaks = (sqrt(seq(from = 30000, to = amax*amax, by = 60000)));yabreaks
ybbreaks = (sqrt(seq(from = 30000, to = bmin*bmin, by = 60000)))*-1;ybbreaks
ybreaks=c(yabreaks,ybbreaks)

## Geometric point #########
p1 = ggplot(data=plot.data) ## alpha=Max.Coef
p1 + geom_point(data=plot.data,aes(x=Pos,y=Abun.sqrt,colour= as.factor(Size),alpha = abs(Abun.2)),size=2)+
  geom_hline(yintercept = c(-1,1), color = "grey50")+
  scale_x_continuous(limits=c(1,maxpos),breaks=xbreaks)+
  scale_y_continuous(limits=c(-450,650),breaks=ybreaks)+
  scale_colour_manual(values = c("21" = "grey50","22" = "grey50","23" = "grey50","24"="indianred1"))+
  
  geom_line(data=plot.data[plot.data$Size == 24 & plot.data$Abun.sqrt > 0,],
            aes(x=Pos,y=Abun.sqrt),size=0.1,fill="coral")+
  geom_line(data=plot.data[plot.data$Size == 24 & plot.data$Abun.sqrt < 0,],
            aes(x=Pos,y=Abun.sqrt),size=0.1,fill="coral")+
  
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=45,hjust = 1),panel.grid.minor.x = element_line(colour='grey100'),
        panel.grid.minor.y = element_line(colour='grey98'),panel.grid.major.x = element_line(colour='grey98'),
        panel.grid.major.y = element_line(colour='grey88'))+
  
  # Add cartesian specific breaks, this is added on second turn, start from P1+ggplot, and rerun from p1+geom_point
  for (i in abreaks){
    p1 = p1+geom_segment(x=i,y=1,xend=i,yend=Inf,linetype = 3,color = "grey88",size=0.1)
    p1 = p1+geom_segment(x=i-3,y=-1,xend=i-3,yend=-Inf,linetype = 3,color = "grey88",size=0.1)
  }
######################################################################################

## PCA Plots- All Affy ##############################################
library("FactoMineR")
library("factoextra")
require("scatterplot3d")
require("car")

setwd("/Users/atul/Google Drive/3.Project-R/1.iSyTE/co-relation")
af.data.all = (read.table("AF_all_exp.txt",sep= "\t",header=TRUE,row.names=1,as.is=TRUE)); af.data.all[1:5,1:5]

## Method 1 , generates a lot of info
af.pca      = PCA(af.data.all[,1:34],  graph = FALSE);af.pca ## Method 1 , generates a lot of info
af.pca.mat  = af.pca$var$cor
fviz_screeplot(af.pca, addlabels = TRUE, ylim = c(0, 50)) ## Shows which PCAs are more informative (usually first 3)

## Method 2
af.prcomp   = prcomp(af.data.all[,1:34],  scale = TRUE);af.prcomp 
af.pca.mat  = af.prcomp$rotation

af.pca.mat <- cbind(af.pca.mat, c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                  2,2,2,2,2,2,2,2,2,2,2,2,3,3,3));af.pca.mat ## Stage-specific labels
# af.pca.mat <- cbind(af.pca.mat, c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,
#                                   3,3,3,3,3,3,3,4,4,4,4,4,5,5,5));af.pca.mat ## refines stage-specific labels
# colors <- brewer.pal(n=8, name="Dark2")
colors = colorRampPalette(c("#2ecc71","#e74c3c","#34495e"))(3)

# af.pca.mat <- cbind(af.pca.mat, c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,5,5,6,6,
#                                   7,7,7,7,8,8,8,9,9,10,10,10,11,11,11));af.pca.mat ## Time-specific labels
# colors <- brewer.pal(n=11, name="Set3")

colnames(af.pca.mat)[35] = c("grp");af.pca.mat

## Plot first three PCA
library("RColorBrewer")

scatter3d(x = af.pca.mat[,1], y = af.pca.mat[,2], z = af.pca.mat[,3], grid = FALSE, surface=F, ellipsoid = F,
          group = as.factor(af.pca.mat[,35]),labels = rownames(af.pca.mat),id.n=nrow(af.pca.mat), level = 0.6,
          xlab = "PCA1", ylab = "PCA2",zlab = "PCA3", axis.scales = FALSE, surface.col = colors,
          ellipsoid.alpha = 0.1, parallel =F)

## w/o labels
scatter3d(x = af.pca.mat[,1], y = af.pca.mat[,2], z = af.pca.mat[,3], grid = FALSE, surface=F, ellipsoid = F,
          group = as.factor(af.pca.mat[,35]), level = 0.6,
          xlab = "PCA1", ylab = "PCA2",zlab = "PCA3", axis.scales = FALSE, surface.col = colors,
          ellipsoid.alpha = 0.1, parallel =F, sphere.size = 1.8)

# Graph of variables: default plot
fviz_pca_var(af.pca,geom = c("point", "text"), col.var = "black")



###############################################
###### Aspa Bean plot - relative enrichment
require(ggplot2)
setwd("/new_data/data2/homes/kakrana/3.ProjectR/10.Asparagus/2.figs")
adata <- read.delim("24Phas.v2.forBoxplot.txt",sep= "\t",header=TRUE,row.names=1,as.is=TRUE)
names(adata);dim(adata);adata[1:5,1:5]

## Extract stages, add factor, and rbind
pre   = adata[,8];pre[1:5] ## Aspa
agrp  = rep(1,length(pre))
pre.data = cbind(pre,agrp)
colnames(pre.data) = c("ratio","grp")
pre.data[1:5,]

mei   = adata[,9];mei[1:5] ## Aspa
bgrp  = rep(2,length(mei))
mei.data = cbind(mei,bgrp)
colnames(mei.data) = c("ratio","grp")
mei.data[1:5,]

post   = adata[,10];post[1:5] ## Aspa
cgrp  = rep(3,length(post))
post.data = cbind(post,cgrp)
colnames(post.data) = c("ratio","grp")
post.data[1:5,]

post2   = adata[,11];post2[1:5] ## Aspa
dgrp  = rep(4,length(post2))
post2.data = cbind(post2,dgrp)
colnames(post2.data) = c("ratio","grp")
post2.data[1:5,]

post3   = adata[,12];post2[1:5] ## Aspa
egrp  = rep(5,length(post3))
post3.data = cbind(post3,egrp)
colnames(post3.data) = c("ratio","grp")
post3.data[1:5,]

plot.data1 = data.frame(rbind(pre.data,mei.data,post.data,post2.data,post3.data))
plot.data1[1:5,];dim(plot.data1)

## Filters
plot.data = plot.data1[plot.data1$ratio > 3,];dim(plot.data)

## Plot
p2 = ggplot(data=plot.data,aes(x = factor(grp), y=ratio))
p2+geom_violin(width=0.7,colour="firebrick",fill="firebrick",Trim = FALSE)+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour='white'),
        panel.grid.major.y = element_line(colour='grey48'))

p2+geom_boxplot(notch = FALSE,width=0.5,colour="firebrick",fill="firebrick")+
  geom_point(shape=1)+
  theme_bw()+
  theme(panel.grid.minor.y = element_line(colour='white'),
        panel.grid.major.y = element_line(colour='grey48'))

