## Process and explore
## single-cell RNA-seq 
## data; the output from
## 10x cellranger aggregate
## function is the input
## for this code

## IMPORTS
library(dplyr)
library(Seurat)
library(ggplot2)

## USER SETTINGS
EXP_NAME  = "E16_P0" ## Give some name to your expereimnt, used to create folders
DATA_FOLD = "/home/anand/0.work/1.scseq/run_aggr_E16_5_P0/outs/count/filtered_feature_bc_matrix"


## ANALYSIS

# # Load the PBMC Example Dataset
# pbmc.data = Read10X(data.dir = "~/GoogleDrive/0.learn/3.omics/1.scrna/filtered_gene_bc_matrices/hg19")
# pbmc      = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc

expr.data = Read10X(data.dir = DATA_FOLD)
expr = CreateSeuratObject(counts = expr.data, project = EXP_NAME, min.cells = 3, min.features = 200)
# expr[,1:10]

## Lets examine a few genes in the first thirty cells
expr.data[c("Mid2", "Rbm41"), 1:30] ## change to any genes that exist in data

## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern = "^mt-") ## Try "^MT-" if all values as same error

## Show QC metrics for the first 5 cells
head(expr@meta.data, 5)

## Visualize QC metrics as a violin plot
VlnPlot(expr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
## for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 = FeatureScatter(expr, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(expr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("feats.png" , plot = last_plot(), dpi = 300, bg = "white") ## set bg = NULL for transparent background
CombinePlots(plots = list(plot1, plot2))
# ggsave("feats.png" , plot = last_plot(), dpi = 300, bg = "white") ## set bg = NULL for transparent background

## filter cells that have unique feature counts over 2,500 or less than 200
## filter cells that have >5% mitochondrial counts
expr = subset(expr, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#### NORMALIZE DATA ############################################
## global-scaling normalization method “LogNormalize” that normalizes the feature 
## expression measurements for each cell by the total expression, multiplies this 
## by a scale factor (10,000 by default), and log-transforms the result
expr  = NormalizeData(expr, normalization.method = "LogNormalize", scale.factor = 10000)

## calculate a subset of features that exhibit high cell-to-cell variation in the dataset
expr = FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2500)

## Identify the 10 most highly variable genes
top20 = head(VariableFeatures(expr), 20)
top20

## plot variable features with and without labels
plot1 <- VariableFeaturePlot(expr)
ggsave("vfp.png" , plot = last_plot(), dpi = 300, bg = "white") ## set bg = NULL for transparent background

plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
ggsave("vfp_labelled.png" , plot = last_plot(), dpi = 300, bg = "white") ## set bg = NULL for transparent background

# CombinePlots(plots = list(plot1, plot2))
# ggsave("vfp_final.png" , plot = last_plot(), dpi = 300, 
#             width=14, height=8, bg = "white") ## set bg = NULL for transparent background

## SCALING THE DATA ############################################
## step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate;
## to ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial
## contamination see https://satijalab.org/seurat/v3.0/sctransform_vignette.html
all.genes = rownames(expr)
expr      = ScaleData(expr, features = all.genes)


#### DIMENSIONALITY REDUCTION ###################################
## perform PCA on the scaled data, and
## examine and visualize PCA results
expr = RunPCA(expr, features = VariableFeatures(object = expr))
print(expr[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(expr, dims = 1:2, reduction = "pca")

DimPlot(expr, reduction = "pca")
ggsave("pca.png" , plot = last_plot(), dpi = 300, bg = "white")


## Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
expr = JackStraw(expr, num.replicate = 100)
expr = ScoreJackStraw(expr, dims = 1:20)

## provides a visualization tool for comparing the distribution of p-values 
## for each PC with a uniform distribution (dashed line)
JackStrawPlot(expr, dims = 1:15)
ggsave("jackStrawPlot.png" , plot = last_plot(), dpi = 300, bg = "white")

## alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based 
## on the percentage of variance explained by each one
ElbowPlot(expr)
ggsave("ElbowPlot.png" , plot = last_plot(), dpi = 300, bg = "white")


#### CLUSTER CELLS #############################################
## embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, 
## with edges drawn between cells with similar feature expression patterns, and then 
## attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’
expr = FindNeighbors(expr, dims = 1:10)
expr = FindClusters(expr, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(expr), 5)

## non-linear dimensional reduction techniques, such as tSNE and UMAP
## Examine and nd visualize tSNE results
expr = RunTSNE(object = expr)
DimPlot(expr, reduction = "tsne")
ggsave("tsne.png" , plot = last_plot(), dpi = 600, bg = "white")

## Examine and nd visualize umap results
expr = RunUMAP(object = expr, dims = 1:10)
DimPlot(expr, reduction = "umap")
ggsave("umap.png" , plot = last_plot(), dpi = 600, bg = "white")


#### FIND DIFFRENTIALLY EXPRESSED FEATURES ###################
## identifes positive and negative markers of a single cluster (specified in ident.1),
## compared to all other cells. ## `min.pct`` argument requires a feature to be detected 
## at a minimum percentage in either of the two groups of cells. `thresh.test`` argument 
## requires a feature to be differentially expressed (on average) by some amount between the two groups.

## Find all markers of provided cluster
clust_x = 6
clusterX.markers <- FindMarkers(expr, ident.1 = clust_x, min.pct = 0.25)
head(clusterX.markers, n = 10)

## Find all markers distinguishing cluster in ident.1 from clusters in ident.2
clust_x = 6
clust_y = c(0, 3)
clusterX.markers <- FindMarkers(expr, ident.1 = clust_x, ident.2 = clust_y, min.pct = 0.25)
head(clusterX.markers, n = 10)

## Find markers for every cluster compared to all remaining cells, report only the positive ones
avg_logFC_cutoff = 1
expr.markers <- FindAllMarkers(expr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
expr.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC_cutoff)
write.csv(expr.markers, file = "top_markers.tsv",row.names=TRUE, sep = "\t")

## Different tests for differential expression
clust_x = 6
clusterX.markers <- FindMarkers(expr, ident.1 = clust_x, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)



## VISUALIZATION
## visualize expression for different genes
VlnPlot(expr, features = c("Mafg", "Hmox1"))
ggsave("VlnPlot.png" , plot = last_plot(),   dpi = 300, bg = "white", width=12, height=8)

## Plot raw counts as well
VlnPlot(expr, features = c("Mafg", "Hmox1"), slot = "counts", log = TRUE)
ggsave("VlnPlotRawCounts.png" , plot = last_plot(), dpi = 300, bg = "white", width=12, height=8)

## Visualizes feature expression on a tSNE or PCA plot
FeaturePlot(expr, features = c("Mafg", "Hmox1"))
ggsave("FeatsExpression.png" , plot = last_plot(), dpi = 300, bg = "white", width=12, height=8)




#### NOTES ###
##############
## Based on https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
## publication: Comprehensive integration of single cell data (2018)

