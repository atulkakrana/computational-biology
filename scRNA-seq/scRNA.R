## hands-on based on https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
## publication: Comprehensive integration of single cell data (2018)

library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data = Read10X(data.dir = "~/GoogleDrive/0.learn/3.omics/1.scrna/filtered_gene_bc_matrices/hg19")
pbmc      = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

## Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
## for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## filter cells that have unique feature counts over 2,500 or less than 200
## filter cells that have >5% mitochondrial counts
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#### NORMALIZE DATA ############################################
## global-scaling normalization method “LogNormalize” that normalizes the feature 
## expression measurements for each cell by the total expression, multiplies this 
## by a scale factor (10,000 by default), and log-transforms the result
pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## calculate a subset of features that exhibit high cell-to-cell variation in the dataset
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

## plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

## SCALING THE DATA ############################################
## step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate;
## to ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial
## contamination see https://satijalab.org/seurat/v3.0/sctransform_vignette.html
all.genes = rownames(pbmc)
pbmc      = ScaleData(pbmc, features = all.genes)


#### DIMENSIONALITY REDUCTION ###################################
## perform PCA on the scaled data, and
## examine and visualize PCA results
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")


## Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc = JackStraw(pbmc, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)

## provides a visualization tool for comparing the distribution of p-values 
## for each PC with a uniform distribution (dashed line)
JackStrawPlot(pbmc, dims = 1:15)

## alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based 
## on the percentage of variance explained by each one
ElbowPlot(pbmc)


#### CLUSTER CELLS #############################################

## embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, 
## with edges drawn between cells with similar feature expression patterns, and then 
## attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

## non-linear dimensional reduction techniques, such as tSNE and UMAP
## Examine and nd visualize tSNE results
pbmc = RunTSNE(object = pbmc)
DimPlot(pbmc, reduction = "tsne")

## Examine and nd visualize umap results
pbmc = RunUMAP(object = pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


#### FIND DIFFRENTIALLY EXPRESSED FEATURES ###################
## dentifes positive and negative markers of a single cluster (specified in ident.1),
## compared to all other cells. ## `min.pct`` argument requires a feature to be detected 
## at a minimum percentage in either of the two groups of cells. `thresh.test`` argument 
## requires a feature to be differentially expressed (on average) by some amount between the two groups.
## Find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# different tests for differential expression
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

## visualize expression for different genes
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
## visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
library(ggplot2)

