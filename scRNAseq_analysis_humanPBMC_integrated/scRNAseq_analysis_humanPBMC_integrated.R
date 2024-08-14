### Installation #####
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)


# Setup the Seurat Objects
library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset
InstallData("ifnb") #human PBMC (Peripheral Blood Mononuclear Cells)

# load dataset
ifnb = LoadData("ifnb")

# split the RNA measurements into two layers one for control cells, one for stimulated cells
ifnb[["RNA"]] = split(ifnb[["RNA"]], f = ifnb$stim)

ifnb #4 layers present: counts.CTRL, counts.STIM, data.CTRL, data.STIM

### Perform Analysis Without Integration ####

## run standard analysis workflow
# perform log-normalization on CTRL and STIM layers
ifnb = NormalizeData(ifnb)
# Identifies the most variable genes across cells
ifnb = FindVariableFeatures(ifnb)
# Centers and scales data matrix, preparation for downstream analysis
ifnb = ScaleData(ifnb)
# perform Principal Component Analysis on scaled data to reduce dimensionality and identify the main axes of variation.
ifnb = RunPCA(ifnb)

# Constructs a graph of nearest neighbors based on the first 30 principal components
ifnb = FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
# Identifies clusters of cells using the Louvain algorithm, with a resolution parameter that determines the granularity of clustering.
ifnb = FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

# Apples UMAP (Uniform Manifold Approximation and Projection) for visualization, using the first 30 principal components.
ifnb = RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# Creates UMAP plot of the data, grouping the cells by either their condition ('stim') or their cluster assignments ('seurat_clusters').
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

## Output UMAP Plot:
# Left Panel: grouping of cells by condition (CTRL vs STIM)
# Right Panel: grouping of cells by cluster (0-25)
ifnb

ifnb = IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] = JoinLayers(ifnb[["RNA"]])

ifnb = FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb = FindClusters(ifnb, resolution = 1)

ifnb = RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
DimPlot(ifnb, reduction = "umap", split.by = "stim")