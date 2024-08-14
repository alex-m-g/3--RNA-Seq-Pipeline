# scRNA-seq Integration
Integration of single-cell sequencing datasets across experimental batches, donors, or conditions.
Integrative analysis helps match shared cell types and cell states across datasets. This is important to boost statistical power and facilitate accurate comparative analysis across datasets.

## Integration Goals
- Identify cell subpopulations that are present in both datasets.
- Obtain cell-type markers that are conserved in both controls and stimulated cells.
- Compare the datasets to find cell-type specific responses to stimulation.

## Dataset
- The dataset comes from human PBMC (Peripheral Blood Mononuclear Cells) from two conditions: (1) interferon-stimulated and (2) control cells.
- We will integrate the two conditions together to jointly identify cell subpopulations across datasets, and then explore how each group differs across conditions.

## Plot Interpretation

### UMAP Axes (umapunintegrated_1 and umapunintegrated_2):
The axes of the plot represent the two-dimensional space created by UMAP, a dimensionality reduction technique. UMAP compresses high-dimensional data (in this case, the gene expression profiles) into two dimensions while preserving the relationships between data points as much as possible.
umapunintegrated_1 and umapunintegrated_2 are the first and second dimensions of the UMAP projection, respectively. Cells that are close to each other in this space have similar gene expression profiles.

### Clusters:
Each point on the plot represents a single cell.
The cells are grouped into clusters, which are color-coded. These clusters were identified using the Louvain algorithm, which finds groups of cells with similar gene expression patterns.
The different colors represent different clusters, with each color corresponding to a cluster labeled with a number (e.g., 0, 1, 2, etc.).
Cluster Labels: The cluster numbers (e.g., 0, 1, 2, etc.) shown in the legend on the right correspond to the colors on the UMAP plot. Each cluster may represent a distinct cell type, cell state, or functional group within the dataset.

### Interpretation of Clusters:
Cells within the same cluster share similar gene expression profiles, suggesting they may have similar biological functions, be of the same cell type, or be in the same state.
Clusters that are spatially close to each other in the UMAP plot may be more similar to each other in terms of gene expression than those that are far apart.

### Biological Insights:
Researchers use these clusters to identify and annotate different cell types or states. For instance, one cluster might represent T cells, another might represent B cells, and so on, depending on the biological context of the experiment.
If specific genes or markers are known to be associated with certain cell types or states, they can be used to further annotate these clusters.

### Example Interpretation:
If, for example, cluster 0 (colored pink) is found to contain cells that express high levels of a gene marker specific to T cells, then this cluster might be annotated as representing T cells.
The proximity or separation of clusters can also provide insight into how different cell populations relate to one another. For example, if clusters 1 and 2 are close together, it might suggest that these cell types are closely related or have similar functions.

### Summary:
The right panel UMAP plot provides a visual representation of how the cells in the dataset group into distinct clusters based on their gene expression profiles. Each cluster likely represents a unique cell population or state, and the spatial relationships between these clusters offer insights into the similarities and differences between these populations.
