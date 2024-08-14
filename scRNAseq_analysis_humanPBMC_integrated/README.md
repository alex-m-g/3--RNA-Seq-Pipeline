### scRNA Sequence Analysis of hPBMC using Integration Method

We will integrate data from the two conditions (CTRL & STIM), so that cells from the same cell type/subpopulation will cluster together.

Procedure is referred to as integration/alignment. When aligning two genome sequences together, identification of shared/homologous regions can help to interpret differences between the sequences as well.

# Goal
To learn shared cell types/states to help enable us to compare control stimulated and control profiles for these individual cell types.

## Seurat v5 Integration
Returns a single dimensional reduction that captures the shared sources of variance across multiple layers, so that cells in a similar biological state will cluster. The method that returns a dimensionsla reduction which can be used for visualization and unsupervised clustering analysis.
For evaluating performance, we use cell type labels.
