import os
import tarfile
import ftplib
import gzip
import shutil
import gc

# Step 1: Download the TAR file via FTP
def download_ftp_file(ftp_url, ftp_directory, filename, local_filename):
    with ftplib.FTP(ftp_url) as ftp:
        ftp.login()  # Anonymous login
        ftp.cwd(ftp_directory)  # Change to the directory where the file is located
        with open(local_filename, 'wb') as f:
            ftp.retrbinary(f"RETR {filename}", f.write)
    print(f"Downloaded {local_filename}")

# Step 2: Extract the TAR file (without GZ mode)
def extract_tar_file(tar_path, extract_folder):
    with tarfile.open(tar_path, "r") as tar:  # No "r:gz", just "r"
        tar.extractall(path=extract_folder)
    print(f"Extracted to {extract_folder}")

# Step 3: Unzip the .csv.gz files
def decompress_gz_files(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.gz'):
                gz_path = os.path.join(root, file)
                csv_path = gz_path[:-3]  # remove .gz extension
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(csv_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                print(f"Decompressed {file} to {csv_path}")

# Your integration function
def integration(folder, cutoff_value):
    out = []
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        
        if file_path.endswith('.csv'):  # Only process CSV files
            processed_data = pp(file_path, cutoff_value)  # Pass the cutoff_value
            out.append(processed_data)
            
            # Clear memory
            del processed_data
            gc.collect()

    # Concatenate all the processed AnnData objects
    adata = sc.concat(out)
    
    # Filter genes and save
    sc.pp.filter_genes(adata, min_cells=10)
    adata.X = csr_matrix(adata.X)  # Compress data for memory efficiency
    adata.write_h5ad('combined.h5ad')
    
    return f"File: 'combined.h5ad' has been written."

def detect_cutoff(df, prediction_col='prediction', diff_col='diff', bins=20, drop_threshold=-30):
    """
    Detects the cutoff point in the 'diff' column where the histogram starts to drop significantly.

    Parameters:
    - df: Pandas DataFrame containing the data.
    - prediction_col: Column name for the doublet/singlet prediction.
    - diff_col: Column name for the difference between doublet and singlet probabilities.
    - bins: Number of bins to use for the histogram.
    - drop_threshold: Percentage drop between bins to detect a significant drop.
    
    Returns:
    - cutoff: The bin edge where a significant drop in counts is detected.
    """

    # Filter the dataframe to only include rows predicted as 'doublet'
    doublet_data = df[df[prediction_col] == 'doublet']

    # Plot the histogram for visualization
    sns.displot(doublet_data, x=diff_col, bins=bins)
    plt.show()

    # Compute the histogram manually to analyze counts and bin edges
    counts, bin_edges = np.histogram(doublet_data[diff_col], bins=bins)

    # Calculate the percentage change between consecutive bins
    diff_counts = np.diff(counts)
    percentage_change = diff_counts / counts[:-1] * 100

    # Detect the first significant drop in percentage change
    drop_index = np.where(percentage_change < drop_threshold)[0]

    if len(drop_index) > 0:
        # Return the corresponding bin edge as cutoff
        cutoff = bin_edges[drop_index[0] + 1]
        return cutoff
    else:
        # If no significant drop is found, return None or a default cutoff
        return None

def detect_pca_elbow(adata, threshold=-1):
    """
    Detects the elbow point in the PCA variance ratio plot.
    
    Parameters:
    - adata: AnnData object after PCA.
    - threshold: A threshold for significant drop in variance ratio.
    
    Returns:
    - Elbow point (integer) to use for nearest neighbor calculation.
    """
    variance_ratios = adata.uns['pca']['variance_ratio']
    diff_ratios = np.diff(variance_ratios)

    # Find the index where the drop in variance ratio is less than the threshold
    elbow_index = np.where(diff_ratios > threshold)[0]

    if len(elbow_index) > 0:
        return elbow_index[0] + 1  # Adding 1 because diff shifts index by 1
    else:
        return len(variance_ratios)  # Use all PCs if no elbow point is found

def parameter_adjust(folder):
    # Loop through the files in the folder
    for file in os.listdir(folder):
        # Construct the full file path
        file_path = os.path.join(folder, file)
        
        # Read the first file into an AnnData object
        adata = sc.read_csv(file_path)
        
        # Transpose the data to have cells as rows and genes as columns
        adata = adata.T
        
        #keep the genes that are found in at least 10 cells
        sc.pp.filter_genes(adata, min_cells=10)

        #keep the top 2000 variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')

        # model setup
        scvi.model.SCVI.setup_anndata(adata)
        # train vae model
        vae = scvi.model.SCVI(adata)
        vae.train()

        #train the solo model, which predicts doublets
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train()
        # pass soft=False to make a new column that is the predicted label (doublet or singlet)
        df = solo.predict()
        df['prediction'] = solo.predict(soft=False)

        # make a new column that is the difference in probability of doublet and singlet columns
        df['diff'] = df.doublet - df.singlet

        sns.displot(df[df.prediction == 'doublet'], x='diff')

        global cutoff_value
        cutoff_value = detect_cutoff(df)

        doublets = df[(df['prediction'] == 'doublet') & (df['diff'] > cutoff_value)]
        adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
        adata = adata[~adata.obs.doublet]

        # mitochondrial genes are annotated as "MT-" in the gene names, lets filter that.
        adata.var['mt'] = adata.var.index.str.startswith('MT-')

        ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
        ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)

        # same as we did with mitochondrial genes, we label it to our main dataframe
        adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

        # calculate qc metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
        
        adata.var.sort_values('n_cells_by_counts')

        sc.pp.filter_genes(adata, min_cells=3)

        adata.var.sort_values('n_cells_by_counts')
        # now every gene listed is in at least 3 cells

        adata.obs.sort_values('total_counts')
        adata = adata[adata.obs.pct_counts_mt < 20] # 20 percent cutoff, will not get rid of any mitochondrial cells but good practice
        adata = adata[adata.obs.pct_counts_ribo < 2]

        sc.pp.normalize_total(adata, target_sum=1e4) # normalize every cell to 10,000 UMI
        sc.pp.log1p(adata) #change to log counts

        adata.raw = adata
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pl.highly_variable_genes(adata)
        # Filter out the non-highly variable genes
        adata = adata[:, adata.var.highly_variable]

        # we're going to regress out the differences that arise due to the total number of counts mitochondrial counts and the ribosomal counts
        #  this will filter out some data that are due to processing and just sample quality, sequencing artifact, etc.
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])

        #normalize each gene to the unit variance of that gene
        sc.pp.scale(adata, max_value=10)

        #Run component analysis to further reduce the dimensions of the data
        sc.tl.pca(adata, svd_solver='arpack') # this calculates 50 pcs

        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

        global elbow_point
        elbow_point = detect_pca_elbow(adata)

        sc.pp.neighbors(adata, n_pcs=elbow_point)

        sc.tl.umap(adata)

        sc.tl.leiden(adata, resolution = 1)

        sc.pl.umap(adata, color=['leiden'])

        # Break after the first file is processed
        break

def processing(combined_file):
    adata = sc.read_h5ad(combined_file)
    adata.obs.groupby('Sample').count()
    sc.pp.filter_genes(adata, min_cells=100)
    # save the raw data into 'counts' layer
    adata.layers['counts'] = adata.X.copy()
    # normalize the counts to 10,000 to every cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = (adata)

    scvi.model.SCVI.setup_anndata(adata, layer='counts',
                              categorical_covariate_keys=["Sample"], #add "batch" if pre-processed in batches; "technology" if pre-processed from multiple computers.
                              continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])
    
    model = scvi.model.SCVI(adata)
    model.train()

    adata.obsm['X_scVI'] = model.get_latent_representation()
    adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)
    sc.pp.neighbors(adata, use_rep = 'X_scVI')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1)
    sc.pl.umap(adata, color=['leiden','Sample'], frameon=False)
    sc.tl.rank_genes_groups(adata, 'leiden')
    markers = sc.get.rank_genes_groups_df(adata, None)
    markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
    markers_scvi = model.differential_expression(groupby = 'leiden')
    markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > 0.5)]
    sc.pl.umap(adata, color = ['leiden'], frameon=False, legend_loc = 'on data')

    ### START ANNOTATION ###
    cell_type_col, annotation_dict = cell_annotation_automatic('combined.h5ad')
    # Add cell type column based on annotation
    adata.obs['cell_type'] = adata.obs.leiden.map(annotation_dict)

    sc.pl.umap(adata, color = ['cell_type'], frameon=False)
    # Save UMAP plots
    umap_path_1 = os.path.join(output_folder, 'umap_cell_type.png')
    sc.pl.umap(adata, color=['cell_type'], frameon=False, save=umap_path_1)

    sc.pl.umap(adata, color = ['cell_type'], frameon=False, legend_fontsize=6, legend_loc='on data')
    umap_path_2 = os.path.join(output_folder, 'umap_cell_type_legend.png')
    sc.pl.umap(adata, color=['cell_type'], frameon=False, legend_fontsize=6, legend_loc='on data', save=umap_path_2)

    adata.uns['scvi_markers'] = markers_scvi
    adata.uns['markers'] = markers

    adata.write_h5ad('integrated.h5ad')
    model.save('model.model')

    return f"Files: 'integrated.h5ad' and 'model.model' have been written."

def analysis(main_file, model_file):
    adata = sc.read_h5ad(main_file)
    def map_condition(x):
        if 'cov' in x:
            return 'COVID19'
        else:
            return 'control'
    adata.obs['condition'] = adata.obs.Sample.map(map_condition)
    num_tot_cells = adata.obs.groupby(['Sample']).count()
    num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
    cell_type_counts = adata.obs.groupby(['Sample', 'condition','cell_type']).count()
    cell_type_counts = cell_type_counts[cell_type_counts.sum(axis=1) > 0].reset_index()
    cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
    cell_type_counts['total_cells'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)
    cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells

    print(f"Cell expression frequency between sample groups")
    # Saving the boxplot
    plt.figure(figsize=(10,4))
    ax = sns.boxplot(data=cell_type_counts, x='cell_type', y='frequency', hue='condition')
    plt.xticks(rotation=35, rotation_mode='anchor', ha='right')
    boxplot_path = os.path.join(output_folder, 'cell_expression_frequency_boxplot.png')
    plt.savefig(boxplot_path)
    plt.close()  # Close the plot to avoid memory issues

    subset = adata[adata.obs['cell_type'].isin(['AT1', 'AT2'])].copy()
    subset.X = subset.X.toarray()
    sc.pp.filter_genes(subset, min_cells=100)
    res = de.test.wald(
        data=subset,
        formula_loc="~ 1 + condition",
        factor_loc_totest="condition"
    )
    dedf = res.summary().sort_values('log2fc', ascending=False).reset_index(drop=True)
    #finding the most upregulated gene expression between the control and covid patient groups
    most_up = dedf.iloc[0].gene
    i = np.where(subset.var_names == most_up)[0][0]

    ctrl = subset[subset.obs.condition == 'control'].X[:, i]
    cov = subset[subset.obs.condition == 'COVID19'].X[:, i]
    print(f'{most_up} is the most upregulated gene expression:')
    print(f'Control: {ctrl.mean()}')
    print(f'COVID19: {cov.mean()}')
    #flip the rows from last to first, first to last
    dedf['log2fc'] = dedf['log2fc']*-1
    dedf = dedf.sort_values('log2fc', ascending=False).reset_index(drop=True)
    dedf = dedf[(dedf.qval < 0.05) & (abs(dedf.log2fc)>0.5)]
    up_down_regulated = dedf[-25:].gene.tolist() + dedf[:25].gene.tolist() #top 25 and bottom 25 from sorted

    # Generate heatmap for up/down-regulated genes
    heatmap_path_1 = os.path.join(output_folder, 'heatmap_condition.png')
    sc.pl.heatmap(subset, up_down_regulated, groupby='condition', swap_axes=True, save=heatmap_path_1)

    model = scvi.model.SCVI.load(model_file, adata)
    scvi_de_at = model.differential_expression(
        idx1 = [adata.obs['cell_type'] == 'AT1'],
        idx2 = [adata.obs['cell_type'] == 'AT2']
    )
    scvi_de_cond = model.differential_expression(
        idx1 = [(adata.obs['cell_type'].isin(['AT1','AT2'])) & (adata.obs.condition == 'COVID19')],
        idx2 = [(adata.obs['cell_type'].isin(['AT1','AT2'])) & (adata.obs.condition == 'control')]
    )
    scvi_de_at = scvi_de_at[(scvi_de_at['is_de_fdr_0.05']) & (abs(scvi_de_at.lfc_mean) > 0.5)]
    scvi_de_at = scvi_de_at.sort_values('lfc_mean')
    scvi_de_cond = scvi_de_cond[(scvi_de_cond['is_de_fdr_0.05']) & (abs(scvi_de_cond.lfc_mean) > 0.5)]
    scvi_de_cond = scvi_de_cond.sort_values('lfc_mean')
    scvi_de_at = scvi_de_at[(scvi_de_at.raw_normalized_mean1 > 0.5) | (scvi_de_at.raw_normalized_mean2 > 0.5)]
    scvi_de_cond = scvi_de_cond[(scvi_de_cond.raw_normalized_mean1 > 0.5) | (scvi_de_cond.raw_normalized_mean2 > 0.5)]
    up_down_regulated1 = scvi_de_at[-25:].index.tolist() + scvi_de_at[:25].index.tolist()
    up_down_regulated2 = scvi_de_cond[-25:].index.tolist() + scvi_de_cond[:25].index.tolist()

    heatmap_path_2 = os.path.join(output_folder, 'heatmap_cell_type.png')
    sc.pl.heatmap(subset, up_down_regulated1, groupby='cell_type', swap_axes=True, layer='scvi_normalized', log=True, save=heatmap_path_2)

    heatmap_path_3 = os.path.join(output_folder, 'heatmap_condition_cell_type.png')
    sc.pl.heatmap(subset, up_down_regulated2, groupby='condition', swap_axes=True, layer='scvi_normalized', log=True, save=heatmap_path_3)

    temp = subset[subset.obs.cell_type == 'AT2']

    
    return f"scRNAseq Analysis concluded. Images saved to 'output' folder"

# Download, extract, and decompress files
ftp_url = "ftp.ncbi.nlm.nih.gov"
ftp_directory = "/geo/series/GSE171nnn/GSE171524/suppl/"
filename = "GSE171524_RAW.tar"
local_filename = "GSE171524_RAW.tar"
extract_folder = "GSE171524_RAW"

def main():
    # Step 1: Run the functions
    download_ftp_file(ftp_url, ftp_directory, filename, local_filename)
    extract_tar_file(local_filename, extract_folder)
    decompress_gz_files(extract_folder)

    output_folder = 'output'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)   
    
    # Step 2: Run parameter adjustments
    parameter_adjust(extract_folder)

    # Step 3: Run integration to process the decompressed files
    integration(extract_folder, cutoff_value)
    gc.collect()

    # Step 4: Run processing
    processing('combined.h5ad')
    gc.collect()

    # Step 5: Run Analysis
    analysis('integrated.h5ad', 'model.model')

if __name__ == "__main__":
    main()
