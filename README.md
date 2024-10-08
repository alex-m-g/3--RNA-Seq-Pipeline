#  A molecular single-cell lung atlas of lethal COVID-19

## Description

Identify lung tissue responses to COVID-19 infection from 26 patients (19 individuals who died from COVID-19, 7 control individuals).

## Dataset

Status:	Public on Apr 29, 2021

Title: Columbia University/NYP COVID-19 Lung Atlas

Organism: Homo sapiens

Experiment type: Expression profiling by high throughput sequencing

Summary:	We profiled 116,314 cells using snRNA-seq of 20 frozen lungs obtained from 19 COVID-19 decedents and seven control patients with short postmortem interval (PMI) autopsies. The COVID-19 cohort comprises seven female and 12 male decedents, including 13 patients of Hispanic ethnicity, with an age range from 58 to >89 years, who had acquired SARS-CoV-2 infection and succumbed to the disease. The average time from symptom onset to death was 27.5 days (range, 4–63 days). After rapid autopsy with a median PMI of 4 hours (range 2–9 hours), collected tissues were either flash-frozen or frozen following OCT (optimal cutting temperature) embedment and subjected to snRNA-seq using a droplet-based platform (10x Genomics). All included patients had underlying hypertensive disorder and frequently one or more additional co-morbidities associated with increased risk for severe COVID-19.
  	
Overall design:	Single-nuclei RNA sequencing of 116,314 cells from 20 frozen lungs obtained from 19 COVID-19 decedents and seven control patients.

Citation: Melms JC, Biermann J, Huang H, Wang Y, et al. A molecular single-cell lung atlas of lethal COVID-19. Nature 2021 Jul;595(7865):114-119.

PMID: 33915568
## Pipeline
[scRNA Sequence Analysis Pipeline_v3](./scRNAseq_Analysis_Pipeline_v3.py)
## Files
1.  [scRNA_seq_Data_Initialize](./scRNA_seq_Data_Analysis_COIVD19.ipynb): Processed a single file to fine-tune the parameters for integration.
2.  [Integration](./scRNA_seq_Analysis_Integration.ipynb): Integrated all patient samples. Output file "combined.h5ad"
3.  [Model Training & Cell Type Identification](./scRNA_seq_Analysis_Cell_Marker_Labeling.ipynb): Processed the "combined.h5ad" file to train SCVI model, clustering, and identify cell types. Output file "integrated.h5ad"
4.  [Visual Analysis Development](./scRNA_seq_Analysis.ipynb): Opened the "integrated.h5ad" file to analyze and develop visualizations identify the effects of COVID-19.

## Results and Analysis
### Cell Marker Identification
![Cell_Marker_Identification](https://github.com/user-attachments/assets/c9137684-ab6a-4fa2-9cdf-2cfc9552ef31)
![Cell_Marker_Identification](https://github.com/user-attachments/assets/c9e2b2e1-e5be-4aff-b992-64be62f32e5a)

### Cell Type Frequency in Control vs COVID Patients
![fig1](https://github.com/user-attachments/assets/c79f2ed6-3817-408c-bcd5-21d08577288e)

### Gene Expression Heatmap in AT1 vs. AT2 (Differential Expression)
![differential_expression_heatmap](https://github.com/user-attachments/assets/ea8b6e44-12f7-429a-ac9c-658fc81bb282)

### Gene Expression Heatmap in COVID-19 Patients vs. Control Patients (Differential Expression)
![heatmap_2](https://github.com/user-attachments/assets/1d04347b-7a1a-4737-89da-479f323a3d44)


### Gene Ontology Enrichment
This violin plot compares the different expressions of gene 'ETV5' in COVID-19 patients and Control patients. The p-value is listed above. 

![violin_plot1](https://github.com/user-attachments/assets/7d304b9d-76a0-4838-8671-77703bb853c2)

## Accomplishments:
- Pre-processing clustering with one sample
- Integration of multiple samples
- Found marker genes
- Plotted said marker genes
- Labeling cells
- Counted the fraction of cells
- Differential expression between two different groups of cells
- Made differential expression heat maps
- Gene ontology
- Keg enrichment
- Comparisons between genes in different conditions
- Statistical testing
- Scored cells based on their expression of gene signatures... and plotted it

## Skills:
- **Programming Languages:** Python (Pandas, Numpy)
- **Data Processing and Analysis:** ScanPy, scVI, SOLO
- **Machine Learning & Deep Learning:** Variational Autoencoders (VAE), scVI model training, doublet detection using SOLO.
- **Bioinformatics:** Gene filtering, Quality control metrics calculations, Differential Expression Analysis.
- **Data Visualization:** Plotting gene markers, differential expression heatmaps.

