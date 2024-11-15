_I_Desiree
========

- Master sample description file: 01_Experimental-design-and-days-of-tissue-sampling.xlsx
- Phenomics variable description file: 02_Phenomics-featuredata.xlsx
- RGB phenotipisation e.g.: 03_Phenomics-RBG_Hyponasty.pdf
- Analysis steps are grouped in Assays within a joint Study

  
<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/presentations/potato.jpg" width="300" class="center">



This protocol describes multi-omics data integration and modelling and systems biology analysis and visualization pipeline in R and Python

***Keywords***:  multiOmics, multi-omics, integrative omics, panOmics, pan-omics

Prior to data analysis and integration, some prerequisites should be met

## Prerequisites and protocol

### Data management framework
All data should be annotated with detailed metadata, including Phenodata – a master sample description table, and preferably structured according to the pISA-tree data management framework (1). Each measurement should be stored as a separate data file, including SampleName (and/or SampleID) and measurements only, preferably as tab-separated text file. If using Excel as primary measurement storage format, prior to analyses, export each data set as tab separated text file. Avoid merged cells and Excel calculations in exported files under any cost.

### Expected measurements 
Omics' strategies include: Hormonomics, Transcriptomics, Proteomics (non-targeted),   Metabolomics, Phenomics and more. Measurements can include single or multiple genotypes with single or multiple tissues under single or combine abiotic or biotic stressors. Preferential measurements are time-series measurements.

### Analysis steps 
#### Step 1: Exprimental design master table
Design Phenodata, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice. Store Phenodata at Investigation level. Define relative path of Phenodata in _INVESTIGATION_METADATA.TXT, as well as in _ASSAY_METADATA.TXT. Phenodata must contain SampleName (and/or SampleID) column, which will be utilized to combine measurements with sample descriptions.
#### Step 2: Data preprocessing and overall inspection 
Prior to analyses, it is expected that the analyst conducted data preprocessing and overall inspection, which might include: i) detection of outliers and faulty measurements, ii) data transformation, iii) interpolation, iv) extrapolation and, v) imputation. For qPCR imputation suggestions see (2). For other steps see suggested packages in the README.md file of this repository.
#### Step 3: Statistical analysis of individual omics data layers
This step focuses on various within-level correlations calculations and visualization, calculation of correlation differences, multidimensional scaling, t-tests and log2FC calculations. Example of input data, consisting of three Omics’ levels can be found within Assay input directory. Experimental design can be inspected from 01_Experimental-design-and-days-of-tissue-sampling.xlsx file, stored at Investigation level. In this scenario, SampleName column was created from condition, time point, and plant replicate number. Plants were exposed to Heat stress (H, later defined as ‘Stress’); measurements were taken at days: 1, 7, 8, and 14; which is denoted under Treatment and SamplingDay columns of Phenodata. For statistical analysis of individual omics data layers prepare data in similar manner, and use scripts from Assay _A_multiOmicsStat-R. Main packages and functions are listed in the README.md file of this repository.
#### Step 4: Correlation based network inference within each omics level 
Correlation-based network can be constructed using ‘Leave-One-Out’ or Lioness approach from all Omics’ levels. Since they are highly depended on thresholds we advise to use thresholding approach from Assay _A_multiOmics-differential-networks-Py. 
#### Step 5: Integration across different omics datasets
This step focuses on Canonical Correlation Analysis and N-Integration Discriminant Analysis with DIABLO using mixOmics (3). Prepare data in a similar manner and run master script from Assay _A_DiABLO-R. Main packages and functions are listed in the README.md file of this repository. Read more about DIABLO at https://mixomics.org/mixDIABLO/.
#### Step 6: Integration of data with prior knowledge 
Log2FC values from Step 3, or standard differential expression analysis, can be visualised in prior knowledge network in a differential network context using Cytoscape. Cystoscape manual is available at https://manual.cytoscape.org/en/stable/. Example of manually created network for this specific data set can be found within the input directory of the Assay _A_multiOmics-differential-networks-Py. For more prior knowledge network for plant species see https://skm.nib.si/biomine/  and https://skm.nib.si/downloads/.
