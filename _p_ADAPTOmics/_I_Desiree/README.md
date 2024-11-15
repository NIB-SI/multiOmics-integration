_I_Desiree
========

In this directory you can find:

- <ins>*_Phenodata_*</ins> - Master sample description file: 01_Experimental-design-and-days-of-tissue-sampling.xlsx
- <ins>*_Featuredata_*</ins> - Phenomics variable description file: 02_Phenomics-featuredata.xlsx
- RGB phenotipisation e.g.: 03_Phenomics-RBG_Hyponasty.pdf

  
Analysis steps are grouped in Assays within a joint Study

___  
<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/presentations/potato.jpg" width="300" class="center">



This protocol describes multi-omics data integration and modelling and systems biology analysis and visualization pipeline in R and Python

***Keywords***:  multiOmics, multi-omics, integrative omics, panOmics, pan-omics

Prior to data analysis and integration, some prerequisites should be met

## Prerequisites and protocol

### Data management framework
- All data should be annotated with detailed metadata, including <ins>*_Phenodata_*</ins> – a master sample description table, and preferably structured according to the pISA-tree data management framework (1)
- Each measurement should be stored as a separate data file, including **SampleName** (and/or **SampleID**) and measurements only, preferably as tab-separated text file
- If using Excel as primary measurement storage format, prior to analyses, export each data set as tab separated text file
- Avoid merged cells and Excel calculations in exported files under any cost

### Expected measurements 
- Omics' strategies include: Hormonomics, Transcriptomics, Proteomics (non-targeted), Metabolomics (leaves and tubers), and Phenomics
- Measurements can include single or multiple genotypes with single or multiple tissues under single or combine abiotic or biotic stressors
- Preferential measurements are time-series measurements

### Analysis steps 
#### Step 1: Exprimental design master table
- Design <ins>*_Phenodata_*</ins>, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice
- Store <ins>*_Phenodata_*</ins> at **_Investigation_** level
- Define relative path of <ins>*_Phenodata_*</ins> in [_INVESTIGATION_METADATA.TXT](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_INVESTIGATION_METADATA.TXT), as well as in <ins>*__ASSAY_METADATA.TXT_*</ins>
- <ins>*_Phenodata_*</ins> must contain **SampleName** (and/or **SampleID**) column, which will be utilized to combine measurements with sample descriptions

#### Step 2: Data preprocessing and overall inspection 
- Prior to analyses, it is expected that the analyst conducted data preprocessing and overall inspection, which might include: 
  * i) detection of outliers and faulty measurements,
  * ii) data transformation,
  * iii) interpolation,
  * iv) extrapolation and,
  * v) imputation
- For qPCR imputation suggestions see (2)
- For other steps see suggested packages in the [README.md](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/README.md) file of this repository

#### Step 3: Statistical analysis of individual omics data layers
- This step focuses on various within-level correlations calculations and visualization, calculation of correlation differences, multidimensional scaling, t-tests and log2FC calculations, variable/feature selection and more
- Example of input data, consisting of various Omics’ levels can be found within **_Assay_** **input** directory
- Experimental design can be inspected from [01_Experimental-design-and-days-of-tissue-sampling.xlsx](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/01_Experimental-design-and-days-of-tissue-sampling.xlsx) file, stored at **_Investigation_** level
- In this scenario, SampleName column was created from condition, time point, and plant replicate number
- Plants were exposed to Heat stress, Drought stress, Waterlogging stress, double- and triple-stress (H, D, W, HD, HDW; later defined as ‘Stress’); measurements considered in this analysis were taken at days: 1, 7, 8, and 14; which is denoted under **Treatment** and **SamplingDay** columns of <ins>*_Phenodata_*</ins>
- For statistical analysis of individual omics data layers prepare data in similar manner, and use scripts from **_Assay_** [_A_multiOmicsStat-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmicsStat-R) and **_Assay_** [_A_multiOmics-FS-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py)
- Main packages and functions are listed in the [README.md](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/README.md) file of this repository

#### Step 4: Correlation based network inference within each omics level 
- Correlation-based network can be constructed using ‘Leave-One-Out’ or Lioness (3) approach from all Omics’ levels
- Since they are highly depended on thresholds we advise to use thresholding approach from **_Assay_** [_A_multiOmics-differential-networks-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-differential-networks-Py) 

#### Step 5: Integration across different omics datasets
- This step focuses on Canonical Correlation Analysis and N-Integration Discriminant Analysis with DIABLO (4) using mixOmics (5)
- Prepare data in a similar manner and run master script from **_Assay_** [_A_DiABLO-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_DiABLO-R)
- Main packages and functions are listed in the [README.md](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/README.md) file of this repository
- Read more about DIABLO at <https://mixomics.org/mixDIABLO/>

#### Step 6: Integration of data with prior knowledge 
- Log2FC values from Step 3, or standard differential expression analysis, can be visualised in prior knowledge network in a differential network context using Cytoscape and [DiNAR](https://github.com/NIB-SI/DiNAR)
- Cystoscape manual is available at <https://manual.cytoscape.org/en/stable/>
- Example of manually created network for this specific data set can be found within the ***input*** directory of the **_Assay_** [_A_multiOmics-differential-networks-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-differential-networks-Py)
- Visualised network are available from **_Assay_** [_A_multiOmics-visualisation-Python](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-visualisation-Python)
- For more prior knowledge network for plant species see <https://skm.nib.si/> and <https://knetminer.com/>

___
(1) Petek, M., Zagorščak, M., Blejec, A. et al. pISA-tree - a data management framework for life science research projects using a standardised directory tree. Sci Data 9, 685 (2022). https://doi.org/10.1038/s41597-022-01805-5

(2) Baebler, Š., Svalina, M., Petek, M. et al. quantGenius: implementation of a decision support system for qPCR-based gene quantification. BMC Bioinformatics 18, 276 (2017). https://doi.org/10.1186/s12859-017-1688-7

(3) Kuijjer ML, Hsieh PH, Quackenbush J, Glass K. lionessR: single sample network inference in R. BMC Cancer. 2019 Oct 25;19(1):1003. doi: 10.1186/s12885-019-6235-7

(4) Amrit Singh, Casey P Shannon, Benoît Gautier, Florian Rohart, Michaël Vacher, Scott J Tebbutt, Kim-Anh Lê Cao, DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays, Bioinformatics, Volume 35, Issue 17, September 2019, Pages 3055–3062, https://doi.org/10.1093/bioinformatics/bty1054

(5) Rohart F, Gautier B, Singh A, Lê Cao KA. mixOmics: An R package for 'omics feature selection and multiple data integration. PLoS Comput Biol. 2017 Nov 3;13(11):e1005752. doi: 10.1371/journal.pcbi.1005752. PMID: 29099853; PMCID: PMC5687754.
