# **‘multiOmics data analysis, integration, and visualisation protocol**
- Data integration and modelling with [R](https://cran.r-project.org/)
- Systems biology analysis and visualization pipeline in [R](https://cran.r-project.org/) and [Python](https://www.python.org/)

## Data management framework
- [pISA-tree on GitHub](https://github.com/NIB-SI/pISA-tree)
  * Petek, M., Zagorščak, M., Blejec, A. et al. pISA-tree - a data management framework for life science research projects using a standardised directory tree. Sci Data 9, 685 (2022). https://doi.org/10.1038/s41597-022-01805-5

## Expected measurements
- one or multiple genotypes
- under single and multiple abiotic/biotic stressors
- experiment duration: hours, days, ... : time-series experimental design
- tissue: single or multiple
- Omics' strategies: 
  * Hormonomics
  * Transcriptomics 
  * Proteomics
  * Metabolomics 
  * Phenomics

## Analysis steps:
1. Design _Phenodata_, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice
2. Data preprocessing and overall inspection
  * detection of outliers and faulty measurements
  * data transformation (if needed)
  * interpolation
    * point-to-point
    * polynomial
  * extrapolation  
    <img src="https://www.statology.org/wp-content/uploads/2021/09/interp3-768x545.png" width=25% height=25%>
  * imputation

3. ***Statistical analysis of individual omics data layers***
4. ***Correlation based network inference within each omics level***
   * Leave-One-Out graphs
   * Lioness
   
5. ***Integration across different omics datasets***
6. ***Integration of data with prior knowledge***
   
## Start of the analysis:
- Data is expected to be arranged within data management framework, with complete and descriptive metadata files, including _Phenodata_ file. 
- 'Omics files are expected to be preprocessed (see suggestions in Step 2). 
- Minimal input files can be found within './input' directory. 
- For Step 3: Statistical analysis of individual omics data layers run script [01_Step3.Rmd](https://github.com/NIB-SI/multiOmics-integration/blob/main/_I_Omics/_S_multiOmics/_A_multiOmics-integration-R/scripts/01_Step3.Rmd)
- For Step 4: Correlation based network inference within/between each omics level run script [02_Step4.Rmd](https://github.com/NIB-SI/multiOmics-integration/blob/main/_I_Omics/_S_multiOmics/_A_multiOmics-integration-R/scripts/02_Step4.Rmd)
- For Step 5: Integration across different omics datasets run script 03_Step5.Rnw

For more info see [multiOmics_data_analysis_Protocol](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_Omics/multiOmics_data_analysis_Protocol.docx)
