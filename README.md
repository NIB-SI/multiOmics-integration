# ADAPT_D5.3
**Potato ‘omics data analysis pipelines**

[ADAPT: Accelerated Development of multiple-stress tolerant potato [Horizon 2020]](https://adapt.univie.ac.at/)

- Deliverable No:	D5.3 
- WP5 Data integration and modelling
- Task 5.2.: Potato systems biology analysis and visualization pipelines

## Data management framework
- [pISA-tree](https://github.com/NIB-SI/pISA-tree)
  * Petek, M., Zagorščak, M., Blejec, A. et al. pISA-tree - a data management framework for life science research projects using a standardised directory tree. Sci Data 9, 685 (2022). https://doi.org/10.1038/s41597-022-01805-5

## Measurements
- cv. Desiree
- under single and multiple abiotic stressors
- experiment duration: 28 days
- tissue: leaves, tubers
- Omics' strategies: 
  * Hormonomics
  * Transcriptomics 
  * Proteomics (non-targeted)
    * representative pan-proteome from Petek, M., Zagorščak, M., Ramšak, Ž. et al. Cultivar-specific transcriptome and pan-transcriptome reconstruction of tetraploid potato. Sci Data 7, 249 (2020). https://doi.org/10.1038/s41597-020-00581-4
  * Metabolomics 
  * Phenomics

## Analysis steps:
1. Design _Phenodata_, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice
2. Data preprocessing and overall inspection
  * interpolation
  * extrapolation 
  * imputation
3. Statistical analysis of individual omics data layers
4. Correlation based network inference within each omics level
5. Integration across different omics datasets
6. Integration of data with prior knowledge
7. 
