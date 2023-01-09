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
    * point-to-point
    * polynomial
  * extrapolation 
  * imputation
    * for qPCR see Baebler, Š., Svalina, M., Petek, M. et al. quantGenius: implementation of a decision support system for qPCR-based gene quantification. BMC Bioinformatics 18, 276 (2017). https://doi.org/10.1186/s12859-017-1688-7
  
3. Statistical analysis of individual omics data layers
   * ggplot {ggplot2}
   * corr.test {psych} - _Find the correlations, sample sizes, and probability values between elements of a matrix or data.frame_
   * cor.plot {psych} - _Create an image plot for a correlation or factor matrix_
   * pairs.panels {psych} - _SPLOM, histograms and correlations for a data matrix_
   * rcorr {Hmisc} - _Matrix of Correlations and P-values_
   * heatmaply_cor {heatmaply} - _Cluster heatmap based on plotly_
   * corrplot {corrplot} - _A visualization of a correlation matrix_
   * pheatmap {pheatmap} - _A function to draw clustered heatmaps_
   * t_test {rstatix}
   * ggdotplot {ggpubr}
   * metaMDS {vegan} - _Nonmetric Multidimensional Scaling with Stable Solution from Random Starts, Axis Scaling and Species Scores_
   * {limma}
     * limma::lmFit
     * limma::makeContrasts
     * limma::contrasts.fit
     * limma::eBayes
     * limma::decideTests
     * limma::topTable
4. Correlation based network inference within each omics level
   * Leave-One-Out graphs
     * qgraph {qgraph}
     * [igraph](https://igraph.org/r/)
5. Integration across different omics datasets
 * [Canonical Correlation Analysis](https://mixomics.org/methods/)
 * [N-Integration Discriminant Analysis with DIABLO](https://mixomics.org/mixDIABLO/)
 * Leave-One-Out graphs
   * qgraph {qgraph}
   * [igraph](https://igraph.org/r/)
6. Integration of data with prior knowledge
   * [igraph](https://igraph.org/r/)
   * [Cytoscape]https://cytoscape.org/)

