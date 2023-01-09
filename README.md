# **‘Omics data analysis pipeline**
- Data integration and modelling
- Systems biology analysis and visualization pipeline

## Data management framework
- [pISA-tree](https://github.com/NIB-SI/pISA-tree)
  * Petek, M., Zagorščak, M., Blejec, A. et al. pISA-tree - a data management framework for life science research projects using a standardised directory tree. Sci Data 9, 685 (2022). https://doi.org/10.1038/s41597-022-01805-5

## Expected measurements
- one or multiple genotypes
- under single and multiple abiotic/biotic stressors
- experiment duration: XY hours, days, ... : time-series experimental design
- tissue: single or multiple
- Omics' strategies: 
  * Hormonomics
  * Transcriptomics 
  * Proteomics (non-targeted)
  * Metabolomics 
  * Phenomics

## Analysis steps:
1. Design _Phenodata_, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice
2. Data preprocessing and overall inspection
  * detection of outliers and faulty measurements
    * [{gplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html)
    * [{rgl}](https://cran.r-project.org/web/packages/rgl/index.html)
  * data transformation (if needed)
    * [{multimode}](https://cran.r-project.org/web/packages/multimode/index.html)
    * [{fitdistrplus}](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    * [{caret}](https://cran.r-project.org/web/packages/caret/index.html)
    * [{glmnet}](https://cran.r-project.org/web/packages/glmnet/index.html)
    * [{MASS}](https://cran.r-project.org/web/packages/MASS/index.html)
    * [{BIGL}](https://cran.r-project.org/web/packages/BIGL/index.html)
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

