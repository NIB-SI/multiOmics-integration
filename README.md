# **‘Omics data analysis pipeline**
- Data integration and modelling with [R](https://cran.r-project.org/)
- Systems biology analysis and visualization pipeline in [R](https://cran.r-project.org/)

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
  * Proteomics
  * Metabolomics 
  * Phenomics

## Analysis steps:
1. Design _Phenodata_, a master experimental design table describing samples for analysis, prior to sample collection according to good data management practice
2. Data preprocessing and overall inspection
  * detection of outliers and faulty measurements
    * [{gplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html)
    * [{rgl}](https://cran.r-project.org/web/packages/rgl/index.html)
    * [{raster}](https://cran.r-project.org/web/packages/raster/index.html)
    * [{multcompView}](https://cran.r-project.org/web/packages/multcompView/index.html)
    * [{gridExtra}](https://cran.r-project.org/web/packages/gridExtra/index.html)
  * data transformation (if needed)
    * [{multimode}](https://cran.r-project.org/web/packages/multimode/index.html)
    * [{fitdistrplus}](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    * [{caret}](https://cran.r-project.org/web/packages/caret/index.html)
    * [{glmnet}](https://cran.r-project.org/web/packages/glmnet/index.html)
    * [{MASS}](https://cran.r-project.org/web/packages/MASS/index.html)
    * [{BIGL}](https://cran.r-project.org/web/packages/BIGL/index.html)
    * [{robustbas}](https://cran.r-project.org/web/packages/robustbase/index.html)
    * [{preprocessCore}](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)
    * [{compositions}](https://cran.r-project.org/web/packages/compositions/index.html)
    * [{mgcv}](https://cran.r-project.org/web/packages/mgcv/index.html)
  * interpolation
    * point-to-point
      * approxfun {stats} - _Returns a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation._
    * polynomial
      * predict {stats} - _A generic function for predictions from the results of various model fitting functions_
      * boxplot.stats {grDevices} - _Box Plot Statistics_
      * mad [{BiocGenerics}](https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html) - _Compute the median absolute deviation for a vector_
      * aq.plot [{mvoutlier}](https://cran.r-project.org/web/packages/mvoutlier/index.html) - _Adjusted quantile plots for multivariate outlier detection _
  * extrapolation
  
    <img src="https://www.statology.org/wp-content/uploads/2021/09/interp3-768x545.png" width=25% height=25%>
  * imputation
    * for qPCR see Baebler, Š., Svalina, M., Petek, M. et al. quantGenius: implementation of a decision support system for qPCR-based gene quantification. BMC Bioinformatics 18, 276 (2017). https://doi.org/10.1186/s12859-017-1688-7
  
3. Statistical analysis of individual omics data layers
   * ggplot [{ggplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html)
   * corr.test [{psych}](https://cran.r-project.org/web/packages/psych/index.html) - _Find the correlations, sample sizes, and probability values between elements of a matrix or data.frame_
   * cor.plot [{psych}](https://cran.r-project.org/web/packages/psych/index.html) - _Create an image plot for a correlation or factor matrix_
   * pairs.panels [{psych}](https://cran.r-project.org/web/packages/psych/index.html) - _SPLOM, histograms and correlations for a data matrix_
   * rcorr [{Hmisc}](https://cran.r-project.org/web/packages/Hmisc/index.html) - _Matrix of Correlations and P-values_
   * heatmaply_cor [{heatmaply}](https://cran.r-project.org/web/packages/heatmaply/index.html) - _Cluster heatmap based on plotly_
   * corrplot [{corrplot}](https://cran.r-project.org/web/packages/corrplot/index.html) - _A visualization of a correlation matrix_
   * pheatmap [{pheatmap}](https://cran.r-project.org/web/packages/pheatmap/index.html) - _A function to draw clustered heatmaps_
   * t_test [{rstatix}](https://cran.r-project.org/web/packages/rstatix/index.html)
   * ggdotplot, ggviolin [{ggpubr}](https://cran.r-project.org/web/packages/ggpubr/index.html)
   * metaMDS [{vegan}](https://cran.r-project.org/web/packages/vegan/index.html) - _Nonmetric Multidimensional Scaling with Stable Solution from Random Starts, Axis Scaling and Species Scores_
   * [{limma}](https://bioconductor.org/packages/release/bioc/html/limma.html) for e.g. non-targeted Proteomics, RNA-seq, ..
     * limma::lmFit
     * limma::makeContrasts
     * limma::contrasts.fit
     * limma::eBayes
     * limma::decideTests
     * limma::topTable
4. Correlation based network inference within each omics level
   * Leave-One-Out graphs
     * qgraph [{qgraph}](https://cran.r-project.org/web/packages/qgraph/)
     * [igraph](https://igraph.org/r/)
5. Integration across different omics datasets
 * [Canonical Correlation Analysis](https://mixomics.org/methods/)
 * [N-Integration Discriminant Analysis with DIABLO](https://mixomics.org/mixDIABLO/)  
   * [{mixOmics}](http://www.bioconductor.org/packages/release/bioc/html/mixOmics.html)
     * block.splsda {mixOmics} _N-integration and feature selection with Projection to Latent Structures models (PLS) with sparse Discriminant Analysis_
     * plotDiablo {mixOmics} _Graphical output for the DIABLO framework_
     * plotVar {mixOmics} _Plot of Variables_
     * plotIndiv {mixOmics} _Plot of Individuals (Experimental Units)_
     * plotArrow {mixOmics} _Arrow sample plot_
     * circosPlot {mixOmics} _circosPlot for DIABLO_
     * cimDiablo {mixOmics} _Clustered Image Maps (CIMs) ("heat maps") for DIABLO_
     * network {mixOmics} _Relevance Network for (r)CCA and (s)PLS regression_
 * Leave-One-Out graphs
   * qgraph [{qgraph}](https://cran.r-project.org/web/packages/qgraph/)
   * [igraph](https://igraph.org/r/)
6. Integration of data with prior knowledge
   * [igraph](https://igraph.org/r/)
   * [DiNAR subApps](https://github.com/NIB-SI/DiNAR/tree/master/subApps)
   * [Cytoscape](https://cytoscape.org/)
   
## The beginning of the analysis:
- Data is expected to be arranged within data management framework, with complete and descriptive metadata files, including _Phenodata_ file. 
- 'Omics files are expected to be preprocessed (see suggestions in Step 2). 
- Minimal input files can be found within './input' directory. 
- For Step 3: Statistical analysis of individual omics data layers run script [01_Step3.Rmd](https://github.com/NIB-SI/multiOmics-integration/blob/main/_I_Omics/_S_multiOmics/_A_multiOmics-integration-R/scripts/01_Step3.Rmd)
- For Step 4: Correlation based network inference within/between each omics level run script [02_Step4.Rmd](https://github.com/NIB-SI/multiOmics-integration/blob/main/_I_Omics/_S_multiOmics/_A_multiOmics-integration-R/scripts/02_Step4.Rmd)
- For Step 5: Integration across different omics datasets run script 03_Step5.Rnw
