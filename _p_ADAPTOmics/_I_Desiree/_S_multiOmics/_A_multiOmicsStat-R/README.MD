Assay _A_multiOmicsStat-R
=========================

***Proteomics***
- ![DEP](https://bioconductor.org/packages/release/bioc/html/DEP.html)
- ![GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

*Denote:* GSEA (2.1) produces multiple files per comparison. Combine information (2.2) and visualise (2.3) in R

***Phenomics***
- Feature selection ([FS](https://github.com/topepo/caret/)) with Random Forest (RF) and Recursive Feature Elimination (RFE) {noninvasive} (3.1)
- Nonparametric statistics {noninvasive, invasive} (3.2, 3.3)
- Various plots {noninvasive, invasive} (3.2, 3.3, 3.4)

*Deote*: hyperspectral (3.4) measurements are not conducted on the same days as RGB, thermal, targeted Omics', etc, therefore are just visualised

***Transcriptomics, hormonomics, metabolomics***
- logFC (4.1)
- t-test: Welch's (4.1) and [Permutation](https://github.com/stamats/MKinfer) t-Test (4.2)


___
# Notes
## R and package versions
- manuscript versions:
    * R version 4.3.1
    * missForest v1.5
    * caret v6.0-94
    * DEP v 1.22.0 ([Bioconductor version 3.17](https://bioconductor.org/packages/3.17/BiocViews.html#___Software))
    * MKinfer v1.1
    * mixOmics v6.24.0 ([Bioconductor version 3.17](https://bioconductor.org/packages/3.17/BiocViews.html#___Software))
- currect versions:
    * R version 4.4.1
    * missForest v1.5
    * caret v6.0-94
    * DEP v 1.28.0 ([Bioconductor version 3.20](https://bioconductor.org/packages/3.20/BiocViews.html#___Software))
    * MKinfer v1.2
    * mixOmics v6.30.0 ([Bioconductor version 3.20](https://bioconductor.org/packages/3.20/BiocViews.html#___Software))

*Denote*: there can be some discrepancies between outputs from different versions

## List of useful R f-ctions and packages
### Data preprocessing and overall inspection
  #### detection of outliers and faulty measurements
  * [{gplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html)
  * [{rgl}](https://cran.r-project.org/web/packages/rgl/index.html)
  * [{raster}](https://cran.r-project.org/web/packages/raster/index.html)
  * [{multcompView}](https://cran.r-project.org/web/packages/multcompView/index.html)
  * [{gridExtra}](https://cran.r-project.org/web/packages/gridExtra/index.html)

  #### data transformation (if needed)
 - [{multimode}](https://cran.r-project.org/web/packages/multimode/index.html)
  - [{fitdistrplus}](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
  - [{caret}](https://cran.r-project.org/web/packages/caret/index.html)
  - [{glmnet}](https://cran.r-project.org/web/packages/glmnet/index.html)
  - [{MASS}](https://cran.r-project.org/web/packages/MASS/index.html)
  - [{BIGL}](https://cran.r-project.org/web/packages/BIGL/index.html)
  - [{robustbas}](https://cran.r-project.org/web/packages/robustbase/index.html)
  - [{preprocessCore}](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)
  - [{compositions}](https://cran.r-project.org/web/packages/compositions/index.html)
  - [{mgcv}](https://cran.r-project.org/web/packages/mgcv/index.html)

  #### interpolation
  ##### point-to-point
  * approxfun {stats} - _Returns a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation._
  ##### polynomial
  * predict {stats} - _A generic function for predictions from the results of various model fitting functions_
  * boxplot.stats {grDevices} - _Box Plot Statistics_
  * mad [{BiocGenerics}](https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html) - _Compute the median absolute deviation for a vector_
  * aq.plot [{mvoutlier}](https://cran.r-project.org/web/packages/mvoutlier/index.html) - _Adjusted quantile plots for multivariate outlier detection_
 
#### extrapolation  
<img src="https://www.statology.org/wp-content/uploads/2021/09/interp3-768x545.png" width=25% height=25%>

#### imputation
* for qPCR see Baebler, Š., Svalina, M., Petek, M. et al. quantGenius: implementation of a decision support system for qPCR-based gene quantification. BMC Bioinformatics 18, 276 (2017). https://doi.org/10.1186/s12859-017-1688-7
* for other levels check [NAguideR](https://github.com/wangshisheng/)
  
### ***Statistical analysis of individual omics data layers***
   * ggplot [{ggplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html) - various plots, <https://r-graphics.org/chapter-ggplot2>
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
     * limma::lmFit - _Linear Model for Series of Arrays_
     * limma::makeContrasts - _Construct Matrix of Custom Contrasts_
     * limma::contrasts.fit - _Compute Contrasts from Linear Model Fit_
     * limma::eBayes - _Empirical Bayes Statistics for Differential Expression_
     * limma::decideTests - _Multiple Testing Across Genes and Contrasts_
     * limma::topTable - _Table of Top Genes from Linear Model Fit_
### ***Correlation based network inference within each omics level***
   * Leave-One-Out graphs
     * qgraph [{qgraph}](https://cran.r-project.org/web/packages/qgraph/)
     * [igraph](https://igraph.org/r/)
   * Lioness
   
   *Denote:* Since results from both methods heavily depend on selected thresholds, Lioness node and edge selection using FDR being even more sensitive on correlation difference cut-off, we suggest to use an automated graph thresholding approach.

### ***Integration across different omics datasets***
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
### ***Integration of data with prior knowledge***
   * [igraph](https://igraph.org/r/)
   * [DiNAR subApps](https://github.com/NIB-SI/DiNAR/tree/master/subApps)
   * [Cytoscape](https://cytoscape.org/)
   * <https://skm.nib.si/>

