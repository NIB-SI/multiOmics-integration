_S_multiOmics
=============

![](https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_Omics/_I_Omics/_S_multiOmics/reports/Pipeline.svg)

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/Pipeline.svg" width="600" class="center">

Run Assays in consecutive steps:

***_A_multiOmicsStat-R***
- Preprocessing
- Analysis of individual omics data layers  

***_A_multiOmics-FS-Py***
- Feature Selection in Python with Scikit-Learn
  
***_A_DiABLO-R***
- Integration across different omics datasets

***_A_multiOmics-differential-networks-Py***
- Automated graph thresholding
- Integration of data with prior knowledge

***_A_multiOmics-visualisation-Python***
- Visualisation of data with prior knowledge

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

## List of useful R f-ctions, packages, and web resources
### Data preprocessing and overall inspection
#### detection of outliers and faulty measurements
- Always plot the data before carrying on with any analysis
- Data visualisation tips at [Friends Don't Let Friends Make Bad Graphs](https://github.com/saeedsiddik/Graph_FriendsDontLetFriends)

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/density.png" height="200"> 
   
<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._NMDS.png" height="200"> <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._PCoA.png" height="200">   
  
- Check vignettes for: [{gplot2}](https://cran.r-project.org/web/packages/ggplot2/index.html),  [{rgl}](https://cran.r-project.org/web/packages/rgl/index.html), [{raster}](https://cran.r-project.org/web/packages/raster/index.html), [{multcompView}](https://cran.r-project.org/web/packages/multcompView/index.html), [{gridExtra}](https://cran.r-project.org/web/packages/gridExtra/index.html)

#### data transformation (if needed)
- Check prerequisites (e.g. distribution, homoscedasticity) to decide between parametric and nonparametric approach
- Check vignettes for: [{multimode}](https://cran.r-project.org/web/packages/multimode/index.html), [{fitdistrplus}](https://cran.r-project.org/web/packages/fitdistrplus/index.html), [{caret}](https://cran.r-project.org/web/packages/caret/index.html), [{glmnet}](https://cran.r-project.org/web/packages/glmnet/index.html), [{MASS}](https://cran.r-project.org/web/packages/MASS/index.html), [{BIGL}](https://cran.r-project.org/web/packages/BIGL/index.html), [{robustbas}](https://cran.r-project.org/web/packages/robustbase/index.html), [{preprocessCore}](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html), [{compositions}](https://cran.r-project.org/web/packages/compositions/index.html), [{mgcv}](https://cran.r-project.org/web/packages/mgcv/index.html)


<img src="https://www.statology.org/wp-content/uploads/2021/09/interp3-768x545.png" width=25% height=25%>  

  
#### interpolation

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._interpolate.png" height="250">  

##### point-to-point
* approxfun {stats} - _Returns a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation._
##### polynomial
* predict {stats} - _A generic function for predictions from the results of various model fitting functions_
* boxplot.stats {grDevices} - _Box Plot Statistics_
* mad [{BiocGenerics}](https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html) - _Compute the median absolute deviation for a vector_
* aq.plot [{mvoutlier}](https://cran.r-project.org/web/packages/mvoutlier/index.html) - _Adjusted quantile plots for multivariate outlier detection_

#### extrapolation

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._extrapolate.png" height="250"> <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._extrapolate_spline.png" height="250">     
  
*Denote:* extrapolate and interpolate only when you can be sure there is a trend, do not generate random or biased synthetic data
  

#### imputation
* for qPCR see Baebler, Š., Svalina, M., Petek, M. et al. quantGenius: implementation of a decision support system for qPCR-based gene quantification. BMC Bioinformatics 18, 276 (2017). https://doi.org/10.1186/s12859-017-1688-7
* for other levels check [NAguideR](https://github.com/wangshisheng/)

*Denote:* understand your data, is it Missing At Random (MAR) or Not Missing At Random (NMAR); do not introduce bias; define Limit of Detection (LOD) and Limit of Quantification (LOQ)
  
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

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._Loo.png" height="300">  

   * Leave-One-Out graphs
     * qgraph [{qgraph}](https://cran.r-project.org/web/packages/qgraph/)
     * [igraph](https://igraph.org/r/)    
   * Lioness [{lionessR}](https://bioconductor.org/packages/release/bioc/vignettes/lionessR/inst/doc/lionessR.html)

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._Lioness.png" height="300">   
   
   *Denote:* Since results from both methods heavily depend on selected thresholds, Lioness node and edge selection using FDR being even more sensitive on correlation difference cut-off, we suggest to use an automated graph thresholding approach.

### ***Integration across different omics datasets***
 * [Canonical Correlation Analysis](https://mixomics.org/methods/)
 * [N-Integration Discriminant Analysis with DIABLO](https://mixomics.org/mixDIABLO/)

   <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._diablo1.png" height="300">  

   <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/e.g._diablo2.png" height="300">  
   
     * block.splsda {mixOmics} - _N-integration and feature selection with Projection to Latent Structures models (PLS) with sparse Discriminant Analysis_
     * plotDiablo {mixOmics} - _Graphical output for the DIABLO framework_
     * plotVar {mixOmics} - _Plot of Variables_
     * plotIndiv {mixOmics} - _Plot of Individuals (Experimental Units)_
     * plotArrow {mixOmics} - _Arrow sample plot_
     * circosPlot {mixOmics} - _circosPlot for DIABLO_
     * cimDiablo {mixOmics} - _Clustered Image Maps (CIMs) ("heat maps") for DIABLO_
     * network {mixOmics} - _Relevance Network for (r)CCA and (s)PLS regression_
### ***Integration of data with prior knowledge and visualisation***
   * [Stress Knowledge Map](https://skm.nib.si/)

   <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/PSS-conceptual.png" height="400"> 
   
   * [SKM tools](https://github.com/NIB-SI/skm-tools)
   * [igraph](https://igraph.org/r/)
   * [DiNAR](https://github.com/NIB-SI/DiNAR)

   <img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/reports/dn.jpg" height="150"> 

   * [Cytoscape](https://cytoscape.org/)


