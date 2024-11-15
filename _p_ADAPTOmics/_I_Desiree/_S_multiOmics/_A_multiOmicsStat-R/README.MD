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
