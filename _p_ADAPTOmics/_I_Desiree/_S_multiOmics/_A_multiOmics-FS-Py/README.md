Assay _A_multiOmics-FS-Py
==============

***Feature Selection in Python with Scikit-Learn***


- Variable selection was conducted on the non-invasive phenomics variable sets 

- The random forest (RF) algorithm from the R package caret v6.0-94 (![_A_multiOmicsStat-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_Omics/_I_Omics/_S_multiOmics/_A_multiOmicsStat-R)) as well as the python package scikit-learn v1.2.0 were used with default settings, as RF showed the best performance out of a selection of algorithms

- Recursive feature elimination was applied in R and multiple importance scores, including mutual information (MI), Anova, RF importance and SHAP values (REF) were computed in Python, showing consistencies between the approaches for the top 5 variables

- The sixth variable, was selected based on expert knowledge

# model comparison

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/model_comparison.png" height="400">

# RF feature importance scores

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/RF_importance_all.png" height="400">
# RF feature improtance, MI

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/RF_importance_mutualInformation.png" height="300">