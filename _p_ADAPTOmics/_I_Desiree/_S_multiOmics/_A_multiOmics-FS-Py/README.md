Assay _A_multiOmics-FS-Py
==============

***Feature Selection in Python with Scikit-Learn***


- Variable selection was conducted on the variable set from multiple imaging sensors

- The random forest (RF) algorithm from the R package caret v6.0-94 ([_A_multiOmicsStat-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmicsStat-R)) as well as the python package scikit-learn v1.2.0 were used with default settings, as RF showed the best performance out of a selection of algorithms

- Recursive feature elimination was applied in R and multiple importance scores, including mutual information (MI), Anova, RF importance and SHAP values (REF) were computed in Python, showing consistencies between the approaches for the top 5 variables

- The sixth variable, was selected based on expert knowledge

# model comparison

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/model comparison.svg" height="400">

# feature importance scores

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/feature_importance_all.svg" height="400">

# average feature importance
<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-FS-Py/reports/avg-drop_cols.svg" height="300">
