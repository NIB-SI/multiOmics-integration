Variable selection was conducted on the non-invasive phenomics variable sets. 

The random forest (RF) algorithm from the R package caret v6.0-94 (Kuhn, 2008) as well as the python package scikit-learn v1.2.0 were used with default settings, as RF showed the best performance out of a selection of algorithms .

Recursive feature elimination was applied in R and multiple importance scores, including mutual information, Anova, RF importance and SHAP values (REF) were computed in Python, showing consistencies between the approaches for the top 5 variables. The sixth variable, was selected based on expert knowledge.
