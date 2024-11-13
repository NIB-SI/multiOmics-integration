Assay _A_015_thresholding-py
==============

***Integration across different omics datasets***

- Correlations between components measured in various Omics’ levels were calculated and visualised using DIABLO (![_A_DiABLO-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_Omics/_I_Omics/_S_multiOmics/_A_DiABLO-R)) as implemented in the mixOmics package. 

- The correlation matrix was calculated separately for each stress as well as for control


***Automated graph thresholding***

- Absolute value thresholding was sequentially applied (from 0 to 2, in increments of 0.01)

- At each increment, edges below the threshold were removed, and subsequently all isolated vertices also removed

- At each threshold increment the network structure was studied

- Finally, per graph thresholds were defined as the largest threshold resulting in a local minima of the density, and applied to obtain thresholded networks

***Integration of data with prior knowledge***

- The background network was merged with thresholded correlation networks and visualised using Cytoscape (![_A_multiOmics-visualisation-Python](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_Omics/_I_Omics/_S_multiOmics/_A_multiOmics-visualisation-Python))
