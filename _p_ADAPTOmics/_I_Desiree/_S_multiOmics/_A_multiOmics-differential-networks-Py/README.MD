# Assay multiOmics-differential-networks-Py


***Integration across different omics datasets***

- Correlations between components measured in various Omics’ levels is calculated and visualised using DIABLO ([_A_DiABLO-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_DiABLO-R)) as implemented in the mixOmics package

- The correlation matrix is calculated separately for each stress as well as for control


***Automated graph thresholding***

- Absolute value thresholding is sequentially applied (from 0 to 2, in increments of 0.01)

- At each increment, edges below the threshold are removed, and subsequently all isolated vertices also removed

- At each threshold increment the network structure is studied

- Finally, per graph thresholds are defined as the largest threshold resulting in a local minima of the density, and applied to obtain thresholded networks

***Integration of data with prior knowledge***

- The background network is merged with thresholded correlation networks and visualised using Cytoscape (![_A_multiOmics-visualisation-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-visualisation-Py))
- LogFC values are taken from [_A_multiOmicsStat-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmicsStat-R)
