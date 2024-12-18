# Assay multiOmics-differential-networks-Py


***Integration across different omics datasets***

- Correlations between components measured in various Omics’ levels is calculated and visualised using DIABLO ([_A_DiABLO-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_DiABLO-R)) as implemented in the mixOmics package

- The correlation matrix is calculated separately for each stress as well as for control


***Automated graph thresholding***

- Absolute value thresholding is sequentially applied (from 0 to 2, in increments of 0.01)

- At each increment, edges below the threshold are removed, and subsequently all isolated vertices also removed

- At each threshold increment the network structure is studied

- Finally, per graph thresholds are defined as the largest threshold resulting in a local minima of the density, and applied to obtain thresholded networks

<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-differential-networks-Py/other/threshold_e.g.png" width="600">

***Integration of data with prior knowledge***

- The background network is merged with thresholded correlation networks and visualised using Cytoscape (![_A_multiOmics-visualisation-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-visualisation-Py))

<p>
<img src="https://github.com/NIB-SI/multiOmics-integration/blob/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-differential-networks-Py/other/e.g.svg" width="600">

<em> e.g. HD-C; line type: control: dots, treatment: equal dash, treatment-control: solid, blue: negative CC diff, red: postitive CC diff </em>
</p>


- LogFC values are taken from [_A_multiOmicsStat-R](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmicsStat-R). For phenomics data relative change % (compared to control) was calculated.

