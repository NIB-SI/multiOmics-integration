# Assay _A_DiABLO-R

## Input

### Protein selection
```
../input/data_Proteomics.txt - data
../input/Proteomics-subset-annotation.txt
```
use protein short names

## Targeted omics
```
../input/data_hormonomics.txt
../input/data_metabolomics.txt
../input/data_qPCR.txt
```

## Phenotyping
```
../input/data_Phenomics.txt
```
4 out of 6 replicates selected



## Sample description
```
../input/phenodata_subset_2023-03-08.txt
```
- HDW - ignore
- C (Ctrl, control): 1, 7, 8, 14, 15, 21, 28
- H: 1, 7, 8, 14
- D: 8, 14
- HD: 8, 14
- W: 1, 7 (severe damage in the second weeks of waterlogging)







# Notes
## sample description
- recovery after H: 15, 21, 28
- recovery after D: 15, 21, 28
- recovery after HD: 15, 21, 28
- recovery after W: 15, 21, 28 (wobbly)

## technical notes
1. for engineered variables that have missing values, due to missing reference at for e.g. initial stage, take consecutive day containing non-negative values
2. for engineered variables that have negative values, and they should not - shift all values by abs(min())


===

The background network is merged with thresholded ([_A_multiOmics-differential-networks-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-differential-networks-Py)) correlation networks and visualised using Cytoscape ([_A_multiOmics-visualisation-Py](https://github.com/NIB-SI/multiOmics-integration/tree/main/_p_ADAPTOmics/_I_Desiree/_S_multiOmics/_A_multiOmics-visualisation-Py))
