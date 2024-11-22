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
- W: 1, 7 (sever damage in later days)



DENOTE: 

1. water.consumption has missing values for day 1, therefore: i) for C, H, D, HD day 2 was taken; ii) for W day 3 was taken (only early day with non-negative values)
2. deltaT had negative values - al values shifted by abs(min())



# Notes
- recovery after H: 15, 21, 28
- recovery after D: 15, 21, 28
- recovery after HD: 15, 21, 28
- recovery after W: 15, 21, 28 (wobbly)





