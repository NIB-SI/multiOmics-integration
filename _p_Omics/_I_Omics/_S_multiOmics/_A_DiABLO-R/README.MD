# Assay _A_DiABLO-R

## Protein selection

input/data_Proteomics.txt
- data


and short names from

input/Proteomics-subset-annotation.txt 
- proteinID	shortName	MapManBins	Use



## Stress and samples

HDW - ignore

C (Ctrl, control): 1, 7, 8, 14, 15, 21, 28
H: 1, 7, 8, 14
D: 8, 14
HD: 8, 14
W: 1, 7 (sever damage after days)


## Phenotyping

input/data_Phenomics.txt

DENOTE: 

1. water.consumption has missing values for day 1, therefore: i) for C, H, D, HD day 2 was taken; ii) for W day 3 was taken (only early day with non-negative values)
2. deltaT had negative values - al values shifted by abs(min())



## random notes


recovery after H: 15, 21, 28
recovery after D: 15, 21, 28
recovery after HD: 15, 21, 28
recovery after W: 15, 21, 28 (wobbly)


logFC
START
H1 vs C1
D8 vs C8
W1 vs C1
STOP
H7 vs C7
HD14 vs C14
W7 vs C7



