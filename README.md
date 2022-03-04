# GNPbaboons
Scripts and analyses performed in *Genomic variation in baboons from central Mozambique unveils complex evolutionary relationships with other Papio species*.

## mtDNA
Please refer to paper for exact methods used to perform BEAST divergence dating analysis using mtDNA, including jModel test. 
GNP mtDNA consensus sequence was constructed using ANGSD. Please see mtdna_ML.sh script for command.

## Y chromosome
Scripts relevant to Y-chromosome analyses can be found in the `ychrom` folder.

## Mapping and Variant Calling
Mapping and Variant Calling scripts can be found in their respective folders. 

## Population Structure
### *D*-Statistics and *f*-statistics
Scripts for both genotype-based *D*-statistics (POPSTAT) and read sampling method (ANGSD) can be found here. 
Details of *f<sub>3</sub>* and *f<sub>4</sub>* 

### PCA
Please see the paper's Material and Methods section for exact parameters used in PLINK and smartPCA.

## Population History
### TreeMix
Scripts for running TreeMix iteratively with different migration edges can be found in the `treemix` folder. 

### PSMC
Snakemake for PSMC can be found in the psmc folder. The `include.bed` file referred to in the Snakemake is the whitelist file (i.e. regions you WANT included). 

### qpGraph, ELAI, and MALDER/ALDER 
For details on how these methods were run, please refer to the paper's Material and Methods section. 
Should you require further clarification on set-up and running of these analyses please contact: Ludovica Molinaro (lu.molinario8@gmail.com) and Luca Pagani (lp.lucapagani@gmail.com).

Scripts and analyses performed in "Genomic variation in baboons from central Mozambique unveils complex evolutionary relationships with other Papio species."
