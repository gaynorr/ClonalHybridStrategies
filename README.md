# ClonalHybridStrategies
This repository contains code and raw data from Labroo et al., 2022. https://link.springer.com/article/10.1007/s00122-023-04377-z

The code folder contains R files for the study code. Each file contains the code for a given scenario as defined in the study. The makeInitialPopulations.R file
generates the initial populations. The exact initial populations used in the study (for direct comparison of results) are available upon request. The PloidyInbreedingDepression.R file demonstrates that autopolyploids have more inbreeding depression than their corresponding diploids produced by genome reduction, as referenced in Supplemental Table 4.

The data folder contains zipped CSV files of the study raw data and a metadata file, "Metadata.csv", that refers to all CSVs except "calcLD.csv". Except "Metadata.csv" and "calcLD.csv", Each CSV file contains the data for a given scenario, and the column names and order of all CSV files are the same. The "calcLD.csv" file contains the initial mean and standard deviation of linkage disequilibrium (D) in the initial population.


