# 	Data and Code for 'Common non-antibiotic drugs enhance selection for antimicrobial resistance in mixture with ciprofloxacin'

- DOI of the paper on biorxiv: https://doi.org/10.1101/2025.03.23.644574
- Sequence data is deposited at ENA - Accession Number PRJEB88784.
- All other data is deposited here as part of this repository and deposited on Zenodo ![image](https://github.com/user-attachments/assets/439bfb2c-e9b4-4e53-9c0a-0d0861c82a4f)


## Outline

This repository contains the datasets and code, analysed as part of the study in the manuscript 'Common non-antibiotic drugs enhance selection for antimicrobial resistance in mixture with ciprofloxacin'. 

It contains the growth OD readings, qPCR assay gene quantities, and processed metagenome sequence data (directly from the AMR++ output). It also contains the code used to analyse these data and produce all figures in this paper. The code for running the bioinformatics pipelines are not present here, but can be run with default parameters. In this study we used AMR++ to analyse the metagenome data. AMR++ is well supported by a GitHub repository - https://github.com/Microbial-Ecology-Group/AMRplusplus. 

### Notes on the code 
All code is written to be run in separate R projects. A project is associated with a group of folders, and all figures, scripts and data are kept within each folder. Personally, I keep my growth data in a folder called ‘Growth Analyses’, with a nested set of folders – data, figures, and scripts. An R project will be associated with these folders. Only the code and data are deposited here, but these can be altered to suit how you work best. 

