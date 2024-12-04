# Code HSPCdynamics
Code used in _Andrea Calabria_ et al, Nature 2024 ([DOI: 10.1038/s41586-024-08250-x](https://www.nature.com/articles/s41586-024-08250-x)).

This repository stores the code and data we used to analyze the clonal dynamics in vivo in patients under Lentiviral Hematopoietic Stem Cell Gene Therapy for the diseases:
- metachromatic leukodystrophy (MLD),
- Wiskott-Aldrich syndrome (WAS),
- β-thalassemia (β-Thal).

## Code
We organized the code in the folder [code](code) as follows: 
- Step 1: run ISAnalytics on raw matrixes to generate the first results per patient and per study. Folder [ISAnalytics](code/1.ISAnalytics_AnalysesPerStudy).
- Step 2: compare studies and patients starting from the results generated in Step 1. All scripts are available in the folder [Compararive Analyses](code/2.Comparative_Analyses). Within this folder, you will also find the code used for confounding factor removal (Bayesian model) and Good Turing model (folder [GoodTuring](code/2.Comparative_Analyses/GoodTuring)).
- [*Notebook*](code/Notebook) here we release the code of our Bayesian regression and Good Turing model, and analytical checks for regressions. HTML (and Rmd) files need to be downloaded for a preview.
- An additional folder [*utils*](code/utils) contains the functions required for the analysis of lineage ourput and commitment (per patient in percentage and at a single-clone resolution).

## Data
All matrixes per patient are available in the folder [data](data), grouped by Disease and then named with the prefix of each patient ID.
We also included the [Accuracy](data/Accuracy) folder containing all results from our bootstrap approach to test the robustness of our results.

## Contacts
For questions about code and data please write an email to Andrea Calabria (calabria.andrea AT hsr.it)
