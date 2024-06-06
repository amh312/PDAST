This repository houses the code for the academic paper **"Personalised antimicrobial susceptibility testing with clinical prediction modelling informs appropriate antibiotic use"**, for the purpose of peer review and subsequent open-sourcing.

If you use this code please cite this repository.

***Instructions for use:***

The source data can be obtained from PhysioNet at https://physionet.org/content/mimiciv/2.2/ once the terms of access are met. The csv filenames used in this code match the following default filenames that can be downloaded from the *hosp* folder at the bottom of the page: *"prescriptions.csv", "diagnoses_icd.csv", "procedures_icd.csv", "labevents.csv", "d_labitems.csv", "poe_detail.csv", "poe.csv", "omr.csv", "admissions.csv", "patients.csv"*, and *"services.csv"*.

*PhysioNet MIMIC-IV citations:
Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark, R. (2023). MIMIC-IV (version 2.2). PhysioNet. https://doi.org/10.13026/6mm1-ek67.

Johnson, A.E.W., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data 10, 1 (2023). https://doi.org/10.1038/s41597-022-01899-x

Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.*

This code was written and run using *R* version 4.3.2 and *Python* version 3.12.0, on a laptop computer running macOS Sonoma version 14.5 with an Apple M1 Pro processor, 16GB random-access memory and 10 cores. The code may need to be run in chunks, depending on application memory. The typical run time of all code was approximately 3-4 hours.

Before running the code, the data and *Python* files should be saved into a local directory - the filepath of this directory should then be substituted into the file in place of **#FILEPATH#** in all R scripts before running the code. The required package versions are included in the *packages.csv* file within this directory. A conda environment was used to run the *Reticulate* interface package - the local environment used should be substituted for **#CONDAENV_FILEPATH#** in the **PDAST_2.R**, **PDAST_2B.R**, and **app.R** scripts.

***Reproducing the study***

1. To reproduce the study, the **PDAST_1.R** script and the **Imports & functions.py** script must be run first.

Then, for the subsequent analyses:

To reproduce the main analysis, fairness analysis, and stability analysis: 
2. Run **UDAST_LR.py**; 
3. Run **PDAST_2.R***

To reproduce the analysis where organism identification is available: 
2. Run **UDAST_LR2.py**

To reproduce the analysis where other AST results are also available: 
2. Run **UDAST_LR3.py**

To reproduce the analysis where all 'I' results are reclassified as 'R': 
2. Run **UDAST_LR4.py**; 
3. Run **PDAST_2B.R***

To reproduce the multinomial analysis: 
2. Run **UDAST_LR5.py**

To run the prototype application: 
2. Run **app.R**, using *"session_urines.csv"* as file upload when prompted.

*This script uses the *Reticulate* interface package to run the *Python* script **Prediction_run.py** within *R*, but the script can alternatively be run within *Python* by running code up to and including creation of the *"daily_urines.csv"* file, running the **Prediction_run.py** file in *Python*, then reading the *"probs_df_overall.csv"* file into *R* as the **probs_df_overall** object, then running the rest of the code without the *Reticulate* code lines.

Please be aware that randomisation processes are used several times throughout the code - results may therefore vary slightly from those presented in the manuscript.

Small simulated datasets have also been provided to test run the code if required - these data are synthetic so will not reproduce the results of the main analysis, and will not run for segments of the code that are designed to be run on large datasets (i.e., the fairness analysis and stability analysis).
