This repository houses the code for the academic paper **"Personalised antimicrobial susceptibility testing with clinical prediction modelling informs appropriate antibiotic use"**, for the purpose of peer review and subsequent open-sourcing.

If you use this code please cite this repository.

***Instructions for use:***

The source data can be obtained from PhysioNet at https://physionet.org/content/mimiciv/2.2/ once the terms of access are met. The csv filenames used in this code match the following default filenames that can be downloaded from the *hosp* folder at the bottom of the page: *"prescriptions.csv", "diagnoses_icd.csv", "procedures_icd.csv", "labevents.csv", "microbiologyevents.csv", "d_labitems.csv", "poe_detail.csv", "poe.csv", "omr.csv", "admissions.csv", "patients.csv"*, and *"services.csv"*.

PhysioNet MIMIC-IV citations:

*Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark, R. (2023). MIMIC-IV (version 2.2). PhysioNet. https://doi.org/10.13026/6mm1-ek67.*

*Johnson, A.E.W., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data 10, 1 (2023). https://doi.org/10.1038/s41597-022-01899-x*

*Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.*

This code was written and run using *R* version 4.3.2 and *Python* version 3.12.0, on a laptop computer running macOS Sonoma version 14.5 with an Apple M1 Pro processor, 16GB random-access memory and 10 cores. The code may need to be run in chunks, depending on application memory. The typical run time of all code was approximately 3-4 hours.

Before running the code, the data and *Python* files should be saved into a secure local directory - the filepath of this directory should then be substituted into the file in place of **#FILEPATH#** in all R scripts before running the code. The required package versions are included in the *packages.csv* file within this directory. A conda environment was used to run the *Reticulate* interface package - the local environment used should be substituted for **#CONDAENV_FILEPATH#** in the **PDAST_2.R**, **PDAST_2B.R**, and **app.R** scripts.

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

Please be aware that randomisation processes are used several times throughout the code - results may therefore vary slightly from those presented in the manuscript.

***Testing the code***

Given the complexity of the analysis, we recommend the process described above to test the full code. However, we have simulated a small synthetic dataset in order to enable the core elements of the model and study design code to be run quickly if required. To run the code:

1. Save the "*Urines5c.csv"*, *"Urines_assess.csv"*, *"patients.csv"* *"microbiologyevents.csv"* synthetic datasets under these filenames in a local directory
2. Run **Imports & functions.py** and **Packages & functions.R**, populating the **#FILEPATH#** and **#CONDAENV_FILEPATH#** aspects as described above
3. Run **UDAST_LR2.py**
4. Run **PDAST_2.R***

The synthetic data is random except for a variable correlation between prior resistance and a resistant result, simulated to yield a credible level of predictive performance for each agent - running code on the synthetic data will not, however, reproduce the results of the main analysis.


*This script uses the *Reticulate* interface package to run the *Python* script **Prediction_run.py** within *R*, but the code can alternatively be run by running **PDAST_2C.R** (instead of running **PDAST_2.R** or **PDAST_2B.R**), running the **Prediction_run.py** file in *Python*, then running **PDAST_3.R**. Remember to substitute in your local directory filepath in place of #FILEPATH#.
