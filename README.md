This repository houses the code for the academic paper "Personalised antimicrobial susceptibility testing with clinical prediction modelling informs appropriate antibiotic use", for the purpose of peer review and subsequent open-sourcing.

If you use this code please cite this repository.

Instructions for use:

The source data can be obtained from PhysioNet at https://physionet.org/content/mimiciv/2.2/ provided the terms of access are met. The csv filenames used in this code match the default filenames that can be downloaded from the 'hosp' folder.

Before running the code, the data and Python files should be saved into a local directory - the filepath of this directory should then be substituted into the file in place of #FILEPATH# before running the code. The required package versions are included in the packages.csv file within this directory.

1. To reproduce the study, PDAST_1.R script must be run first.

Then, for the subsequent analyses:

To reproduce the main analysis
2. Run UDAST_LR.py
3. Run PDAST_2.R*

To reproduce the analysis where organism identification is available
2. Run UDAST_LR2.py

To reproduce the analysis where other AST results are also available
2. Run UDAST_LR3.py

To reproduce the analysis where all 'I' results are reclassified as 'R'
2. Run UDAST_LR4.py
3. Run PDAST_2B.R*

To reproduce the multinomial analysis
2. Run UDAST_LR5.py

To run the prototype application
2. Run app.R, using "session_urines.csv" as file upload when prompted

*This script uses reticulate to run the python script 'Prediction_run.py' within R, but the script can alternatively be run within python by running code up to and including creation of the 'daily_urines.csv' file, running the Prediction_run.py file in Python, then reading the "probs_df_overall.csv" file into R as the "probs_df_overall" object, then running the rest of the code.



