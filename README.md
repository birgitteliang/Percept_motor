# Percept_motor
Plot rest and motor task recordings of Percept DBS data

percept_motor.ipynb is the working script.
To run the script you need to make a file "percept_log_template.xlsx" where you fill out the subject-id, hemisphere, and the names of the data files.
Epochs should be filled out in the provided epochs template file "SubX_epochs.xlsx"

The script calls the different script
1)import_subject_data
This script contains a function to import data.
It needs to be updated with your own path to your data.

2) dat_preproc
The script contains functions to analyze the percept data.
The script is written by Varvara Mathiopoulou.

3) accelerometer
Contains functions to analyze accelerometer data.

4) powerspectra
Contains functions to plot time-frequency plots and powerspectra.

Furthermore the script:
baseline_correction
is provided in the folder, but is not called from the percept_motor script.
