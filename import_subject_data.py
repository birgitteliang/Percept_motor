import os
import pandas as pd
import mne

def import_files(percept_ID):
    # Read the Excel file containing the data
    pathFiles = '/Users/birgittethomsen/Library/CloudStorage/OneDrive-Personal/Dokumenter/PhD/Charité/data/'
    df = pd.read_excel(os.path.join(pathFiles, 'percept_log.xlsx'))
    
    # Find the row corresponding to the given percept_ID
    subject_data = df.loc[df['percept_ID'] == percept_ID]

    if subject_data.empty:
        print(f"Percept_ID '{percept_ID}' not found in the data.")
        return None
    
    # Extract the file names for the subject
    ltpfile = subject_data['ltpfile'].iloc[0]
    accfile = subject_data['accfile'].iloc[0]
    epochsfile = subject_data['epochsfile'].iloc[0]
    hemisphere = subject_data['hemisphere'].iloc[0]

    # Read the data files based on the extracted file names
    raw = mne.io.read_raw_fieldtrip(os.path.join(pathFiles, ltpfile), info=None)
    acc = pd.read_csv(os.path.join(pathFiles, accfile), header=None)
    epochs = pd.read_excel(os.path.join(pathFiles, epochsfile))

    return raw, acc, epochs, hemisphere