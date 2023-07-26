import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne

# FFT_epoch
# Finds time-frequency (fast fourier transformation) for one epoch
# The shape of the result is two-dimensional
# Inputs
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) Sxx: FFT output (np array (3 dimensiontal - stn x frequency x time))
# 3) hemisphere: The hemisphere of the STN electrode of interest (string) - Right: 'R', Left: 'L'
# 4) start: Start time of epoch (integer) e.g. 20
# 5) end: End time of epoch (integer) e.g. 30
    # e.g. FFT_epoch(Sxx, 'L', 20, 30)
def FFT_epoch(subid, Sxx, hemisphere, start, end):
    if hemisphere == 'L':
        side = 0
    elif hemisphere == 'R':
        side = 1
    
    FFT = Sxx[side, :, start:end]

    fig, axs = plt.subplots(figsize=(8, 6))
    axs.pcolormesh(FFT, cmap='viridis', vmin=0, vmax=0.2)
    axs.set_title(subid)
    axs.set_ylim([0, 120])  # Set y-axis limit to 0-120
    im = axs.pcolormesh(FFT, cmap='viridis', vmin=0, vmax=0.2)
    fig.colorbar(im)
    plt.tight_layout()
    plt.savefig(subid + ' Time-frequency plot', dpi = 300)
    plt.show()

    return FFT
    
# CALCULATE_POWERSPECTRUM
# Calculate the powerspectrum for one epoch
# Inputs
# FFT: The fft of one epoch (ndarray (2 dimensional - frequncy x time))
# The dimension of the result is one-dimensional (frequency)
    #e.g. calculate_powerspectrum(rs_off_fft)
def calculate_powerspectrum(FFT):
    return np.mean(FFT, 1)

# CALCULATE_MEAN_PS_BLOCKS
# Calculate the mean power spectrum of several epochs e.g. tapping blocks
# FFT: The fft of each epoch (ndarray (2 dimensional - frequncy x time))
# The dimension of the result is one-dimensional (frequency)
    #e.g. calculate_mean_ps_blocks(tap_off_1_fft,tap_off_2_fft)
def calculate_mean_ps_blocks(*FFTs):
    block_all = np.array(FFTs)
    block_mean = np.mean(block_all, axis=0)
    block_ps = np.mean(block_mean, 1)
    return block_ps

# PLOT_ALL_FFTS
# Plots time-frequency plots (fast fourier transformation) for all defined epochs
# Create a new folder to save the figures
# Inputs
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) Sxx: FFT output (np array (3 dimensiontal - stn x frequency x time))
# 3) hemisphere: The hemisphere of the STN electrode of interest (string) - Right: 'R', Left: 'L'
# 4) epochs: dataframe with times and names of epochs (size (2,10))
    #e.g. plot_all_FFTs('Sub999',Sxx,'L',epochs)
def plot_all_FFTs(subid, Sxx, hemisphere, epochs):
    from powerspectra import FFT_epoch

    # Create subject folder
    folder_name = subid
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Create a dictionary to store the FFTs of all epochs
    ffts = {}

    # Iterate through the columns of the 'epochs' DataFrame
    for column_name, column_values in epochs.iteritems():
        # Extract the start and stop times from the first and second rows
        start_time = column_values[0]
        stop_time = column_values[1]

        # Calculate the FFT using the extracted start and stop times
        FFT = FFT_epoch(subid, Sxx, hemisphere, start_time, stop_time)

        # Plot Time-frequency plot
        fig, axs = plt.subplots(1, 1, figsize=(8, 6))
        im = axs.pcolormesh(FFT, cmap='viridis', vmin=0, vmax=0.2)
        axs.set_title(column_name)
        axs.set_ylim([0, 120])  # Set y-axis limit to 0-120
        fig.colorbar(im)

        plt.tight_layout()

        # Save the figure in new subject folder
        fig.savefig(os.path.join(folder_name, f'{subid}_{column_name}.png'), dpi=300)
        plt.close(fig)  # Close the figure to free up memory

        # Store the FFT in the ffts dictionary with the epoch name as the key
        ffts[column_name] = FFT

    return ffts

# PLOT_POWER
# Plots powerspectra for all epochs
# Makes two subplots. One in the beta range (5-35 Hz) and one in the gamma range (60-90)
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) Sxx: FFT output (np array (3 dimensiontal - stn x frequency x time))
# 3) hemisphere: The hemisphere of the STN electrode of interest (string) - Right: 'R', Left: 'L'
# 4) epochs: dataframe with times and names of epochs (size (2,10))
    #e.g. plot_power('Sub999',tap_off_ps)
from powerspectra import plot_all_FFTs, calculate_powerspectrum, calculate_mean_ps_blocks

def plot_power(subid, Sxx, hemisphere, epochs):
    # Make dictionary with FFTs of all epochs
    ffts = plot_all_FFTs(subid, Sxx, hemisphere, epochs)
    items = list(ffts.items())
    ffts_keys = ffts.keys()
    # Create a new dictionary to store the powerspectra
    powerspectra_all = {}

    # Calculate powerspectra of the resting state blocks
    for key, value in ffts.items():
        if key.startswith('rs_'):
            # Create the new variable name with the suffix '_ps'
            new_key = key + '_ps'

            # Perform the calculation using the existing variable
            new_value = calculate_powerspectrum(value)

            # Store the new key-value pair in the new dictionary
            powerspectra_all[new_key] = new_value

    # Find epochs for tapping blocks
    tap_off_keys = [item for item in ffts.keys() if 'tap_off' in item]
    tap_off_values = [ffts[item] for item in tap_off_keys]
    tap_clin_keys = [item for item in ffts.keys() if 'tap_clin' in item]
    tap_clin_values = [ffts[item] for item in tap_clin_keys]
    tap_thres_keys = [item for item in ffts.keys() if 'tap_thres' in item]
    tap_thres_values = [ffts[item] for item in tap_thres_keys]

    # Calculate powerspectra of tapping blocks
    powerspectra_all['tap_off_ps'] = calculate_mean_ps_blocks(*tap_off_values)
    powerspectra_all['tap_clin_ps'] = calculate_mean_ps_blocks(*tap_clin_values)
    powerspectra_all['tap_thres_ps'] = calculate_mean_ps_blocks(*tap_thres_values)

    # Plot the powerspectra in beta and gamma ranges
    plt.figure(figsize=(10, 12))

    # Define subplot specifications
    subplots = [
        (3, 2, 1, {'xlim': (5, 35), 'ylim': (0, 3), 'ylabel': '0 mA', 'title': 'beta'}),
        (3, 2, 2, {'xlim': (60, 90), 'ylim': (0, 0.12), 'title': 'gamma', 'legend': True}),
        (3, 2, 3, {'xlim': (5, 35), 'ylim': (0, 3), 'ylabel': 'Clinical-0.5 mA'}),
        (3, 2, 4, {'xlim': (60, 90), 'ylim': (0, 0.12)}),
        (3, 2, 5, {'xlim': (5, 35), 'ylim': (0, 3), 'xlabel': 'Frequency [Hz]', 'ylabel': 'Threshold'}),
        (3, 2, 6, {'xlim': (60, 90), 'ylim': (0, 0.12), 'xlabel': 'Frequency [Hz]'})
    ]

    for i, (rows, cols, idx, settings) in enumerate(subplots, 1):
        plt.subplot(rows, cols, idx)
        key_rest = ['rs_off_ps', 'rs_off_ps', 'rs_clin_ps', 'rs_clin_ps', 'rs_thres_ps', 'rs_thres_ps'][i-1]
        key_tap = ['tap_off_ps', 'tap_off_ps', 'tap_clin_ps', 'tap_clin_ps', 'tap_thres_ps', 'tap_thres_ps'][i-1]
        plt.plot(np.arange(0, 126), powerspectra_all[key_rest], color='blue', label='Rest')
        plt.plot(np.arange(0, 126), powerspectra_all[key_tap], color='red', label='Tapping')

        plt.xlim(settings['xlim'])
        plt.ylim(settings['ylim'])
        if 'xlabel' in settings:
            plt.xlabel(settings['xlabel'])
        if 'ylabel' in settings:
            plt.ylabel(settings['ylabel'])
        if 'title' in settings:
            plt.title(settings['title'])
        if 'legend' in settings and settings['legend']:
            plt.legend(loc='upper right')

    plt.suptitle(subid)
    plt.tight_layout()
    plt.show()

    # Save the figure in new subject folder
    folder_name = subid
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    plt.savefig(os.path.join(folder_name, f'{subid}_powerspectra.png'), dpi=300)
    

