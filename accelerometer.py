import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import mne

# ACC_PLOT
# Plots accelerometer data of the three dimensions and saves the plot in the subject folder.
# Input
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) acc: accelerometer data (dataframe (2 dimensions axis x samples))
    #e.g. acc('Sub999',acc)
def acc_plot(subid, acc):
    time = acc.shape[1]
    pivot_table = acc.pivot_table(index=None, columns =['x','y','z'])
    pivot_table.insert(0, 'Time',range(0,time))

    #Plot accelerometer data
    #250 samples per second therefore divided by 250
    plt.figure(figsize=(10, 12))
    plt.subplot(3, 1, 1)
    plt.plot(np.array(pivot_table['Time']) / 250, pivot_table['x'], '.-', linewidth=0.5, markersize=0.5, color='blue')
    plt.title('Accelerometer data')
    plt.ylabel('X acceleration')

    plt.subplot(3, 1, 2)
    plt.plot(np.array(pivot_table['Time']) / 250, pivot_table['y'], '.-', linewidth=0.5, markersize=0.5, color='red')
    plt.ylabel('Y acceleration')

    plt.subplot(3, 1, 3)
    plt.plot(np.array(pivot_table['Time']) / 250, pivot_table['z'], '.-', linewidth=0.5, markersize=0.5, color='green')
    plt.xlabel('Time (s)')
    plt.ylabel('Z acceleration')

    plt.show()

    #Save figure 
    folder_name = subid
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    plt.savefig(os.path.join(folder_name, f'{subid}_Accelerometer_data.png'), dpi=300)

# ACC_STIM
# Plots the z accelerometer data and the stimulation amplitude on top and saves the plot in the subject folder.
# Can be used to find times of epochs of interest.
# Input
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) raw: the raw percept data (mne object (rawarray))
# 3) acc: accelerometer data (dataframe (2 dimensions axis x measures))
# 4) hemisphere: The hemisphere of the STN electrode of interest (string) - Right: 'R', Left: 'L'
    #e.g. acc_stim('Sub999',raw,acc,'L')
def acc_stim(subid, raw, acc, hemisphere):
    if hemisphere == 'L':
        side = 4
    elif hemisphere == 'R':
        side = 5
    time = acc.shape[1]
    pivot_table = acc.pivot_table(index=None, columns =['x','y','z'])
    pivot_table.insert(0, 'Time',range(0,time))
    
    #Stim channel data
    y = raw.get_data(reject_by_annotation = 'omit',picks=[side])
    # Create the figure and axes
    fig, ax1 = plt.subplots()
    # Set the x-axis label
    ax1.set_xlabel('Time (s)')
    # Plot the z acceleration on the left y-axis
    ax1.set_ylabel('Z acceleration')
    ax1.plot(np.array(pivot_table['Time']) / 250, pivot_table['z'], '.-', linewidth=0.5, markersize=0.5, color='blue')
    # Create a twin axis with a different scale for the right y-axis
    ax2 = ax1.twinx()
    # Plot the stimulation current on the right y-axis
    ax2.set_ylabel('Stimulation current (mA)')
    ax2.plot(raw.times, np.squeeze(y)/3, 'w', linewidth = 1.5, color='red')
    plt.title(subid + ' Accelerometer and Stimulation Amplitude')
    plt.show()

    #Save figure 
    folder_name = subid
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    plt.savefig(os.path.join(folder_name, f'{subid}_Accelerometer and Stimulation Amplitude.png'), dpi=300)

# ACC_EPOCHS
# Creates a plot of the z accelerometer data, the stimulation amplitude, and the epochs.
# and saves the plot in the subject folder.
# Input
# 1) subid: The subject id (string) e.g. "Sub-999"
# 2) raw: the raw percept data (mne object (rawarray))
# 3) acc: accelerometer data (dataframe (2 dimensions axis x measures))
# 4) epochs: dataframe with times and names of epochs (size (2,10))
# 5) hemisphere: The hemisphere of the STN electrode of interest (string) - Right: 'R', Left: 'L'
    #e.g. acc_epochs('Sub999',raw, acc, epochs, 'L')
def acc_epochs(subid, raw, acc, epochs, hemisphere):
    if hemisphere == 'L':
        side = 4
    elif hemisphere == 'R':
        side = 5
    time = acc.shape[1]
    pivot_table = acc.pivot_table(index=None, columns =['x','y','z'])
    pivot_table.insert(0, 'Time',range(0,time))

    #Stim channel data
    y = raw.get_data(reject_by_annotation = 'omit',picks=[side])
    # Create the figure and axes
    fig, ax1 = plt.subplots()
    # Set the x-axis label
    ax1.set_xlabel('Time (s)')

    # Plot the z acceleration on the left y-axis
    ax1.set_ylabel('Z acceleration')
    ax1.plot(np.array(pivot_table['Time']) / 250, pivot_table['z'], '.-', linewidth=0.5, markersize=0.5, color='blue')

    # Create a twin axis with a different scale for the right y-axis
    ax2 = ax1.twinx()

    # Plot the stimulation current on the right y-axis
    ax2.set_ylabel('Stimulation current (mA)')
    ax2.plot(raw.times, np.squeeze(y) / 3, 'w', linewidth=1.5, color='red')

    for epoch_name, epoch_range in epochs.items():
        ax1.axvspan(epoch_range[0], epoch_range[1], alpha=0.3, color='gray')

    # Show the plot
    plt.title(subid + ' Epochs')
    plt.show()

    #Save figure 
    folder_name = subid
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    plt.savefig(os.path.join(folder_name, f'{subid}_Accelerometer and Epochs.png'), dpi=300)