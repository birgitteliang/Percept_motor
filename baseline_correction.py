import os
import numpy as np
import re
import pandas as pd
from matplotlib import pyplot as plt
import mne
from mne.baseline import rescale
from scipy.signal import hann
from importlib import reload
import dat_preproc

#BL_CORRECT
#Calculates the baseline corrected time-frequency plot.
#Resting state with the same stimulation amplitude is used as baseline.
#Every second of tapping is corrected for the averages power of resting state.
#Inputs 1)tapping time-frequency data (two-dimensional shape) 2)resting state power spectrum (one dimensional shape) 3)Title of your plot
    #e.g. blcorrect(tap_off_fft, rs_off_ps, 'Sub999_Baseline_corrected_0mA')

def bl_correct(subid,):
    Sxx_norm = np.empty(shape=(126, 10))
    Sxx_norm[:] = np.nan

    for j in range(tap.shape[1]):
        Sxx_thisnorm = ((tap[:, j] - rs) / rs) * 100
        Sxx_norm[:, j] = Sxx_thisnorm

    plt.pcolormesh(Sxx_norm, cmap='bwr', vmin=-300, vmax=300)
    plt.xlabel('Time[s]')
    plt.ylabel('Frequency [Hz]')
    plt.title(title)
    cbar = plt.colorbar()
    cbar.set_label('Percent change')
    plt.show()