from matplotlib import pyplot as plt
import mne
import scipy
import numpy as np
#from scipy.fft import fft, fftfreq
from scipy.signal import spectrogram, hann, butter, filtfilt

#DATA FILTERING

#A. LOW/HIGH PASS FILTER
def low_highpass_filter(data, frequency_cutoff_low, frequency_cutoff_high):
        """
        low_highpass_filter applies low (100Hz) and high pass (5Hz) filter, with a 
        5th order butterworth filter.

        input <-- data:
        np.array of shape channels x freqs x time e.g. data = raw.get_data(picks=[0,1])

        output <-- pass_filtered_dat:
        np.array with the same shape as input array

        requires: scipy.signal
        """
        filter_order = 5 
        #frequency_cutoff_low = 5 
        #frequency_cutoff_high = 100 
        fs = 250 
            
        # create the filter
        b, a = scipy.signal.butter(filter_order, (frequency_cutoff_low, frequency_cutoff_high), 
                btype='bandpass', output='ba', fs=fs)              
        pass_filtered_dat = scipy.signal.filtfilt(b, a, data) 
    
        return pass_filtered_dat

#B. BANDSTOP FILTER
def bandstop_filter(lowcut, highcut, data):
        """
        bandstop_filter applies a bandstop filter around the low/highcut given

        input <-- data:
        np.array of shape channels x freqs x time e.g. data = raw.get_data(picks=[0,1])

        output <-- stop_filtered_dat:
        np.array with the same shape as input array

        requires: scipy.signal
        """      
        order = 4
        nyq = 0.5 * 250 #sampling rate
        low = lowcut / nyq
        high = highcut / nyq
        
        #create the filter
        b, a = scipy.signal.butter(order, [low, high], btype='bandstop')
        
        stop_filtered_dat = scipy.signal.filtfilt(b,a,data)

        return stop_filtered_dat

#FFT Transformation and Spectrogram Plotting
def fft_rawviz(raw, x, win_samp, noverlap):
        """
        fft_rawviz performs a Fast Fourier Transformation to data and creates TF plots 
        with stimulation amplitude on top

        # Input:
        #x = filt_dat
        #win_samp = window for fft in samples, e.g. 250 for 1 sec
        #noverlap e.g. 0.25 (for 25%)
        """

        fs = 250
        window = hann(win_samp, sym=False)
        f, t, Sxx = scipy.signal.spectrogram(x = x, fs = fs, window = window, noverlap = noverlap)
         
        #Plot Spectrograms of both STNs
        fig, axes = plt.subplots(1,2, figsize = (18,6))
        fig.suptitle('FFT Transformations')

        ax_c = 0
        stim = 4
        for kj in np.array([0,1]):
                
                ax2 = axes[kj].twinx() #make right axis linked to the left one
                if kj == 1:
                        stim_data = (raw.get_data(picks = stim)[0,:]) #define stim channel
                elif kj == 0:
                        stim_data = (raw.get_data(picks = stim)[0,:])
                
                #Plot LFP data
                axes[ax_c].specgram(x = x[kj,:], Fs = fs, noverlap = noverlap, cmap = 'viridis',
                        vmin = -25, vmax = 10)
                axes[ax_c].set_ylim(bottom = 5,top = 100)
                axes[ax_c].set_xlim(0,raw.n_times/250)
                
                amplitude = stim_data/3
                #Plot stim channel on top
                ax2.plot(raw.times, amplitude, 'w', linewidth = 1.5)
                ax2.set_yticks(np.arange(0,4,0.25))

                #Right y axis label only for second plot to avoid crowd
                if kj == 1:
                        ax2.set_ylabel('Stimulation Amplitude [mA]')
                
                axes[ax_c].set_ylabel('Frequency [Hz]')
                axes[ax_c].set_xlabel('Time [sec]')
                axes[ax_c].set_title(raw.ch_names[kj])

                ax_c += 1
                stim += 1

        
        plt.show(block = False)

        #np.save(new_fname, Sxx)
        return f, t, Sxx


def fft_transform(x, win_samp, noverlap):
        """
        fft_transform performs a Fast Fourier Transformation to data without plotting them

        input:
        - x = filt_dat as np.array of shape channel x freq x time
        - win_samp = window for fft in samples, e.g. 250 for 1 sec
        - noverlap e.g. 0.25 (for 25%)

        output:
        - f (frequencies), t (time), Sxx: transformed data of shape same as input data
        """
        fs = 250
        window = hann(win_samp, sym=False)
        
        f, t, Sxx = scipy.signal.spectrogram(x = x, fs = fs, window = window, noverlap = noverlap)

        plt.specgram(x = x, Fs = fs, noverlap = noverlap, cmap = 'viridis',
                        vmin = -25, vmax = 10)
        
        plt.ylim(5,100)
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')

        return f, t, Sxx

#EPOCH_PS
def epoch_ps(filt_dat, epoch_df, window, noverlap, side, ylim2, title):
        """A dummy docstring."""
        ps = np.empty([epoch_df.shape[0],126])

        for index, row in epoch_df.iterrows():
                this_onset = epoch_df.onset[index]*250
                this_offset = this_onset + (epoch_df.duration[index]*250)
                
                #window = hann(250, sym=False)

                ff, Pxx = scipy.signal.welch(filt_dat[side,this_onset:this_offset], fs = 250, 
                        nperseg = window, noverlap = noverlap)
                
                #xnew = np.linspace(ff.min(), ff.max(), 50)
                #spl = make_interp_spline(ff, Pxx, k = 3)
                #y_smooth = spl(xnew)
                
                
                #plt.plot(xnew, y_smooth, label = key)
                plt.plot(ff, Pxx, label = epoch_df.description[index])
                ps[index] = Pxx
        
        plt.xlim([50, 100])
        plt.ylim([0, ylim2])

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [$V^{2}$/Hz]')
        plt.legend()
        plt.title(title)
        plt.show(block = False)

        return ps

from scipy.interpolate import make_interp_spline, BSpline
from scipy import stats

def mypower(ps):
        
        x = np.arange(1,127)
        xvals = np.linspace(1,127,1250)

        if ps.ndim > 2:
                y = np.mean(ps,1)
                spl = make_interp_spline(x,y, k=3)  # type: BSpline
                power_smooth = spl(xvals)
                sem = stats.sem(ps,1)
                spl_sem = make_interp_spline(x,sem,k=3)  # type: BSpline
                sem_smooth = spl_sem(xvals) 

                plt.plot(xvals, power_smooth)
                plt.fill_between(xvals, power_smooth - sem_smooth, power_smooth + sem_smooth)
        else:
                y = ps
                spl = make_interp_spline(x,y, k=3)  # type: BSpline
                power_smooth = spl(xvals)
                plt.plot(xvals, power_smooth)

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Spectral Power')

# FFT_EPOCHS
# FFT Transformation and Spectrogram Plotting
# Plotting of epochs on top of the time-frequncy (fast fourier transportion) of the full data set combined with sitm data.
# Inputs
# raw=raw data
# x=filtered data
# win_samp = window for fft in samples, e.g. 250 for 1 sec
# noverlap = 0.25 (for 25%)
# epochs = epoch data
# hemisphere = Right: 'R', Left: 'L'
def fft_epochs(raw, x, win_samp, noverlap, epochs, hemisphere):
        """
        fft_rawviz performs a Fast Fourier Transformation to data and creates TF plots 
        with stimulation amplitude on top

        # Input:
        #x = filt_dat
        #win_samp = window for fft in samples, e.g. 250 for 1 sec
        #noverlap e.g. 0.25 (for 25%)
        """

        fs = 250
        window = hann(win_samp, sym=False)
        f, t, Sxx = scipy.signal.spectrogram(x = x, fs = fs, window = window, noverlap = noverlap)
         
        #Plot Spectrograms of both STNs
        fig, axes = plt.subplots(1,2, figsize = (18,6))
        fig.suptitle('FFT Transformations')

        ax_c = 0
        stim = 4
        for kj in np.array([0,1]):
                
                ax2 = axes[kj].twinx() #make right axis linked to the left one
                if kj == 1:
                        stim_data = (raw.get_data(picks = stim)[0,:]) #define stim channel
                elif kj == 0:
                        stim_data = (raw.get_data(picks = stim)[0,:])
                
                #Plot LFP data
                axes[ax_c].specgram(x = x[kj,:], Fs = fs, noverlap = noverlap, cmap = 'viridis',
                        vmin = -25, vmax = 10)
                axes[ax_c].set_ylim(bottom = 5,top = 100)
                axes[ax_c].set_xlim(0,raw.n_times/250)
                
                amplitude = stim_data/3
                #Plot stim channel on top
                ax2.plot(raw.times, amplitude, 'w', linewidth = 1.5)
                ax2.set_yticks(np.arange(0,4,0.25))

                #Right y axis label only for second plot to avoid crowd
                if kj == 1:
                        ax2.set_ylabel('Stimulation Amplitude [mA]')
                
                axes[ax_c].set_ylabel('Frequency [Hz]')
                axes[ax_c].set_xlabel('Time [sec]')
                axes[ax_c].set_title(raw.ch_names[kj])

                ax_c += 1
                stim += 1

        if hemisphere == 'L':
                side = 0
        elif hemisphere == 'R':
                side = 1
        for epoch_name, epoch_range in epochs.items():
                axes[side].axvspan(epoch_range[0], epoch_range[1], alpha=0.3, color='red')
        
        plt.show(block = False)

        #np.save(new_fname, Sxx)
        return f, t, Sxx