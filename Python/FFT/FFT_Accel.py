import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter
import scipy.fftpack


def pwc(x, N=30, Th=25):
    # x : signal in
    # N : window size in samples
    # Th : threshold above which we think it looks like a step function
    out = np.zeros(x.shape)
    xw = np.zeros(x.shape)

    for i in range(2, len(x) - 1):
        iL = list(range(max(0, i - N), i))  # for python 3
        iR = list(range(i, min(len(x), i + N)))

        mL = np.average(x[iL])
        mR = np.average(x[iR])

        if abs(mR - mL) > Th:
            vL = np.var(x[iL])
            vR = np.var(x[iR])

            out[iL] = out[iL] + mL / vL
            out[iR] = out[iR] + mR / vR

            xw[iL] = xw[iL] + 1 / vL
            xw[iR] = xw[iR] + 1 / vR
        else:
            v = np.var(x[iL + iR])
            out[iL + iR] = out[iL + iR] + (mL + mR) / (2 * v)
            xw[iL + iR] = xw[iL + iR] + 1 / v

    # fill in edges
    out[0] = out[1]
    out[-1] = out[-2]
    return out / xw


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
#**************************************
pos = 3
Fs = 300.0 # sampling rate
plot_time = 60
start_time = 30
#**************************************
Ts = 1/Fs # sampling interval
sample_skip = 0 # number of FFT steps to delete for removing DC
cnts_to_kg = 0.009313226
num_to_plot = int(plot_time * Fs)
start_plot = int(start_time * Fs)
#fileNameLC = r'C:\Users\patmcc\Desktop\TESTDATA_1029\LWP_BIN_SPACER_WEIGHT_POS5_SLOTB_1KG_25.2_20171029_300hz_2.csv'
#fileName = r'C:\Users\patmcc\Desktop\LWP Location\retake\LWP_BIN_SPACER_WEIGHT_POS5_SLOTC_1KG_25.2_20171026_100hz.csv'
#fileName = r'C:\Users\patmcc\Desktop\LWP Location\retake\cleaned\LWP_LOC_BIN12_1KG_POS3_WEIGHT_25.2G.csv'
fileNameLC = r'C:\Users\patmcc\Desktop\CSV_CLEANER\cleaned\Shelf7027_test17SS_LC.csv'
fileNameACCEL = r'C:\Users\patmcc\Desktop\CSV_CLEANER\cleaned\Shelf7027_test17SS_ACCEL.csv'
#fileName = r'C:\Users\patmcc\Desktop\LAB5\LWP_BIN_SPACER_WEIGHT_POS1_1KG_25.2_20171029_120hz_Lab5.csv'

# read csv file
data_lc1 = np.genfromtxt(fileNameLC, skip_header=1+start_plot, dtype = float, delimiter = ',', max_rows = num_to_plot, usecols=(1 + pos * 2)) # select column to process by 'usecols=n'
data_lc2 = np.genfromtxt(fileNameLC, skip_header=1+start_plot, dtype = float, delimiter = ',', max_rows = num_to_plot, usecols=(4 + pos * 2)) # select column to process by 'usecols=n'
data_ACCEL = np.genfromtxt(fileNameACCEL, skip_header=1+start_plot, dtype = float, delimiter = ',', max_rows = num_to_plot, usecols=17) # select column to process by 'usecols=n'
data_ACCEL = data_ACCEL - np.average(data_ACCEL[0:100])
data_lc1 = ((data_lc1 - np.average(data_lc1[0:100]))) * cnts_to_kg
data_lc2 = (data_lc2 - np.average(data_lc2[0:100])) * cnts_to_kg
data = (data_lc1 + data_lc2)
# data_filter = butter_lowpass_filter(data, 5, Fs, 12)
# data_filter = butter_lowpass_filter(data_filter, .5, Fs, 2)
data_filter = pwc(data)
accel_filter = pwc(data_ACCEL)

n = len(data) # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
frq = frq[range(n//2)] # one side frequency range

#window = signal.kaiser(n, beta=14)
#window = 1

Y= 20*np.log10(scipy.fft(data))
Y= Y[range(n//2)]

Y_filtered = 20*np.log10(scipy.fft(data_filter))
Y_filtered = Y_filtered[range(n//2)]

# Y_1= 20*np.log10(scipy.fft(data_lc1))
# Y_1= Y_1[range(n//2)]
#
# Y_2= 20*np.log10(scipy.fft(data_lc2))
# Y_2= Y_2[range(n//2)]

Y_1= 20*np.log10(scipy.fft(data_ACCEL))
Y_1= Y_1[range(n//2)]

dc_delete = np.arange(0,sample_skip,1)
Y = np.delete(Y, dc_delete)
Y_filtered = np.delete(Y_filtered, dc_delete)
frq = np.delete(frq, dc_delete)

t = np.arange(0,T,Ts) # time vector
fig, ax = plt.subplots(2, 2)
fig.canvas.set_window_title(fileNameLC)

ax[0,0].plot(t,data, linewidth=0.5) # plotting the amplitude in counts
ax[0,0].plot(t,data_filter, linewidth=1) # plotting the amplitude in grams
ax[0,0].set_xlabel('Time (s)')
ax[0,0].set_ylabel('Amplitude (grams)')
ax[1,0].plot(frq,abs(Y), linewidth=0.5) # plotting the spectrum
ax[1,0].plot(frq,abs(Y_filtered), linewidth=0.5,alpha=0.8) # plotting the spectrum
ax[1,0].set_xlabel('Freq (Hz)')
ax[1,0].set_ylabel('|Y(freq)| (dB)')
ax[0,0].set_ylim(-100,100)
ax[1,0].set_ylim(0,130)
ax[0,1].plot(t,data_ACCEL, linewidth=0.5) # plotting the amplitude in counts
ax[0,1].plot(t,accel_filter, linewidth=0.5) # plotting the amplitude in counts
ax[0,1].set_xlabel('Time (s)')
ax[0,1].set_ylabel('Amplitude (g)')
ax[1,1].plot(frq,abs(Y_1), linewidth=0.5) # plotting the spectrum
ax[1,1].set_xlabel('Freq (Hz)')
ax[1,1].set_ylabel('|Y(freq)| (dB)')
ax[0,1].set_ylim(-0.05,0.05)
# ax[1,1].set_ylim(0,130)
# ax[0,1].plot(t,data_lc1, linewidth=0.5) # plotting the amplitude in counts
# ax[0,1].plot(t,data_lc2, linewidth=0.5, alpha=0.8) # plotting the amplitude in grams
# ax[0,1].set_xlabel('Time (s)')
# ax[0,1].set_ylabel('Amplitude (g)')
# ax[1,1].plot(frq,abs(Y_1), linewidth=0.5) # plotting the spectrum
# ax[1,1].plot(frq,abs(Y_2), linewidth=0.5, alpha=0.8) # plotting the spectrum
# ax[1,1].set_xlabel('Freq (Hz)')
# ax[1,1].set_ylabel('|Y(freq)| (dB)')
# ax[0,1].set_ylim(-100,100)
# ax[1,1].set_ylim(0,130)

rms = np.sqrt(np.mean(data**2))
rms_filtered = np.sqrt(np.mean(data_filter**2))
print(str(round(rms, 2)) + ' grms')
print(str(round(rms_filtered, 2)) + ' grms Filtered')

plt.tight_layout()
plt.show()