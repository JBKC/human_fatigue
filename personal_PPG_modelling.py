'''
Determining fatigue levels from pre-recorded PPG sensor inputs
'''

import serial
import numpy as np
import matplotlib, matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, medfilt, spectrogram
import pywt
import pandas as pd

def extract_data(file):
    with open(file, "r") as f:
        lines = f.readlines()
        total_time = int(lines[1].strip())
        signal = lines[2:-1]
        signal = [int(line.strip()) for line in signal]  # Convert lines from string to float
        signal = np.array(signal)
        fs = len(signal) / (total_time/1000)
        print(fs)
    return signal, fs

def normalise(signal):
    mean = np.mean(signal, axis=0)
    normalised = signal - mean
    return normalised

def multilevel_decomposition(signal):

    '''
    Visualise the DWT multilevel decomposition
    '''

    wavelet = 'db2'
    max_levels = 12

    plt.figure(figsize=(10, 6))
    plt.subplot(max_levels + 1, 1, 1)
    plt.plot(signal)
    plt.title("Original Signal")
    plt.ylabel("Amplitude")
    for level in range(1, max_levels + 1):
        coeffs = pywt.wavedec(signal, wavelet, mode="per", level=level)
        coeff = coeffs[0]
        plt.subplot(max_levels + 1, 1, level + 1)
        # for i, coeff in enumerate(coeffs):
        plt.plot(coeff, label=f"Level {level} - Coeff {0}")
        plt.title(f"Decomposition at Level {level}")
        plt.ylabel("Coeffs")
    plt.tight_layout()
    plt.show()

def reconstruct_signal(signal):

    wavelet = 'db2'

    def high_noise(signal):
        '''
        Which high-frequency noise to remove (detail space)
        '''

        level = 2
        coeffs = pywt.wavedec(signal, wavelet, mode='per', level=level)

        for level in range(1, level +1):
            coeffs[level][:] = 0

        high_denoised_signal = pywt.waverec(coeffs, wavelet, mode='per')

        return high_denoised_signal

    def low_noise(signal):
        '''
        Which low-frequency noise to remove (approximation space)
        '''

        level = 8
        coeffs = pywt.wavedec(signal, wavelet, mode='per', level=level)
        coeffs[0][:] = 0
        print(coeffs)
        low_denoised_signal = pywt.waverec(coeffs, wavelet, mode='per')

        return low_denoised_signal

    denoised_signal = high_noise(signal)
    denoised_signal = low_noise(denoised_signal)

    plt.figure(figsize=(10, 4))
    plt.plot(denoised_signal, 'red')
    plt.title('Denoised Signal')

    plt.plot(ppg_signal, 'black')
    plt.title('Original Signal')

    plt.tight_layout()
    plt.show()

    return denoised_signal

def gabor_transform(signal, fs):
    window_length = 1                 # in seconds
    nperseg = int(window_length*fs)
    noverlap = nperseg // 4             # overlap between consecutive windows

    frequencies, times, Sxx = spectrogram(signal, fs=fs, nperseg=nperseg, noverlap=noverlap, scaling='spectrum')
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='auto', cmap='jet', vmin=-60, vmax=0)
    print(frequencies)
    plt.colorbar(label='Power (dB)')
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.xlim(0, max(times))  # Set x-axis limits from 0 to the maximum time value
    plt.ylim(0, 40)  # Set y-axis limits from 0 to the maximum frequency value
    plt.show()

file = '/Users/jamborghini/Documents/ARDUINO/JB DATA (top of wrist, 115600baud)/training data/post-training/26:10:23, late evening, post-training.txt'
ppg_signal, fs = extract_data(file)
ppg_signal = normalise(ppg_signal)

# multilevel_decomposition(ppg_signal)
reconstruct_signal(ppg_signal)
#gabor_transform(ppg_signal, fs)