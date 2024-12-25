import Lab3Functions as lf3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, filtfilt, hilbert
# Funktionen zur EMG-Verarbeitung
def remove_offset(emg):
    """Entfernt Mittel um Offset zu eliminieren"""
    return emg - np.mean(emg)

def bandpass_filter(emg, lowcut, highcut, fs, order=4):
    """
    Butterworth Filter für Frequenzen zwischen lowcut und highcut
    """
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, emg)

def rectify_signal(emg):
    """Gleichrichtung des Signals"""
    return np.abs(emg)

def calculate_envelope(emg, fs, cutoff=3, order=2):
    """Berechnung der Hüllkurve mit Hilbert-Transformation"""
    nyquist = 0.5 * fs
    low = cutoff / nyquist
    b, a = butter(order, low, btype='low')
    analytic_signal = hilbert(emg)
    envelope = np.abs(analytic_signal)
    return filtfilt(b, a, envelope)

def process_emg(data, fs):
    """
    Verarbeitet EMG-Daten: Offset, Filterung, Gleichrichtung, Hüllkurve
    """
    emg_no_offset = remove_offset(data['emg'])
    emg_filtered = bandpass_filter(emg_no_offset, 20, 190, fs)
    emg_rectified = rectify_signal(emg_filtered)
    emg_envelope = calculate_envelope(emg_filtered, fs)
    
    result = pd.DataFrame({
        't': data['t'],
        'emg_no_offset': emg_no_offset,
        'emg_filtered': emg_filtered,
        'emg_rectified': emg_rectified,
        'emg_envelope': emg_envelope
    })
    return result

# MVC-Wert (Maximum Voluntary Contraction)
MVC = 18.157144482993942  # mV
sampling_frequency = 411.76  # Hz
weights, mvc, fatigue = lf3.import_data(';')
mvc_processed = process_emg(mvc, sampling_frequency)
weights_processed = process_emg(weights, sampling_frequency)
fatigue_processed = process_emg(fatigue, sampling_frequency)


def process_emg(data, fs):
    """
    Verarbeitet EMG-Daten: Offset, Filterung, Gleichrichtung, Hüllkurve
    """
    emg_no_offset = remove_offset(data['emg'])
    emg_filtered = bandpass_filter(emg_no_offset, 20, 190, fs)
    emg_rectified = rectify_signal(emg_filtered)
    emg_envelope = calculate_envelope(emg_filtered, fs)
    
    result = pd.DataFrame({
        't': data['t'],
        'emg_no_offset': emg_no_offset,
        'emg_filtered': emg_filtered,
        'emg_rectified': emg_rectified,
        'emg_envelope': emg_envelope
    })
    return result

#wheights plotting
plt.figure()

plt.subplot(2, 1, 1)

plt.plot(weights_processed['t'], weights_processed['emg_filtered'], label='Filtered EMG')

plt.legend()
plt.title('Weights Experiment')
plt.xlabel('Time [s]')
plt.ylabel('EMG [mV]')
plt.grid()
plt.show()
