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

sampling_frequency = 411.76 # Hz
weights, mvc, fatigue = lf3.import_data(';')
mvc_processed = process_emg(mvc, sampling_frequency)
weights_processed = process_emg(weights, sampling_frequency)
fatigue_processed = process_emg(fatigue, sampling_frequency)


# Plot fatigue
'''plt.figure()
plt.plot(fatigue_processed['t'], fatigue_processed['emg_filtered'], label='EMG Hüllkurve')
plt.xlabel('Zeit [s]')
plt.ylabel('EMG Hüllkurve [mV]')
plt.show()
print(len(fatigue_processed['emg_filtered']))'''


# Fatiguedaten einlesen
fatigue1 = pd.read_csv('Fatigue1.txt', delimiter=';')  # 1
fatigue2 = pd.read_csv('Fatigue2.txt', delimiter=';')  # 2
fatigue3 = pd.read_csv('Fatigue2.txt', delimiter=';')  # 3

# Indexlängen der Fatiguedaten finden
fatigue1_length = len(fatigue1)
fatigue2_length = len(fatigue2)
fatigue3_length = len(fatigue3)

# Fatiguedaten verarbeiten
fatigue1_processed = process_emg(fatigue.iloc[:fatigue1_length], sampling_frequency)
fatigue2_processed = process_emg(fatigue.iloc[fatigue1_length:fatigue1_length + fatigue2_length], sampling_frequency)
fatigue3_processed = process_emg(fatigue.iloc[fatigue1_length + fatigue2_length:], sampling_frequency)


# Zeit von Millisekunden auf Sekunden umwandeln
fatigue1_processed['t'] = (fatigue1_processed['t'] - fatigue1_processed['t'].iloc[0]) / 1000
fatigue2_processed['t'] = (fatigue2_processed['t'] - fatigue2_processed['t'].iloc[0]) / 1000
fatigue3_processed['t'] = (fatigue3_processed['t'] - fatigue3_processed['t'].iloc[0]) / 1000

# Daten auf die Haltephase beschränken
fatigue1_processed = fatigue1_processed[(fatigue1_processed['t'] >= 3) & (fatigue1_processed['t'] <= 17)]
fatigue2_processed = fatigue2_processed[(fatigue2_processed['t'] >= 4) & (fatigue2_processed['t'] <= 18)]
fatigue3_processed = fatigue3_processed[(fatigue3_processed['t'] >= 1) & (fatigue3_processed['t'] <= 15)]

# Zeit von fatigue1-3 jeweils von 0 starten lassen
fatigue1_processed['t'] = fatigue1_processed['t'] - fatigue1_processed['t'].iloc[0]
fatigue2_processed['t'] = fatigue2_processed['t'] - fatigue2_processed['t'].iloc[0]
fatigue3_processed['t'] = fatigue3_processed['t'] - fatigue3_processed['t'].iloc[0]



