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


# MVC-Wert (Maximum Voluntary Contraction)
MVC = 18.157144482993942  # mV
sampling_frequency = 411.76  # Hz
weights, mvc, fatigue = lf3.import_data(';')
mvc_processed = process_emg(mvc, sampling_frequency)
weights_processed = process_emg(weights, sampling_frequency)
fatigue_processed = process_emg(fatigue, sampling_frequency)

# Normalize the EMG envelope to MVC
weights_processed['emg_envelope_normalized'] = weights_processed['emg_envelope'] / MVC * 100

# Weights plotting
"""plt.figure()

plt.subplot(2, 1, 1)
plt.plot(weights_processed['t'], weights_processed['emg_envelope_normalized'], label='Filtered EMG (Normalized)')
plt.legend()
plt.title('Weights Experiment')
plt.xlabel('Time [s]')
plt.ylabel('EMG [% MVC]')
plt.grid()
plt.show()"""

import Lab3Functions as lf3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# MVC-Wert (Maximum Voluntary Contraction)
MVC = 18.157144482993942  # mV
sampling_frequency = 411.76  # Hz

# Gewichtsdaten einlesen
weights1 = pd.read_csv('Weight1.txt', delimiter=';')  # 2.5 kg
weights2 = pd.read_csv('Weight2.txt', delimiter=';')  # 5 kg
weights3 = pd.read_csv('Weight3.txt', delimiter=';')  # 10 kg

# Indexlängen der Gewichtsdaten finden
weights1_length = len(weights1)
weights2_length = len(weights2)
weights3_length = len(weights3)
weights1_processed = weights_processed[:weights1_length]
weights2_processed = weights_processed[weights1_length:weights1_length+weights2_length]
weights3_processed = weights_processed[weights1_length+weights2_length:weights1_length+weights2_length+weights3_length]

# Zeit von weights1-3 jeweils von 0 starten lassen
weights1_processed['t'] = (weights1_processed['t'] - weights1_processed['t'].iloc[0]) / 1000
weights2_processed['t'] = (weights2_processed['t'] - weights2_processed['t'].iloc[0]) / 1000
weights3_processed['t'] = (weights3_processed['t'] - weights3_processed['t'].iloc[0]) / 1000

# Nur die ersten 10 Sekunden der Daten behalten
weights1_processed = weights1_processed[weights1_processed['t'] <= 10]S
weights2_processed = weights2_processed[weights2_processed['t'] <= 10]
weights3_processed = weights3_processed[weights3_processed['t'] <= 10]

# Plotten von weights1
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(weights1_processed['t'], weights1_processed['emg_envelope_normalized'], label='2,5 kg')
plt.title('Experiment - 2,5 kg')
plt.xlabel('Zeit in s')
plt.ylabel('Anteil von MVC in %')
plt.grid()

# Plotten von weights2
plt.subplot(3, 1, 2)
plt.plot(weights2_processed['t'], weights2_processed['emg_envelope_normalized'], label='5 kg')
plt.title('Experiment - 5 kg')
plt.xlabel('Zeit in s')
plt.ylabel('Anteil von MVC in %')
plt.grid()

# Plotten von weights3
plt.subplot(3, 1, 3)
plt.plot(weights3_processed['t'], weights3_processed['emg_envelope_normalized'], label='10 kg')
plt.title('Experiment - 10 kg')
plt.xlabel('Zeit in s')
plt.ylabel('Anteil von MVC in %')
plt.grid()

plt.tight_layout()
plt.show()
