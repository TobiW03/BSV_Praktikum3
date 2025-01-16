# Imports
import Lab3Functions as lf3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import butter, filtfilt, hilbert
from plotly.subplots import make_subplots
import plotly.graph_objects as go

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

# Sampling-Frequenz
sampling_frequency = 411.76  # Hz

# Daten laden und verarbeiten
weights, mvc, fatigue = lf3.import_data(';')
mvc_processed = process_emg(mvc, sampling_frequency)
weights_processed = process_emg(weights, sampling_frequency)
fatigue_processed = process_emg(fatigue, sampling_frequency)
"""
# Interaktive Auswahl der Bursts
Smvc_s, mvc_e, weights_s, weights_e, fatigue_s, fatigue_e = lf3.get_bursts(
    mvc_processed['emg_filtered'],
    weights_processed['emg_filtered'],
    fatigue_processed['emg_filtered']
)

# Berechnung des persönlichen MVC
mvc_burst_means = []
for start, end in zip(mvc_s, mvc_e):
    burst_segment = mvc_processed['emg_envelope'][start:end]
    burst_mean = np.mean(burst_segment)
    mvc_burst_means.append(burst_mean)

# Mitteln der drei Bursts
personal_mvc = np.mean(mvc_burst_means)

# Ausgabe
print("Startpunkte der Bursts:", mvc_s)
print("Endpunkte der Bursts:", mvc_e)
print("Mittelwerte der Bursts:", mvc_burst_means)
print("Persönlicher MVC:", personal_mvc,"mV")
"""
# Visualisierung der Verarbeitungsstufen
fig = make_subplots(
    rows=2, cols=2
)

fig.add_trace(go.Scatter(x=(mvc['t'] / 1000), y=mvc_processed['emg_no_offset'], mode='lines'), row=1, col=1)
fig.add_trace(go.Scatter(x=(mvc['t'] / 1000), y=mvc_processed['emg_filtered'], mode='lines'), row=1, col=2)
fig.add_trace(go.Scatter(x=(mvc['t'] / 1000), y=mvc_processed['emg_rectified'], mode='lines'), row=2, col=1)
fig.add_trace(go.Scatter(x=(mvc['t'] / 1000), y=mvc_processed['emg_envelope'], mode='lines'), row=2, col=2)

fig.update_layout(
    height=1200,
    width=1000,
    font=dict(size=14),
    title_font=dict(size=18),
    showlegend=False,
    template="plotly_white",
    xaxis_title="Zeit in s",
    yaxis_title="EMG-Signal in mV"
)

fig.update_xaxes(showgrid=False)
fig.update_yaxes(showgrid=False)
fig.update_yaxes(title="EMG-Signal in mV", row=1, col=1)
fig.update_yaxes(title="EMG-Signal in mV", row=1, col=2)
fig.update_yaxes(title="EMG-Signal in mV", row=2, col=1)
fig.update_yaxes(title="EMG-Signal in mV", row=2, col=2)
fig.update_xaxes(title="Zeit in s", row=1, col=1)
fig.update_xaxes(title="Zeit in s", row=1, col=2)
fig.update_xaxes(title="Zeit in s", row=2, col=1)
fig.update_xaxes(title="Zeit in s", row=2, col=2)

fig.show()