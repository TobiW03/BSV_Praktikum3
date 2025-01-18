import scipy
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


def calculate_median_frequency(power, frequencies):
    area_freq = scipy.integrate.cumtrapz(power, frequencies, initial=0)
    total_power = area_freq[-1]
    median_freq = frequencies[np.where(area_freq >= total_power / 2)[0][0]]
    return median_freq


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

# Define beginning, middle, and end segments with fixed 0.5 second intervals
def get_fixed_segment(data, start_time, duration=0.5):
    end_time = start_time + duration
    return data[(data['t'] >= start_time) & (data['t'] < end_time)]

# Anfang, Mitte, Ende von Bursts
fatigue1_processed_beginning = get_fixed_segment(fatigue1_processed, 1) # 0-0.5s
fatigue2_processed_beginning = get_fixed_segment(fatigue2_processed, 1)
fatigue3_processed_beginning = get_fixed_segment(fatigue3_processed, 1)

fatigue1_processed_middle = get_fixed_segment(fatigue1_processed, 7)   # 7-7.5s
fatigue2_processed_middle = get_fixed_segment(fatigue2_processed, 7)
fatigue3_processed_middle = get_fixed_segment(fatigue3_processed, 7)

fatigue1_processed_end = get_fixed_segment(fatigue1_processed, 13)    # 13-13.5s
fatigue2_processed_end = get_fixed_segment(fatigue2_processed, 13)
fatigue3_processed_end = get_fixed_segment(fatigue3_processed, 13)

print(fatigue1_processed_beginning['emg_filtered'])
power1_beginning,frequencies1_beginning= lf3.get_power(fatigue1_processed_beginning['emg_filtered'], 411.76)
power2_beginning,frequencies2_beginning= lf3.get_power(fatigue2_processed_beginning['emg_filtered'], 411.76)
power3_beginning,frequencies3_beginning= lf3.get_power(fatigue3_processed_beginning['emg_filtered'], 411.76)

power1_middle,frequencies1_middle= lf3.get_power(fatigue1_processed_middle['emg_filtered'], 411.76)
power2_middle,frequencies2_middle= lf3.get_power(fatigue2_processed_middle['emg_filtered'], 411.76)
power3_middle,frequencies3_middle= lf3.get_power(fatigue3_processed_middle['emg_filtered'], 411.76)

power1_end,frequencies1_end= lf3.get_power(fatigue1_processed_end['emg_filtered'], 411.76)
power2_end,frequencies2_end= lf3.get_power(fatigue2_processed_end['emg_filtered'], 411.76)
power3_end,frequencies3_end= lf3.get_power(fatigue3_processed_end['emg_filtered'], 411.76)

#### Aufgabe 9 Plot ####
# Plot frequency1 for Beginning segment only
plt.figure()

# Raw power spectrum for Beginning segment
plt.plot(frequencies1_beginning, power1_beginning, label='Rohes Leistungsspektrum', linestyle='-')

# Filtered power spectrum using Butterworth filter with cutoff frequency of 40Hz
def butter_lowpass_filter(data, cutoff, fs, order=2):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)

filtered_power1_beginning = butter_lowpass_filter(power1_beginning, 40, sampling_frequency)

plt.plot(frequencies1_beginning, filtered_power1_beginning, label='Gefiltertes Leistungsspektrum', linestyle='--')

# Average frequency of filtered power spectrum for Beginning segment
average_freq = np.sum(frequencies1_beginning * filtered_power1_beginning) / np.sum(filtered_power1_beginning)

plt.axvline(average_freq, color='r', linestyle='-.', label='Durchschnittsfrequenz')

plt.xlabel('Frequenz in Hz')
plt.ylabel('Leistung in a.u.')
plt.legend()
plt.show()

###Aufgabe 9 Plot####

# Subplot für Anfang, Mitte, Ende für die drei Segmente
fig, axs = plt.subplots(3, 3, figsize=(15, 10))

segments = ['Anfang', 'Mitte', 'Ende']
fatigue_data = [
    (fatigue1_processed_beginning, fatigue1_processed_middle, fatigue1_processed_end),
    (fatigue2_processed_beginning, fatigue2_processed_middle, fatigue2_processed_end),
    (fatigue3_processed_beginning, fatigue3_processed_middle, fatigue3_processed_end)
]

for i, (beginning, middle, end) in enumerate(fatigue_data):
    for j, segment in enumerate([beginning, middle, end]):
        power, frequencies = lf3.get_power(segment['emg_filtered'], sampling_frequency)
        filtered_power = butter_lowpass_filter(power, 40, sampling_frequency)
        
        # Calculate average frequency
        average_freq = np.sum(frequencies * filtered_power) / np.sum(filtered_power)

        axs[i, j].plot(frequencies, power, label='Rohes Leistungsspektrum', linestyle='-')
        axs[i, j].plot(frequencies, filtered_power, label='Gefiltertes Leistungsspektrum', linestyle='--')
        axs[i, j].axvline(average_freq, color='r', linestyle='-.', label='Durchschnittsfrequenz')
        axs[i, j].set_title(f'Segment {segments[j]} - Fatigue {i+1}')
        axs[i, j].set_xlabel('Frequenz in Hz')
        axs[i, j].set_ylabel('Leistung in a.u.')
        axs[i, j].legend()

plt.tight_layout()
plt.show()


#### Aufgabe 10 Plot ####
# Plot frequency1 for Beginning segment only
plt.figure()

filtered_power1_beginning = butter_lowpass_filter(power1_beginning, 40, sampling_frequency)

plt.plot(frequencies1_beginning, filtered_power1_beginning, label='Gefiltertes Leistungsspektrum', linestyle='-')

# Marking the fiber type regions with hatching patterns for better visibility in black and white prints

plt.xlabel('Frequenz in Hz')
plt.ylabel('Leistung in a.u.')
plt.show()

###Aufgabe 10 Plot####


####Aufgabe 11 Plot###
# Define beginning, middle, and end segments with fixed 0.5 second interval
# Analyze fatigue data and calculate median frequencies
median_freq1_beginning = calculate_median_frequency(power1_beginning, frequencies1_beginning)
median_freq1_middle = calculate_median_frequency(power1_middle, frequencies1_middle)
median_freq1_end = calculate_median_frequency(power1_end, frequencies1_end)

median_freeq2_beginning = calculate_median_frequency(power2_beginning, frequencies2_beginning)
median_freeq2_middle = calculate_median_frequency(power2_middle, frequencies2_middle)
median_freeq2_end = calculate_median_frequency(power2_end, frequencies2_end)

median_freeq3_beginning = calculate_median_frequency(power3_beginning, frequencies3_beginning)
median_freeq3_middle = calculate_median_frequency(power3_middle, frequencies3_middle)
median_freeq3_end = calculate_median_frequency(power3_end, frequencies3_end)

# Combine median frequencies for all fatigues into a list
median_freqs = [median_freq1_beginning, median_freq1_middle, median_freq1_end]
median_freqs2 = [median_freeq2_beginning, median_freeq2_middle, median_freeq2_end]
median_freqs3 = [median_freeq3_beginning, median_freeq3_middle, median_freeq3_end]

# Plot median frequencies for all fatigues with different markers
plt.figure()
plt.plot([1/14*100, 50, 13/14*100], median_freqs, label='Ermüdung 1', marker='o')
plt.plot([1/14*100, 50, 13/14*100], median_freqs2, label='Ermüdung 2', marker='s')
plt.plot([1/14*100, 50, 13/14*100], median_freqs3, label='Ermüdung 3', marker='^')
plt.xlabel('Zeitpunkt der Messung in %')
plt.ylabel('Medianfrequenz in Hz')
plt.xlim(0, 100)  # Set y-axis limit to 100
plt.legend()
plt.show()
####Aufgabe 11 Plot###





