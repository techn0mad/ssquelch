import numpy as np
from scipy.signal import butter, lfilter
import sounddevice as sd

# Define the sampling rate and frequency range of interest
fs = 44100 # Hz
f_low = 500 # Hz
f_high = 3000 # Hz

# Design a bandpass filter to isolate the voice band energy
b, a = butter(4, [f_low / (fs/2), f_high / (fs/2)], btype='bandpass')

# Initialize variables for tracking the energy and rate of change
energy = 0
prev_energy = 0
delta_energy = 0
delta_energy_smooth = 0
rate = 0

# Define the callback function to read audio samples from the input stream
def audio_callback(indata, frames, time, status):
    global energy, prev_energy, delta_energy, delta_energy_smooth, rate

    # Apply the bandpass filter to the input samples
    y = lfilter(b, a, indata[:, 0])

    # Compute the energy of the filtered signal
    energy = np.sum(np.square(y))

    # Compute the rate of change of the energy using a first-order difference
    delta_energy = energy - prev_energy

    # Smooth the rate of change using a low-pass filter with a cutoff frequency of 3.5 Hz
    alpha = 1 - np.exp(-2*np.pi*3.5/fs)
    delta_energy_smooth = alpha*delta_energy_smooth + (1-alpha)*delta_energy

    # Update the previous energy value
    prev_energy = energy

    # Check if the rate of change exceeds the threshold of 0.4 Hz
    if delta_energy_smooth > 0.4:
        # Speech activity detected, set rate to the smoothed delta_energy value
        rate = delta_energy_smooth
        print("rate: ", rate, " energy: ", energy, " delta_energy: ", delta_energy)
    else:
        # No speech activity detected, set rate to zero
        rate = 0

# Start the input stream with the specified callback function
with sd.InputStream(callback=audio_callback, blocksize=8192, samplerate=fs):
    # Run indefinitely
    while True:
        pass
