# Syllabic rate squelch/voice activity detection (VAD)

## Introduction

Based on the existence of a number of (Motorola, Harris, Rockwell, Codan)
hardware-based implementations of squelch circuits that depend on the
characteristics of human speech to detect a voice signal and open a squelch,
rather than the much older process of simply detecting RF carrier energy
and using it to open a squelch, it was considered that DSP techniques should
be able to be used to accomplish the same thing in software.

After many years of just thinking about it, and the recent fame of the
ChatGPT AI, it was decided to see how far it could get with the problem:

> Q: Write DSP code to detect voice-band energy (500-3,000 Hz) which is
> varying at a rate of 0.4 to 3.5 Hz

> Sure, here’s some example DSP code in Python that could be used to
> detect voice-band energy in the frequency range of 500-3,000 Hz and
> determine the rate at which the energy is varying:

```
import numpy as np
from scipy.signal import butter, lfilter

# Define the sampling rate and frequency range of interest
fs = 8000 # Hz
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

# Loop over incoming audio samples
for i in range(num_samples):

    # Read the next audio sample
    x = read_sample()

    # Apply the bandpass filter
    y = lfilter(b, a, x)

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
    else:
        # No speech activity detected, set rate to zero
        rate = 0
```

> This code applies a fourth-order Butterworth bandpass filter to
> isolate the voice-band energy in the frequency range of 500-3,000
> Hz. The energy of the filtered signal is then computed and compared
> to the previous energy value to determine the rate of change. This
> rate of change is then smoothed using a low-pass filter with a
> cutoff frequency of 3.5 Hz.

> If the smoothed rate of change exceeds a threshold of 0.4 Hz, speech
> activity is detected and the rate value is set to the smoothed
> delta_energy value. If no speech activity is detected, the rate
> value is set to zero.

> Note that this is just one example of how to detect voice-band
> energy and determine the rate of change, and there are many
> variations and improvements that could be made depending on the
> specific application and requirements.

> Q: How would the code above be altered to sample audio continuously
> from an audio interface?

> To sample audio continuously from an audio interface, you would need
> to modify the code to read samples from the interface in a loop
> rather than assuming a fixed number of samples.

> Here’s an updated version of the code that continuously reads audio
> samples from an input stream:

```
import numpy as np
from scipy.signal import butter, lfilter
import sounddevice as sd

# Define the sampling rate and frequency range of interest
fs = 8000 # Hz
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
    else:
        # No speech activity detected, set rate to zero
        rate = 0

# Start the input stream with the specified callback function
with sd.InputStream(callback=audio_callback, blocksize=1024, samplerate=fs):
    # Run indefinitely
    while True:
        pass
```

> This code uses the sounddevice package to create an input stream
> with a specified sample rate of 8000 Hz and a block size of 1024
> samples. The audio_callback function is called for each block of
> samples, and applies the bandpass filter to the input samples,
> computes the energy and rate of change, and updates the global
> variables.

> The while True loop at the end of the code runs indefinitely and
> allows the input stream to continue reading samples until the
> program is interrupted or terminated.

> Q: What limits the voice-band energy detection sensitivity in this
> algorithm?

> The voice-band energy detection sensitivity in this algorithm is
> limited by several factors:

> 1. Sampling rate: The frequency resolution of the algorithm is limited
> by the sampling rate. According to the Nyquist-Shannon sampling
> theorem, the highest frequency that can be accurately sampled is
> half the sampling rate. Therefore, the sensitivity to voice-band
> energy at higher frequencies may be limited by the sampling rate.

> 2. Bandpass filter design: The sensitivity to voice-band energy within
> the 500-3000 Hz range is largely determined by the design of the
> bandpass filter. In this algorithm, a fourth-order Butterworth
> filter is used, which provides a flat response within the passband
> and steep attenuation outside the passband. However, other filter
> designs may provide different sensitivities to voice-band energy.

> 3. Noise level: The sensitivity to voice-band energy is also affected
> by the level of background noise. If the noise level is high
> relative to the speech signal, the algorithm may have difficulty
> detecting the speech energy and distinguishing it from the noise.

> 4. Threshold level: The sensitivity of the algorithm can also be
> controlled by the threshold level set for the rate of change of
> energy. In this algorithm, a threshold of 0.4 Hz is used to detect
> speech activity. If the threshold level is set too low, the
> algorithm may generate false positives, while if the threshold level
> is set too high, it may fail to detect speech activity.

> 5. Time window: The sensitivity of the algorithm may also depend on the
> duration of the time window used to compute the energy and rate of
> change. In this algorithm, a single block of samples is used for
> each calculation, which corresponds to a time window of about 128
> ms. Using a longer time window may provide a more accurate estimate
> of the energy and rate of change, but may also increase the delay in
> detecting speech activity.

Surprising, some casual testing of this code (see test.py) using a microphone
for input seems to show that it does function in a rudimentary way.

## References

1. ["Smart" Squelch for SSB](https://www.repeater-builder.com/projects/jpl-vox-sq/smart-squelch-for-ssb.pdf)
2. [The Jet Propulsion Laboratory's "Smart VOX"](https://www.repeater-builder.com/projects/jpl-vox-sq/ssb-squelch.html)
3. [A VAD/VOX Algorithm for Amateur Radio Applications](https://www.thinkmind.org/articles/tele_v7_n12_2014_4.pdf)
4. [A Hybrid VOX System Using Emulated Hardware Behaviors](http://www.thinkmind.org/articles/icdt_2012_5_50_20080.pdf)
5. [Audio Frequency Squelch System](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/3873926)
6. [Audio Frequency Squelch System](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/3873925)
7. [Squelch Circuit](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/4414689)
8. [Multiband Squelch Method and Apparatus](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6397050)







