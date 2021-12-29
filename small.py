import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
from playsound import playsound as play
import IPython
from scipy.signal import spectrogram, lfilter, freqz, tf2zpk
from scipy.io import wavfile
seq = np.array([[1,2,3],[3,4,5]])
print(seq.T)
print(seq)
