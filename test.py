import functions as f
import matplotlib.pyplot as plt
import numpy as np

time, voltage, meta = f.load_fid('acquisition_2267.txt', "BNCSynthesizer.frequency", "BNCSynthesizer.offset")
frequency, offset = meta['BNCSynthesizer.frequency'], meta['BNCSynthesizer.offset']
fig, ax = plt.subplots()

freqs, fft = f.fft(time, voltage, remove_if=True)


ax.plot(freqs / 1e6 + frequency - offset*1e-3, np.abs(fft))
#ax.plot(time, voltage)

plt.show()