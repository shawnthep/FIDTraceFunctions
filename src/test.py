import functions as f
import matplotlib.pyplot as plt
import numpy as np
time, voltage = f.load_fid('acquisition_2267.txt')


fig, ax = plt.subplots()

freqs, fft = f.fft(time, voltage)

ax.plot(freqs, np.abs(fft))

plt.savefig('test_fft.png')
