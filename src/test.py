import FIDTraceFunctions.FIDTraceFunctions as f
import matplotlib.pyplot as plt
import numpy as np
time, voltage = f.load_fid('src/acquisition_1090.txt')


popt, pcov, perr = f.fit_time_trace(time, voltage)


fig, ax = plt.subplots()

freqs, fft = f.fft(time, voltage)

ax.plot(freqs, np.abs(fft))

plt.savefig('test_fft.png')
