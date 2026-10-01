import FIDTraceFunctions.FIDTraceFunctions as f
import matplotlib.pyplot as plt
import numpy as np
<<<<<<< HEAD:test.py
=======
time, voltage = f.load_fid('src/acquisition_1090.txt')


popt, pcov, perr = f.fit_time_trace(time, voltage)

>>>>>>> 46215c4b0161a9dfbef3bfc39c5f8ea86a2437a5:src/test.py

time, voltage, meta = f.load_fid('acquisition_2267.txt', "BNCSynthesizer.frequency", "BNCSynthesizer.offset")
frequency, offset = meta['BNCSynthesizer.frequency'], meta['BNCSynthesizer.offset']
fig, ax = plt.subplots()

freqs, fft = f.fft(time, voltage, remove_if=True)


ax.plot(freqs / 1e6 + frequency - offset*1e-3, np.abs(fft))
#ax.plot(time, voltage)

plt.show()