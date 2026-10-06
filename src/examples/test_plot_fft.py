import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import FIDTraceFunctions.FIDTraceFunctions as f 
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(str(Path(__file__).parent))

time, voltage, meta = f.load_fid('src/examples/acquisition_1090.txt', "BNCSynthesizer.frequency", "BNCSynthesizer.offset")
frequency, offset = meta['BNCSynthesizer.frequency'], meta['BNCSynthesizer.offset']

mask = (time > 20e-6) 

time, voltage = time[mask], voltage[mask]
def func_to_fit(t, a, w1, T21, b, w2, T22, p1, p2, offset):

        return a*np.cos(2*np.pi*w1*t + p1)*np.exp(-t/T21) + b*np.cos(2*np.pi*w2*t + p2)*np.exp(-t/T22) + offset
'''fig, ax = plt.subplots()

freqs, fft = f.fft(time, voltage, remove_if=True)

#ax.plot(freqs / 1e6 , np.abs(fft))
ax.plot(time, voltage)
ax.plot(time, func_to_fit(time, 0.003, 47.335e6, 200e-6,0.003, 47.46e6, 200e-6, 0,0, 0))

plt.show()'''

