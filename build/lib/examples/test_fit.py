import FIDTraceFunctions.FIDTraceFunctions as f
import matplotlib.pyplot as plt
import numpy as np

time, voltage, meta = f.load_fid('src/acquisition_1099.txt', "BNCSynthesizer.frequency", "BNCSynthesizer.offset")
frequency, offset = meta['BNCSynthesizer.frequency'], meta['BNCSynthesizer.offset']



popt, pcov = f.fit_time_trace(time, voltage, return_fig=True)
