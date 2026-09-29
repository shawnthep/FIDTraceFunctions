import functions as f
import matplotlib.pyplot as plt
time, voltage = f.load_fid('acquisition_2267.txt')


fig, ax = plt.subplots()

ax.plot(time, voltage)

plt.savefig('test_plot.png')