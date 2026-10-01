import configparser
import numpy as np
from scipy.signal import get_window
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def load_fid(path, *fields):
    #fields here are represented as "gage.samplerate", "composer.chCwidth", etc.
    
    header_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            clean = line.strip()
            if "--- DATA ---" in clean:
                break
            if clean.startswith("#"):
                content = clean.lstrip("#").strip()
                if content:
                    header_lines.append(content)

    # Parse INI header
    config = configparser.ConfigParser(delimiters="=", allow_no_value=True)
    config.read_string("\n".join(header_lines))


    sample_rate = config["GaGe"].getint("samplerate")

    # Load raw voltage data and compute time vector
    voltage = np.loadtxt(path)
    time_axis = np.linspace(
        0, len(voltage) / sample_rate, len(voltage), endpoint=False
    )
    meta = {}
    for field in fields:
        section, _, key = field.partition(".") #uses gage.samplerate and converts it to "gage" and "samplerate"

        if config.has_option(section, key):

            meta[field] = _auto_type(config[section][key])
        else:

            warnings.warn(f"{path}: {field!r} not found in config")
    
    if fields:
    
        return time_axis, voltage, meta
    else:
        return time_axis, voltage

def _auto_type(value):

    for cast in (int, float, str):

        try:
            return cast(value)
        except ValueError:
            pass
    return value


def fft(time, voltage, window="box", remove_if=False, if_freq=47e6, fit_range=(500e-6, 900e-6), range = (20e6)):
    t = np.asarray(time, dtype=float)
    data = np.asarray(voltage, dtype=float)
    dt = np.median(np.diff(t))


    if remove_if:
        basis = np.column_stack([np.sin(2*np.pi*if_freq*t), np.cos(2*np.pi*if_freq*t)])
        i0, i1 = np.searchsorted(t, fit_range)
        coef, *_ = np.linalg.lstsq(basis[i0:i1], data[i0:i1], rcond=None)
        data = data - basis @ coef

    freqs = np.fft.rfftfreq(len(data), 5e-9)
    print(dt)
    #only returns the range of interest ie. 20MHz.

    mask = (freqs >= if_freq - range) & (freqs <= if_freq + range)

    #return frequencies in Hz
    return freqs[mask], np.fft.rfft(data * get_window(window, len(data)))[mask]
def fit_time_trace(time, voltage, start_time = 50e-6, end_time = 300e-6, return_fig = False):
    basis = np.column_stack([np.sin(2*np.pi*47e6*time), np.cos(2*np.pi*47e6*time)])
    i0, i1 = np.searchsorted(time, (500e-6, 1000e-6))
    coef, *_ = np.linalg.lstsq(basis[i0:i1], voltage[i0:i1], rcond=None)
    voltage = voltage - basis @ coef
    mask = (time > start_time) & (time < end_time)
    time_slice, voltage_slice = time[mask], voltage[mask]

    def func_to_fit(t, a, w1, T21, b, w2, T22, p1, p2, offset):

        return a*np.cos(2*np.pi*w1*t + p1)*np.exp(-t/T21) + b*np.cos(2*np.pi*w2*t + p2)*np.exp(-t/T22) + offset
    amplitude_guess = np.max(np.abs(voltage_slice[0:40]))
    initial_guess = [amplitude_guess, 47.335e6, 200e-6,amplitude_guess, 47.46e6, 200e-6, 0,0, 0]

    popt, pcov = curve_fit(func_to_fit, time_slice, voltage_slice, p0 = initial_guess)

    a, w1, T21, b, w2, T22, p1, p2, offset = popt

    print(f'Frequencies = {w1/1e6}MHz and {w2/1e6}Mhz\n')
    print(f'Difference in frequencys = {(w2-w1) / 1e3}kHz\n')
    #print(f'Amplitudes = {a, b}a.u')
    print(f'Dephasing => {T21 / 1e-6}us ({1 / (np.pi*T21) / 1e3}kHz) and {T22 / 1e-6}us ({1 / (np.pi*T22) / 1e3}kHz)')
    if return_fig:
        fig, ax = plt.subplots()
        ax.plot(time_slice / 1e-6, voltage_slice)

        ax.plot(time_slice / 1e-6, func_to_fit(time_slice, *popt))
        ax.set_ylabel('Amplitude (a.u.)')
        ax.set_xlabel('Time (us)')
        plt.show()
    return popt, pcov