import configparser
import numpy as np
from scipy.signal import get_window
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


def fft(time, voltage, window="hann", remove_if=False, if_freq=47e6,
        fit_range=(100e-6, 900e-6), range = (20e6)):
    t = np.asarray(time, dtype=float)
    data = np.asarray(voltage, dtype=float)
    dt = np.median(np.diff(t))

    if remove_if:
        basis = np.column_stack([np.sin(2*np.pi*if_freq*t), np.cos(2*np.pi*if_freq*t)])
        i0, i1 = np.searchsorted(t, fit_range)
        coef, *_ = np.linalg.lstsq(basis[i0:i1], data[i0:i1], rcond=None)
        data = data - basis @ coef

    freqs = np.fft.rfftfreq(len(data), dt)
    mask = (freqs > if_freq - range) & (freqs < if_freq + range)
    return freqs[mask], np.fft.rfft(data * get_window(window, len(data)))[mask]