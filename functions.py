import configparser
import numpy as np
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

