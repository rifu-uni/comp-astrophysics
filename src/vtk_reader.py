import numpy as np
import os

def read_pluto_vtk(filepath):
    with open(filepath, 'rb') as f:
        raw = f.read()

    lines = raw.decode('latin-1').split('\n')
    for line in lines:
        if line.startswith('X_COORDINATES'):
            x = np.array([float(v) for v in line.split()[3:]])
            break

    fields = {}
    pos = raw.find(b'SCALARS')
    while pos != -1:
        eol = raw.find(b'\n', pos)
        parts = raw[pos:eol].split()
        name = parts[1].decode()
        lut_eol = raw.find(b'\n', eol + 1)
        data_start = lut_eol + 1
        n_cells = len(x) - 1 if len(x) > 1 else 128
        field_raw = raw[data_start:data_start + n_cells * 4]
        fields[name] = np.frombuffer(field_raw, dtype='>f4').copy()
        pos = raw.find(b'SCALARS', data_start + n_cells * 4 + 1)

    x_active = np.linspace(1/(2*len(x)-2), 1 - 1/(2*len(x)-2), len(x)-1)

    return x_active, fields

def read_pluto_vtk_coarse(filepath):
    with open(filepath, 'rb') as f:
        raw = f.read()

    fields = {}
    pos = raw.find(b'SCALARS')
    while pos != -1:
        eol = raw.find(b'\n', pos)
        parts = raw[pos:eol].split()
        name = parts[1].decode()
        lut_eol = raw.find(b'\n', eol + 1)
        data_start = lut_eol + 1
        field_raw = raw[data_start:data_start + 128 * 4]
        fields[name] = np.frombuffer(field_raw, dtype='>f4').copy()
        pos = raw.find(b'SCALARS', data_start + 128 * 4 + 1)

    return fields
