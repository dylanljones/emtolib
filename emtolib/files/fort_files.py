# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2024-02-01

from pathlib import Path

import numpy as np
from scipy import constants as cc

ry2ev = cc.value("Rydberg constant times hc in eV")  # 13.605693122994


def read_sigma_z(path, unit="ry"):
    path = Path(path)
    unit = unit.lower()
    if unit not in ("ry", "ev"):
        raise ValueError(f"Invalid unit {unit}!")
    if not path.exists():
        raise FileNotFoundError(f"File {path.name} does not exist!")
    data = np.loadtxt(path)

    z = data[:, 0]
    sig_real = data[:, 1::2]
    sig_imag = data[:, 2::2]
    sig = sig_real + 1j * sig_imag
    nskip = 0
    for i, _z in enumerate(z):
        if _z != 0:
            nskip = i
            break
    # Remove first rows
    z, sig = z[nskip:], sig[nskip:]
    # Unpack spin channels
    _, idx = np.unique(z, return_index=True)
    n = len(z) // len(idx)
    z = np.split(z, n)[0]
    sig = np.array(np.split(sig, n))
    # Reshape to (n_channels, n_orbitals, n_energy)
    sig = np.swapaxes(sig, 1, 2)
    if unit == "ev":
        z = z * ry2ev
        sig = sig * ry2ev
    return z, sig


def read_sigma_iw(path, unit="ry"):
    path = Path(path)
    unit = unit.lower()
    if unit not in ("ry", "ev"):
        raise ValueError(f"Invalid unit {unit}!")
    if not path.exists():
        raise FileNotFoundError(f"File {path.name} does not exist!")
    data = np.loadtxt(path)

    iw = data[:, 0]
    sig_real = data[:, 1::2]
    sig_imag = data[:, 2::2]
    sig = sig_real + 1j * sig_imag
    # Unpack spin channels
    _, idx = np.unique(iw, return_index=True)
    n = len(iw) // len(idx)
    iw = np.split(iw, n)[0]
    sig = np.array(np.split(sig, n))
    # Reshape to (n_channels, n_orbitals, n_energy)
    sig = np.swapaxes(sig, 1, 2)
    # Channels are (at1_up, at1_dn, at2_up, at2_dn, ...)
    # Reshape to (natom, nspin, n_orbitals, n_energy)
    shape = sig.shape
    nchannels = shape[0]
    natom = nchannels // 2
    nspin = 2
    sig = sig.reshape((natom, nspin, *shape[1:]))

    if unit == "ev":
        sig = sig * ry2ev
        iw = iw * ry2ev
    return iw, sig
