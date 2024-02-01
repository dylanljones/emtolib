# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2024-02-01

from pathlib import Path

import numpy as np


def read_fermi_surface_file(path):
    """Load the data from a fermi-surface output file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File {path.name} does not exist!")

    fs_up = dict()
    fs_dn = dict()
    for line in path.read_text().splitlines():
        if not line:
            continue
        x, y, z = line.split()
        x, y, z = float(x), float(y), float(z)
        k = (x, y)
        if k in fs_up:
            fs_dn[k] = z
        else:
            fs_up[k] = z

    kpoints = np.array([list(p) for p in fs_up.keys()])
    kx = np.unique(kpoints[:, 0])
    ky = np.unique(kpoints[:, 1])
    fs = np.zeros((2, len(kx), len(ky)))
    for i, x in enumerate(kx):
        for j, y in enumerate(ky):
            k = (x, y)
            fs[0, i, j] = fs_up[k]
            fs[1, i, j] = fs_dn[k]

    return kx, ky, fs
