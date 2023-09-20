# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-09-20

from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial
from scipy import optimize
import h5py
from .directory import walk_emtodirs


def extract_alat_etots(root):
    sws, alat, etot = list(), list(), list()
    for folder in walk_emtodirs(root):
        prn = folder.prn
        if prn is None or not prn.converged:
            continue
        _sws, _, _alat = prn.get_lattice_constants()
        _etot = prn.get_etot()
        sws.append(_sws)
        alat.append(_alat)
        etot.append(_etot)
    idx = np.argsort(sws)
    sws = np.array(sws)[idx]
    alat = np.array(alat)[idx]
    etot = np.array(etot)[idx]
    return sws, alat, etot


def fit_etots(x, etot, deg=3):
    poly = Polynomial.fit(x, etot, deg=deg)
    sol = optimize.minimize(poly, x0=3)
    assert sol.success
    return poly, sol.x


def save_alat_etots(root, deg=3, filename="sws.hdf5"):
    root = Path(root)
    with h5py.File(root / filename, "w") as h5:
        for folder in root.iterdir():
            if folder.is_dir():
                sws, alat, etot = extract_alat_etots(folder)
                ds = h5.create_dataset(folder.name, data=np.array([sws, alat, etot]))
                poly_sws, sws_opt = fit_etots(sws, etot, deg=deg)
                poly_alat, alat_opt = fit_etots(alat, etot, deg=deg)
                ds.attrs["sws_opt"] = sws_opt
                ds.attrs["alat_opt"] = alat_opt
                ds.attrs["poly_sws"] = poly_sws.coef
                ds.attrs["poly_alat"] = poly_alat.coef


def read_data(root, key, quantity="sws", filename="sws.hdf5"):
    with h5py.File(root / filename, "r") as h5:
        ds = h5[key]
        if quantity == "sws":
            x = ds[0]
            etot = ds[2]
            coeffs = ds.attrs["poly_sws"]
            popt = ds.attrs["sws_opt"]
        elif quantity == "alat":
            x = ds[1]
            etot = ds[2]
            coeffs = ds.attrs["poly_alat"]
            popt = ds.attrs["alat_opt"]
        else:
            raise ValueError(f"Unknown quantity: {quantity}")

        poly = Polynomial(coeffs, domain=[x[0], x[-1]])
        return x, etot, poly, popt


def read_opt_latt(root, filename="sws.hdf5"):
    root = Path(root)
    sws = dict()
    alat = dict()
    with h5py.File(root / filename, "r") as h5:
        for key in h5.keys():
            ds = h5[key]
            sws_opt = ds.attrs["sws_opt"]
            alat_opt = ds.attrs["alat_opt"]
            sws[key] = sws_opt[0]
            alat[key] = alat_opt[0]
    return sws, alat
