# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

import time
import shutil
from pathlib import Path
import numpy as np
import xarray as xr
from typing import Union, Callable
from .directory import EmtoDirectory, walk_emtodirs
from .common import elements


def save_dataset(ds: xr.Dataset, location: Union[Path, str] = ".", *, name: str, **kw):
    location = Path(location)
    location.mkdir(parents=True, exist_ok=True)
    dir_name = location / name
    if not dir_name.suffix:
        dir_name = dir_name.with_suffix(".zarr")
    ds.to_zarr(dir_name, mode="w", **kw)
    return dir_name


def construct_dataset(folder: EmtoDirectory):
    dat = folder.dat
    # Input/Output file
    prn = folder.get_prn()
    hopfield = prn.get_hopfield()
    dos_ef = prn.get_dos_ef()

    atoms = [atom.symbol for atom in dat.atoms]
    hopfields = np.array([hopfield[atom] for atom in atoms])
    concs = np.array([atom.conc for atom in dat.atoms])
    masses = np.array([elements[atom]["mass"] for atom in atoms])
    debye = np.array([elements[atom]["debye_0K"] for atom in atoms])
    mass_avg = np.sum(concs * masses)
    debye_avg = np.sum(concs * debye)

    # DOS file
    dosfile = folder.get_dos()
    energy, dos = dosfile.get_total_dos()

    data = dict(
        energy=("n", energy),
        dos=("n", dos),
        hopfields=("m", hopfields),
    )
    coords = dict(
        concs=concs,
        masses=masses,
        debye=debye,
    )
    attrs = dict(
        atoms=atoms,
        mass=mass_avg,
        debye_temp=debye_avg,
        dos_ef=dos_ef,
    )

    # Create dataset
    ds = xr.Dataset(data, coords=coords, attrs=attrs)
    return ds


def rmdir(path, maxtry=10):

    e = None
    if path.exists():
        try:
            shutil.rmtree(path)
        except PermissionError:
            pass
        i = 1
        while path.exists():
            time.sleep(0.1)
            try:
                shutil.rmtree(path)
            except PermissionError:
                pass
            if i >= maxtry:
                raise e
            i += 1


def update_datasets(
    src: Union[Path, str], location: Union[Path, str] = None, force=False, exclude=()
):
    emto_root = Path(src)
    xarr_root = emto_root / "xarr" if location is None else Path(location)
    if xarr_root.exists():
        root_mtime = emto_root.stat().st_mtime
        xarr_mtime = xarr_root.stat().st_mtime
        if (xarr_root.exists() and root_mtime <= xarr_mtime) and not force:
            # print("Datasets are up-to date.")
            return xarr_root
        rmdir(xarr_root)
    else:
        xarr_root.mkdir(parents=True)

    print("=" * 50)
    print("Updating datasets", xarr_root)
    print("=" * 50)

    for folder in walk_emtodirs(emto_root):
        if exclude and folder.path.name in exclude:
            continue
        print(folder)
        try:
            ds = construct_dataset(folder)
            save_dataset(ds, xarr_root, name=folder.path.name)
        except Exception as e:
            print("Error", e)

    print("-" * 50)
    return xarr_root


def load_dataset(dir_name: Union[Path, str], **kw):
    if not dir_name.suffix:
        dir_name = dir_name.with_suffix(".zarr")
    return xr.open_zarr(dir_name, **kw)


def load_datasets(location: Union[Path, str], keyfunc: Callable):
    xarr_root = Path(location)
    datasets = dict()
    for folder in xarr_root.iterdir():
        ds = load_dataset(folder)
        key = keyfunc(ds)
        datasets[key] = ds
    return datasets
