# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from mplstyles import use_mplstyle
from emtolib.xarr import load_dataset, update_datasets

FIGDIR = Path("figs")


def load_datasets(root: Union[Path, str], xarr_dir: str = "xarr"):
    xarr_root = Path(root) / xarr_dir
    datasets = dict()
    for folder in xarr_root.iterdir():
        pressure = int(folder.name.split("_")[0])
        ds = load_dataset(folder)
        datasets[pressure] = ds
    return datasets


def dos_cpa_subplots(xmin=-0.5, xmax=0.5, ymax1=80, ymax2=60):
    # Figure with two subplots in first row and one in second row
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    gs.update(wspace=0.025, hspace=0.05, bottom=0.15, left=0.15, right=0.95, top=0.95)

    ax1.xaxis.set_tick_params(labelbottom=False)
    ax2.xaxis.set_tick_params(labelbottom=False)
    ax2.yaxis.set_tick_params(labelleft=False)

    ax1.axvline(0, color="dimgray", ls="--", lw=0.5)
    ax2.axvline(0, color="dimgray", ls="--", lw=0.5)
    ax3.axvline(0, color="dimgray", ls="--", lw=0.5)

    ax1.set_xlim(xmin, xmax)
    ax2.set_xlim(xmin, xmax)
    ax3.set_xlim(xmin, xmax)
    ax1.set_ylim(0, ymax1)
    ax2.set_ylim(0, ymax1)
    ax3.set_ylim(0, ymax2)
    return fig, (ax1, ax2, ax3)


def table(datasets):
    ds100 = datasets[1]
    ds0 = datasets[0]
    ds25 = datasets[0.25]
    ds50 = datasets[0.50]
    ds75 = datasets[0.75]

    debye = [
        ds0.attrs["debye_temp"],
        ds25.attrs["debye_temp"],
        ds50.attrs["debye_temp"],
        ds75.attrs["debye_temp"],
        ds100.attrs["debye_temp"],
    ]
    hopfield = [
        ds0.hopfields.values[0],
        np.sum(ds25.hopfields.values * ds25.coords["concs"].values),
        np.sum(ds50.hopfields.values * ds50.coords["concs"].values),
        np.sum(ds75.hopfields.values * ds75.coords["concs"].values),
        ds100.hopfields.values[0],
    ]
    print("Debye temperature")
    print(debye)
    print("Hopfield parameters")
    print(hopfield)


def generate_figure(datasets):
    dsv = datasets[0]
    dsnb = datasets[1]
    ds25 = datasets[0.25]
    ds50 = datasets[0.50]
    ds75 = datasets[0.75]

    # Plot DOS
    use_mplstyle("figure", "aps")
    fig, (ax1, ax2, ax3) = dos_cpa_subplots(ymax1=80, ymax2=45)
    lw = 0.8
    ax1.plot(dsnb.energy, dsnb.dos, label="Nb", lw=lw)
    ax2.plot(dsv.energy, dsv.dos, label="V", lw=lw)
    ax3.plot(ds25.energy, ds25.dos, label=r"$c_{Nb}=0.25$", lw=lw)
    ax3.plot(ds50.energy, ds50.dos, label=r"$c_{Nb}=0.50$", lw=lw)
    ax3.plot(ds75.energy, ds75.dos, label=r"$c_{Nb}=0.75$", lw=lw)

    ax3.set_xlabel("$E-E_F$ (Ry)")
    ax3.set_ylabel("States / Ry")
    ax3.yaxis.set_label_coords(-.1, 1.03)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    fig.savefig(FIGDIR / "NbV_dos_cpa.png")


def main():
    root = Path("app") / "Nb-Ti" / "CPA" / "Nb44"
    update_datasets(root, force=False)
    datasets = load_datasets(root)
    print(datasets)

    fig, ax = plt.subplots()
    for p, ds in datasets.items():
        ax.plot(ds.energy, ds.dos, label=f"{p} GPa")
    ax.legend()
    plt.show()



if __name__ == "__main__":
    main()