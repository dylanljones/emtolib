# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
import matplotlib.pyplot as plt
from typing import Union
from mplstyles import use_mplstyle
from scipy import constants as const
from emtolib.xarr import load_dataset, update_datasets

FIGDIR = Path("figs")

ry_to_ev = const.physical_constants["Rydberg constant times hc in eV"][0]


def load_datasets(root: Union[Path, str], xarr_dir: str = "xarr"):
    xarr_root = Path(root) / xarr_dir
    datasets = dict()
    for folder in xarr_root.iterdir():
        ds = load_dataset(folder)
        try:
            i = ds.attrs["atoms"].index("Nb")
            c_nb = ds.coords["concs"].values[i]
        except ValueError:
            c_nb = 0
        datasets[c_nb] = ds
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


def generate_figure(datasets):
    fac = 1  # / (2 * ry_to_ev)
    dsv = datasets[0]
    dsnb = datasets[1]
    ds25 = datasets[0.25]
    ds50 = datasets[0.50]
    ds75 = datasets[0.75]

    # Plot DOS
    use_mplstyle("figure", "aps", color_cycle="seaborn-colorblind")
    fig, (ax1, ax2, ax3) = dos_cpa_subplots(ymax1=80, ymax2=45)
    lw = 0.8
    ax1.plot(dsnb.energy, dsnb.dos * fac, label="Nb", lw=lw)
    ax2.plot(dsv.energy, dsv.dos * fac, label="V", lw=lw)
    ax3.plot(ds25.energy, ds25.dos * fac, label=r"$c_{Nb}=0.25$", lw=lw)
    ax3.plot(ds50.energy, ds50.dos * fac, label=r"$c_{Nb}=0.50$", lw=lw)
    ax3.plot(ds75.energy, ds75.dos * fac, label=r"$c_{Nb}=0.75$", lw=lw)

    ax3.set_xlabel("$E-E_F$ (Ry)")
    # ax3.set_ylabel("States / eV atoms")
    ax3.set_ylabel("States / Ry")
    ax3.yaxis.set_label_coords(-.1, 1.03)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    return fig, (ax1, ax2, ax3)


def main():
    root = Path("app") / "Nb-V" / "CPA2"
    update_datasets(root, force=False)
    datasets = load_datasets(root)
    fig, axs = generate_figure(datasets)
    fig.savefig(FIGDIR / "NbV_dos_cpa.png")
    plt.show()


if __name__ == "__main__":
    main()
