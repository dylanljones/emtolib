# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-08-21

import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mplstyles import use_mplstyle, mplstyle_context
from emtolib import Path, CONFIG, walk_emtodirs, EmtoDirectory, elements
from emtolib.mcmillan import phonon_coupling, mcmillan

ry2ev = constants.value("Rydberg constant times hc in eV")  # 13.605693122994

ROOT = CONFIG["app"] / "Ti-V"
FIGS = CONFIG["figs"]
Vanadium = elements["V"]
Titan = elements["Ti"]
DEBYE_TI = Titan.debye_0K
DEBYE_V = Vanadium.debye_0K
MASS_TI = Titan.mass
MASS_V = Vanadium.mass


def read_fig_data(path):
    delim = ","
    text = Path(path).read_text()
    lines = text.splitlines()[2:]
    cols = len(lines[0].split(delim))
    if cols % 2 != 0:
        raise ValueError("Number of columns must be odd")
    dsets = [list() for _ in range(cols // 2)]

    for line in lines:
        vals = line.split(delim)
        for i in range(len(vals) // 2):
            xval = vals[2 * i + 0]
            yval = vals[2 * i + 1]
            if xval and yval:
                dsets[i].append((float(xval), float(yval)))

    dsets = [np.array(dset) for dset in dsets]
    return dsets


def _extract_tc_u_v(root):
    uu, etas = [], []
    sws = 0
    for folder in walk_emtodirs(root):
        dat = folder.dat
        sws = dat.sws
        prn = folder.prn
        u = dat.get_atom("V").u[2]
        if u in [2.2, 2.3, 2.4]:
            continue
        uu.append(u)
        hopfield = prn.get_sublat_hopfields(dat).sum(axis=1)[0]
        etas.append(hopfield)

    idx = np.argsort(uu)
    uu = np.array(uu)[idx]
    etas = np.array(etas)[idx]

    mass = Vanadium.mass
    theta = Vanadium.debye_0K

    lambd = phonon_coupling(2 * etas, mass, theta)
    tc = mcmillan(theta, lambd, mu_star=0.13)
    return uu, tc, sws


def _extract_tc_c_tiv(root, ckey):
    cc = list()
    hopfields = list()
    for folder in walk_emtodirs(root):
        dat = folder.dat
        prn = folder.prn
        if prn is None:
            print(f"Skipping {folder.name}")
            continue
        try:
            c = dat.get_concentration(ckey)
            eta = prn.get_sublat_hopfields(dat).sum(axis=1)[0]
            cc.append(c)
            hopfields.append(eta)
        except Exception as e:
            print(f"Error in {folder.name}")
            print(e)
            continue

    idx = np.argsort(cc)
    cc = np.array(cc)[idx]
    hopfields = np.array(hopfields)[idx]

    mass = cc * MASS_V + (1 - cc) * MASS_TI
    debye = cc * DEBYE_V + (1 - cc) * DEBYE_TI
    lamb = phonon_coupling(2 * hopfields, mass, debye)
    t_c = mcmillan(debye, lamb, mu_star=0.13)
    return cc, t_c


def plot_tc_conc_tiv(save=False):
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig, ax = plt.subplots()
    ckey = "Ti"
    for u in [2, 4]:
        root = ROOT / f"nl3_u{u}"
        cc, tc = _extract_tc_c_tiv(root, ckey)
        print(cc, tc)
        ax.plot(cc, tc, "-o", lw=0.8, ms=2, label=f"U={u} eV")

    regions = [1 - 0.335, 1 - 0.145, 1]

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 9)
    ax.set_xlabel(r"x")
    ax.set_ylabel("$T_c$ (K)")
    ax.fill_between([regions[-2] - 0.02, regions[-1]], 0, 9, color="k", alpha=0.1, lw=0)
    ax.fill_between([regions[-3], regions[-2] + 0.02], 0, 9, color="red", alpha=0.1, lw=0)
    x1 = regions[-1] + 0.5 * (regions[-2] - regions[-1])
    x2 = regions[-2] + 0.5 * (regions[-3] - regions[-2])
    y = 0.9
    ax.text(x1, y, r"$\alpha$", ha="center", va="center")
    ax.text(x2, y, r"$\beta + \omega$", ha="center", va="center")
    ax.text(0.55, y, r"$\beta$ (bcc)", ha="center", va="center")
    ax.annotate(
        r"$\alpha + \beta + \omega$",
        xy=(regions[-2], y),
        xytext=(regions[-2], -1.2),
        ha="center",
        arrowprops=dict(width=0.1, lw=0.5, color="k", headwidth=3, headlength=3, shrink=0.0),
    )

    path = Path("exp") / "TiV_tc_Matin_fig1.csv"
    dsets = read_fig_data(path)
    x, y = dsets[1].T
    x /= 100
    x = 1 - x
    mask = x < regions[0]
    # x, y = x[mask], y[mask]
    ax.plot(x, y, "-o", lw=0.8, ms=2, label="exp", zorder=1)

    ax.legend(loc="lower left", frameon=True)
    if save:
        fig.savefig(FIGS / "TiV_conc_tc.png", dpi=900)


def plot_dos_cpa(save=False):
    xlim = -8, +8
    root = ROOT / "CPA" / "u2_400K"
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig, ax = plt.subplots()  # figsize=[3.375, 1.0 * 2.531])
    ax.grid(axis="x", zorder=-1)
    ax.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    names = ["Ti10", "Ti35", "Ti60"]
    labels = ["$x=0.10$", "$x=0.35$", "$x=0.60$"]
    colors = ["C0", "C1", "C3"]

    kwargs = dict(lw=0.7, zorder=2, ls="-")
    for name, label, color in zip(names, labels, colors):
        folder = EmtoDirectory(root / name)
        dosfile = folder.dos
        energy, dos = dosfile.get_total_dos()
        ax.plot(energy * ry2ev, dos, color=color, label=label, **kwargs)

    ax.set_xlabel("$E - E_F$ (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, None)
    ax.legend(frameon=True)
    if save:
        fig.savefig(FIGS / "TiV_DOS_CPA.png", dpi=900)


def main():
    save = False
    plot_tc_conc_tiv(save)
    # plot_dos_cpa(save)
    plt.show()


if __name__ == "__main__":
    main()
