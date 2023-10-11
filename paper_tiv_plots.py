# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-08-21

import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mplstyles import use_mplstyle, mplstyle_context
from emtolib import Path, CONFIG, walk_emtodirs, EmtoDirectory, elements
from emtolib.mcmillan import phonon_coupling, mcmillan
from emtolib.meff import deriv_iw, effective_mass

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
    colors = ["C3", "C0", "C1"]
    for u, c in zip([0, 2, 4], colors):
        root = ROOT / "CPA" / f"u{u}_400K"  # / f"nl3_u{u}"
        cc, tc = _extract_tc_c_tiv(root, ckey)
        ax.plot(cc, tc, "-o", lw=0.8, ms=2, color=c, label=f"U={u} eV")

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
    ax.plot(x, y, "-o", lw=0.8, ms=2, color="C2", label="exp", zorder=1)

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


def plot_sigma_iw(save=False):
    print("---- Sigma(iω) ----")
    inset = False
    u = 2
    c = 0.35
    root = ROOT / "CPA"
    xlim = 0, +10
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    folder = EmtoDirectory(root, f"u{u}_400K", f"Ti{int(c*100)}")
    iw, sig_iw = folder.get_sigma_iw(unit="ev")
    sig_iw_tot = np.average(sig_iw, weights=[c, 1 - c], axis=0)

    deriv = deriv_iw(iw, sig_iw_tot)
    meff = effective_mass(deriv)
    print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.3f}")
    print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.3f}")

    y0 = sig_iw_tot.imag[:, :, 0] - deriv * iw[0]
    poly_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    poly_eg = Polynomial([y0[0, 2], deriv[0, 2]])

    fig, ax = plt.subplots()

    ax.set_xlabel(r"$i \omega_n$ (eV)")
    ax.set_ylabel(r"Im $\Sigma$ (eV)")
    ax.plot(iw, sig_iw_tot.imag[0, 0], "-o", color="C0", label=r"$t_{2g}$", ms=2, zorder=2)
    ax.plot(iw, sig_iw_tot.imag[0, 2], "-o", color="C1", label=r"$e_{g}$", ms=2, zorder=2)
    c = "C0"
    ax.plot(iw, sig_iw.imag[0, 0, 0], "s", color=c, label=r"Ti $t_{2g}$", ms=1., zorder=2)
    ax.plot(iw, sig_iw.imag[1, 0, 0], "D", color=c, label=r"V $t_{2g}$", ms=1., zorder=2)
    c = "C1"
    ax.plot(iw, sig_iw.imag[0, 0, 2], "s", color=c, label=r"Ti $e_{g}$", ms=1., zorder=2)
    ax.plot(iw, sig_iw.imag[1, 0, 2], "D", color=c, label=r"V $e_{g}$", ms=1., zorder=2)
    ax.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    ax.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)

    ax.set_xlim(*xlim)
    ax.set_ylim(-0.42, 0.01)
    ax.legend(frameon=True)
    if inset:
        # These are in unitless percentages of the figure size. (0,0 is bottom left)
        left, bottom, width, height = [0.38, 0.6, 0.4, 0.3]
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.plot(iw, sig_iw_tot.imag[0, 0], "o", label=r"$t_{2g}$", ms=2, zorder=2)
        ax2.plot(iw, sig_iw_tot.imag[0, 2], "o", label=r"$e_{g}$", ms=2, zorder=2)
        ax2.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
        ax2.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)
        ax2.set_xlim(0, 0.6)
        ax2.set_ylim(-0.1, 0)

    if save:
        fig.savefig(FIGS / "TiV_selfiw.png", dpi=900)


def plot_sigma_iw2(save=False):
    print("---- Sigma(iω) ----")
    inset = False
    u = 2
    c = 0.35
    root = ROOT / "CPA"
    xlim = 0, +4.5
    ylim = -0.42, 0.01
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    folder = EmtoDirectory(root, f"u{u}_400K", f"Ti{int(c*100)}")
    iw, sig_iw = folder.get_sigma_iw(unit="ev")
    sig_iw_tot = np.average(sig_iw, weights=[c, 1 - c], axis=0)

    deriv = deriv_iw(iw, sig_iw_tot)
    meff = effective_mass(deriv)
    print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.3f}")
    print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.3f}")

    y0 = sig_iw_tot.imag[:, :, 0] - deriv * iw[0]
    poly_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    poly_eg = Polynomial([y0[0, 2], deriv[0, 2]])

    fig = plt.figure()  # figsize=[3.375, 1.0 * 2.531])
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_yticklabels([])
    ax1.text(0.3, 0.95, "t$_{2g}$", transform=ax1.transAxes, ha="center", va="center")
    ax2.text(0.3, 0.95, "e$_{g}$", transform=ax2.transAxes, ha="center", va="center")
    # plot_tc_conc_tiv(save)
    ax1.set_xlabel(r"$i \omega_n$ (eV)")
    ax2.set_xlabel(r"$i \omega_n$ (eV)")
    ax1.set_ylabel(r"Im $\Sigma$ (eV)")

    ax1.plot(iw, sig_iw_tot.imag[0, 0], "-o", color="C0", label=r"Total", ms=2, zorder=2)
    ax2.plot(iw, sig_iw_tot.imag[0, 2], "-o", color="C1", label=r"Total", ms=2, zorder=2)
    ax1.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    ax2.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)

    ax1.plot(iw, sig_iw.imag[0, 0, 0], "s", color="C0", label=r"Ti", ms=2., zorder=2)
    ax2.plot(iw, sig_iw.imag[0, 0, 2], "s", color="C1", label=r"Ti", ms=2., zorder=2)

    ax1.plot(iw, sig_iw.imag[1, 0, 0], "x", color="C0", label=r"V", ms=2., zorder=2)
    ax2.plot(iw, sig_iw.imag[1, 0, 2], "x", color="C1", label=r"V", ms=2., zorder=2)

    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax2.set_ylim(*ylim)
    ax1.grid(axis="y")
    ax2.grid(axis="y")
    ax1.legend(frameon=True)
    ax2.legend(frameon=True)
    if save:
        fig.savefig(FIGS / "TiV_selfiw2.png", dpi=900)


def extract_meffs(root):
    s, t2g, eg = 0, 0, 2
    cc = []
    meffs = []
    for folder in walk_emtodirs(root):
        print(folder)
        dat = folder.dat
        c = dat.get_concentration("Ti")
        try:
            if c > 0.7:
                continue
            iw, sig_iw = folder.get_sigma_iw(unit="ev")
            if c == 0:
                sig_iw = sig_iw[0]
            else:
                sig_iw = np.average(sig_iw, weights=[c, 1 - c], axis=0)
            deriv = deriv_iw(iw, sig_iw)
            meff = 1 - deriv
            print(meff[s])

            cc.append(c)
            meffs.append([meff[s, t2g], meff[s, eg]])
        except FileNotFoundError:
            continue

    idx = np.argsort(cc)
    cc = np.array(cc)[idx]
    meffs = np.array(meffs)[idx]
    return cc, meffs


def plot_meff(save=False):
    print("---- m*(U) ----")
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    root = ROOT / "CPA"

    fig, ax = plt.subplots(figsize=[3.375, 0.75 * 2.531])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$m^* / m$")

    u = 2
    cc, meffs = extract_meffs(root / f"u{u}_400K")
    meff_total = (3 * meffs[:, 0] + 2 * meffs[:, 1]) / 5
    ax.plot(cc, meff_total, "o--", color="C0", ms=2, label=f"$U={u}$")

    u = 4
    cc, meffs = extract_meffs(root / f"u{u}_400K")
    meff_total = (3 * meffs[:, 0] + 2 * meffs[:, 1]) / 5
    #ax.plot(cc, meff_total, "s-.", color="C1", ms=2, label=f"$U={u}$")

    # uu, meffs = extract_meffs(root / "400K")
    # ax.plot(uu, meffs[:, 0], "s", color="C0", ms=2, label="t2g")
    # ax.plot(uu, meffs[:, 1], "s", color="C1", ms=2, label="eg")
    # meff_total = (3 * meffs[:, 0] + 2 * meffs[:, 1]) / 5
    # ax.plot(uu, meff_total, "s-.", color="C1", ms=2, label="400K")
    # ax.set_xlim(0.45, 4.05)
    #ax.set_ylim(1, 2)
    ax.legend(frameon=True)
    ax.grid(axis="y")
    if save:
        fig.savefig(FIGS / "TiV_c_meff.png", dpi=900)


def main():
    save = True

    # plot_dos_cpa(save)
    # plot_sigma_iw(save)
    # plot_meff(save)
    plot_sigma_iw2(save)
    plt.show()


if __name__ == "__main__":
    main()
