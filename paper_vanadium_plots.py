# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-08-24

import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mplstyles import use_mplstyle, mplstyle_context, colors
from emtolib import Path, CONFIG, walk_emtodirs, EmtoDirectory, elements
from emtolib.mcmillan import phonon_coupling, mcmillan

ry2ev = constants.value("Rydberg constant times hc in eV")  # 13.605693122994

FIGS = CONFIG["figs"]
ROOT = CONFIG["app"] / "V"
Vanadium = elements["V"]
DEBYE_V = Vanadium.debye_0K
MASS_V = Vanadium.mass


def _extract_tc_u(root):
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


def plot_tc_u(save=False):
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig, ax = plt.subplots()

    root = ROOT / "nl3_57"
    uu, tc, sws = _extract_tc_u(root)
    ax.plot(uu, tc, "-o", label=r"$a=2.9958 \AA$")

    root = ROOT / "nl3_57_2"
    uu, tc, sws = _extract_tc_u(root)
    ax.plot(uu, tc, "-o", label=r"$a=3.0233 \AA$")

    ax.set_xlabel("U (eV)")
    ax.set_ylabel("Tc (K)")
    # ax.set_ylim(3.15, 4.1)
    ax.set_xlim(-0.1, 6.1)

    ax.axhline(5.3, ls="--", color="r", label="Exp.")
    ax.grid()
    ax.legend()
    if save:
        fig.savefig(FIGS / "vanadium_tc_u.png", dpi=900)


def plot_dos_2(save=False):
    xlim = -8, +8
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    data_xps = np.loadtxt(Path("exp", "XPS_png.dat"))

    folder = EmtoDirectory(ROOT, "nl3_57_2", "u23")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()

    fig, ax = plt.subplots()
    ax.plot(energy * ry2ev, dos, label="DFT+DMFT")

    x, y = data_bis.T
    ax.plot(x, y / 3, "o", ms=1.5, color="k", label="BIS + XPS (a.u.)")
    x, y = data_xps.T
    ax.plot(x, y / 3, "o", ms=1.5, color="k")

    ax.set_xlabel("$E - E_F$ (eV)")
    ax.set_ylabel("TDOS (states/eV)")
    ax.axvline(0, color="k", ls="--", lw=0.5)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 35)

    ax.legend()
    if save:
        fig.savefig(FIGS / "vanadium_dos_u=23.png", dpi=900)


def addtext(ax, x, y, text, ha="center", va="center"):
    ax.text(x, y, text, ha=ha, va=va, transform=ax.transAxes)


def plot_dos(save=False):
    xlim = -8, +8
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig = plt.figure(figsize=[3.375, 1.0 * 2.531])
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 2])
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax1.grid(axis="x")
    ax2.grid(axis="x")
    ax3.set_axisbelow(True)
    ax3.grid(axis="x", zorder=-1)
    ax3.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    data_xps = np.loadtxt(Path("exp", "XPS_png.dat"))
    x, y = data_xps.T
    ax1.plot(x, y, "o", ms=1., color="k", label="BIS + XPS (a.u.)")
    addtext(ax1, 0.05, 0.9, "XPS", ha="left", va="top")

    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    x, y = data_bis.T
    ax2.plot(x, y, "o", ms=1., color="k")
    addtext(ax2, 0.95, 0.9, "BIS", ha="right", va="top")

    # Panel 2

    root = ROOT / "sws_opt"

    names = ["u00", "u20", "u40"]
    labels = ["U=0.0 eV", "U=2.0 eV", "U=4.0 eV", "U=3.0 eV"]
    cols = ["C0", "C1", "C3", "C4"]
    for name, label, col in zip(names, labels, cols):
        folder = EmtoDirectory(root, name)
        dosfile = folder.dos
        energy, dos = dosfile.get_total_dos()
        ax3.plot(energy * ry2ev, dos, lw=0.7, color=col, label=label)

    # Styling

    ax1.set_xlim(xlim[0], 0)
    ax2.set_xlim(0, xlim[1])
    ax3.set_xlim(*xlim)

    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 43)
    # ax4.set_ylim(0, 45)

    ax1.set_ylabel("Exp. (a.u.)")
    ax3.set_xlabel("$E - E_F$ (eV)")
    ax3.set_ylabel("DOS (states/eV)")
    # ax4.set_ylabel("DOS (states/eV)")

    # addtext(ax3, 0.95, 0.9, r"$a=2.9958 \AA$", ha="right", va="top")
    # addtext(ax4, 0.95, 0.9, r"$a=3.0233 \AA$", ha="right", va="top")

    ax3.legend(loc="upper left", frameon=True, fontsize="6", )

    if save:
        fig.savefig(FIGS / "vanadium_dos.png", dpi=900)


def plot_dos_3(save=False):
    xlim = -8, +8
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig = plt.figure(figsize=[3.375, 1.3 * 2.531])
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 2, 2])
    gs.update(left=0.15, bottom=0.10, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2, :])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=0)
    ax4.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=0)

    data_xps = np.loadtxt(Path("exp", "XPS_png.dat"))
    x, y = data_xps.T
    ax1.plot(x, y, "o", ms=1., color="k", label="BIS + XPS (a.u.)")
    addtext(ax1, 0.05, 0.9, "XPS", ha="left", va="top")

    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    x, y = data_bis.T
    ax2.plot(x, y, "o", ms=1., color="k")
    addtext(ax2, 0.95, 0.9, "BIS", ha="right", va="top")

    # Panel 2

    folder = EmtoDirectory(ROOT, "nl3_57", "u10")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=1.0 eV")

    folder = EmtoDirectory(ROOT, "nl3_57", "u23")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=2.3 eV")

    folder = EmtoDirectory(ROOT, "nl3_57", "u40")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=4.0 eV", color="C3")

    # Panel 3

    folder = EmtoDirectory(ROOT, "nl3_57", "u10")
    dosfile = folder.dos
    energy, dos1 = dosfile.get_total_dos()
    folder = EmtoDirectory(ROOT, "nl3_57_2", "u10")
    dosfile = folder.dos
    energy, dos2 = dosfile.get_total_dos()
    ax4.plot(energy * ry2ev, np.abs(dos1 - dos2), lw=0.7, label="U=1.0 eV")

    folder = EmtoDirectory(ROOT, "nl3_57", "u23")
    dosfile = folder.dos
    energy, dos1 = dosfile.get_total_dos()
    folder = EmtoDirectory(ROOT, "nl3_57_2", "u23")
    dosfile = folder.dos
    energy, dos2 = dosfile.get_total_dos()
    ax4.plot(energy * ry2ev, np.abs(dos1 - dos2), lw=0.7, label="U=2.3 eV")

    folder = EmtoDirectory(ROOT, "nl3_57", "u40")
    dosfile = folder.dos
    energy, dos1 = dosfile.get_total_dos()
    folder = EmtoDirectory(ROOT, "nl3_57_2", "u40")
    dosfile = folder.dos
    energy, dos2 = dosfile.get_total_dos()
    ax4.plot(energy * ry2ev, np.abs(dos1 - dos2), lw=0.7, label="U=4.0 eV", color="C3")

    ax1.set_xlim(xlim[0], 0)
    ax2.set_xlim(0, xlim[1])
    ax3.set_xlim(*xlim)
    ax4.set_xlim(*xlim)
    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 45)
    # ax4.set_ylim(0, 45)

    ax1.set_ylabel("Exp. (a.u.)")
    ax4.set_xlabel("$E - E_F$ (eV)")
    ax3.set_ylabel("DOS (states/eV)")
    ax4.set_ylabel("DOS (states/eV)")

    addtext(ax3, 0.95, 0.9, r"$a=2.9958 \AA$", ha="right", va="top")
    addtext(ax4, 0.95, 0.9, r"$a=3.0233 \AA$", ha="right", va="top")

    ax1.grid(axis="x")
    ax2.grid(axis="x")
    ax3.grid(axis="x")
    ax4.grid(axis="x")
    ax3.legend(frameon=True)

    if save:
        fig.savefig(FIGS / "vanadium_dos_w_err.png", dpi=900)


def plot_dos_4(save=False):
    xlim = -8, +8
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig = plt.figure(figsize=[3.375, 1.2 * 2.531])
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 2])
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax1.grid(axis="x")
    ax2.grid(axis="x")
    ax3.set_axisbelow(True)
    ax3.grid(axis="x", zorder=-1)
    ax3.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    data_xps = np.loadtxt(Path("exp", "XPS_png.dat"))
    x, y = data_xps.T
    ax1.plot(x, y, "o", ms=1., color="k", label="BIS + XPS (a.u.)")
    addtext(ax1, 0.05, 0.9, "XPS", ha="left", va="top")

    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    x, y = data_bis.T
    ax2.plot(x, y, "o", ms=1., color="k")
    addtext(ax2, 0.95, 0.9, "BIS", ha="right", va="top")

    # Panel 2

    root = ROOT / "nl3_57"
    folder = EmtoDirectory(root, "u10")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=1.0 eV")

    folder = EmtoDirectory(root, "u20")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=2.0 eV")

    folder = EmtoDirectory(root, "u30")
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    ax3.plot(energy * ry2ev, dos, lw=0.7, label="U=3.0 eV", color="C3")

    # Errors
    root2 = ROOT / "nl3_57_2"

    folder1 = EmtoDirectory(root, "u10")
    folder2 = EmtoDirectory(root2, "u10")
    dosfile = folder1.dos
    energy, dos1 = dosfile.get_total_dos()
    dosfile = folder2.dos
    _, dos2 = dosfile.get_total_dos()
    err1 = np.abs(dos1 - dos2)

    folder1 = EmtoDirectory(root, "u20")
    folder2 = EmtoDirectory(root2, "u20")
    dosfile = folder1.dos
    energy, dos1 = dosfile.get_total_dos()
    dosfile = folder2.dos
    _, dos2 = dosfile.get_total_dos()
    err2 = np.abs(dos1 - dos2)

    folder1 = EmtoDirectory(root, "u30")
    folder2 = EmtoDirectory(root2, "u30")
    dosfile = folder1.dos
    energy, dos1 = dosfile.get_total_dos()
    dosfile = folder2.dos
    _, dos2 = dosfile.get_total_dos()
    err3 = np.abs(dos1 - dos2)

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.22, 0.5, 0.3, 0.15]
    ax4 = fig.add_axes([left, bottom, width, height])
    ax4.plot(energy * ry2ev, err1, color="C0", lw=0.5)
    ax4.plot(energy * ry2ev, err2, color="C1", lw=0.5)
    ax4.plot(energy * ry2ev, err3, color="C3", lw=0.5)
    ax4.set_xlim(*xlim)
    ax4.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax4.tick_params(axis='y', which='major', labelsize=5)
    ax4.set_xticklabels([])
    ax4.set_ylim(0, 10)
    # Styling

    ax1.set_xlim(xlim[0], 0)
    ax2.set_xlim(0, xlim[1])
    ax3.set_xlim(*xlim)

    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 43)
    # ax4.set_ylim(0, 45)

    ax1.set_ylabel("Exp. (a.u.)")
    ax3.set_xlabel("$E - E_F$ (eV)")
    ax3.set_ylabel("DOS (states/eV)")
    # ax4.set_ylabel("DOS (states/eV)")

    # addtext(ax3, 0.95, 0.9, r"$a=2.9958 \AA$", ha="right", va="top")
    # addtext(ax4, 0.95, 0.9, r"$a=3.0233 \AA$", ha="right", va="top")

    ax3.legend(loc="upper right", frameon=True, fontsize="6", )

    if save:
        fig.savefig(FIGS / "vanadium_dos_two_sws.png", dpi=900)


def plot_sigma_z(save=False):
    print("---- Sigma(z) ----")
    xlim = -4, +4
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    folder = EmtoDirectory(ROOT, "nl3_57_2", "u20")
    z, sig_z = folder.get_sigma_z()
    z *= ry2ev
    sig_z *= ry2ev

    # Linear fit of Re Σ(0) for quasiparticle weight
    idx = np.argmin(np.abs(z))
    i1 = idx - 3
    i2 = idx + 3
    assert i1 < i2
    deriv = (sig_z[:, :, i2] - sig_z[:, :, i1]).real / (z[i2] - z[i1])
    meff = 1 - deriv

    print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.6f}")
    print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.6f}")

    domain = [-3, 3]
    lin_t2g = Polynomial([sig_z[0, 0, idx].real, deriv[0, 0]])
    lin_eg = Polynomial([sig_z[0, 2, idx].real, deriv[0, 2]])

    # Fit Imaginary part of Σ(0)

    i0, i1 = idx-3, idx+3
    poly2_t2g = Polynomial.fit(z[i0:i1], sig_z[0, 0, i0:i1].imag, 2)
    poly2_eg = Polynomial.fit(z[i0:i1], sig_z[0, 2, i0:i1].imag, 2)

    fig = plt.figure(figsize=[3.375, 1.3 * 2.531])
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    gs.update(left=0.15, wspace=0.025, hspace=0.02)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax1.set_xticklabels([])
    ax1.set_axisbelow(True)
    ax1.grid(axis="x", zorder=-1)
    ax2.set_axisbelow(True)
    ax2.grid(axis="x", zorder=-1)

    ax1.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax1.axhline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax2.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax2.axhline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    ax1.plot(z, sig_z.imag[0, 0], label=r"$t_{2g}$", zorder=2)
    ax1.plot(z, sig_z.imag[0, 2], label=r"$e_{g}$", zorder=2)
    ax1.plot(*poly2_t2g.linspace(domain=domain), ls="--", lw=0.7, color="C0", zorder=1)
    ax1.plot(*poly2_eg.linspace(domain=domain), ls="--", lw=0.7, color="C1", zorder=1)

    ax1.legend(frameon=True)

    ax2.plot(z, sig_z.real[0, 0], label=r"$t_{2g}$", zorder=2)
    ax2.plot(*lin_t2g.linspace(domain=domain), ls="--", lw=0.7, color="C0", zorder=1)

    ax2.plot(z, sig_z.real[0, 2], label=r"$e_{g}$", zorder=2)
    ax2.plot(*lin_eg.linspace(domain=domain), ls="--", lw=0.7, color="C1", zorder=1)

    ax1.set_ylabel(r"Im $\Sigma$ (eV)")
    ax2.set_ylabel(r"Re $\Sigma$ (eV)")
    ax2.set_xlabel("$E - E_F$ (eV)")
    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax1.set_ylim(-0.85, 0.05)
    ax2.set_ylim(-0.6, 0.45)

    fig.subplots_adjust(wspace=0, hspace=0)
    if save:
        fig.savefig(FIGS / "vanadium_selfz.png", dpi=900)


def plot_sigma_iw(save=False):
    print("---- Sigma(iω) ----")
    xlim = 0, +10
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    folder = EmtoDirectory(ROOT, "nl3_57_2", "u20")
    iw, sig_iw = folder.get_sigma_iw(unit="ev")
    iw *= ry2ev
    sig_iw = sig_iw[0]

    deriv = (sig_iw[:, :, 1] - sig_iw[:, :, 0]).imag / (iw[1] - iw[0])
    meff = 1 - deriv
    print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.6f}")
    print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.6f}")

    y0 = sig_iw.imag[:, :, 0] - deriv * iw[0]
    poly_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    poly_eg = Polynomial([y0[0, 2], deriv[0, 2]])

    fig, ax = plt.subplots()

    ax.set_xlabel(r"$i \omega_n$ (eV)")
    ax.set_ylabel(r"Im $\Sigma$ (eV)")
    ax.plot(iw, sig_iw.imag[0, 0], "-o", label=r"$t_{2g}$", ms=2, zorder=2)
    ax.plot(iw, sig_iw.imag[0, 2], "-o", label=r"$e_{g}$", ms=2, zorder=2)
    ax.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    ax.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)

    ax.set_xlim(*xlim)
    ax.set_ylim(-0.38, 0.01)
    ax.legend(frameon=True)

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.38, 0.6, 0.4, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(iw, sig_iw.imag[0, 0], "-o", label=r"$t_{2g}$", ms=2, zorder=2)
    ax2.plot(iw, sig_iw.imag[0, 2], "-o", label=r"$e_{g}$", ms=2, zorder=2)
    ax2.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    ax2.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)
    ax2.set_xlim(0, 0.6)
    ax2.set_ylim(-0.1, 0)

    if save:
        fig.savefig(FIGS / "vanadium_selfiw.png", dpi=900)


def plot_sws_lambda(save=False):
    xlim = -0.1, 5.1
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    root = ROOT / "sws_opt"

    uu = list()
    sws = list()
    alat = list()
    etas = list()
    for folder in walk_emtodirs(root):
        dat = folder.dat
        prn = folder.prn
        if prn is None or not prn.converged:
            continue
        v = dat.get_atom("V")
        u = v.get_u(2)
        if u > 5.0:
            continue
        _sws, _, _alat = prn.get_lattice_constants()
        hopfield = prn.get_sublat_hopfields(dat).sum(axis=1)[0]

        uu.append(u)
        sws.append(_sws)
        alat.append(_alat)
        etas.append(hopfield)

    idx = np.argsort(uu)
    uu = np.array(uu)[idx]
    sws = np.array(sws)[idx]
    alat = np.array(alat)[idx]
    etas = np.array(etas)[idx]

    mass = Vanadium.mass
    theta = Vanadium.debye_0K
    lambd = phonon_coupling(2 * etas, mass, theta)

    # Plot

    fig = plt.figure(figsize=[3.375, 1.3 * 2.531])
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 2])
    gs.update(left=0.15, bottom=0.10, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # ax1.axhline(3.024, color="k", ls="-.", lw=0.5, label="$a$ (exp)")
    ax1.plot(uu, alat, "o-", ms=3, lw=0.8)
    ax2.plot(uu, lambd, "o-", ms=3, lw=0.8)

    def label(mu):
        return r"$\mu^* = " + f"{mu:.2f}" + "$"

    mustar = 0.13
    tc = mcmillan(theta, lambd, mu_star=mustar)
    ax3.plot(uu, tc, "o-", ms=3, lw=0.8, label=label(mustar))

    mustar = 0.14
    tc = mcmillan(theta, lambd, mu_star=mustar)
    ax3.plot(uu, tc, "o-", ms=3, lw=0.8, label=label(mustar))

    mustar = 0.15
    tc = mcmillan(theta, lambd, mu_star=mustar)
    # ax3.plot(uu, tc, "o-", ms=4, lw=0.8, color="C3", label=label(mustar))

    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax3.set_xlim(*xlim)
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax1.set_ylabel(r"$a_{eq}$ ($\AA$)")
    ax2.set_ylabel(r"$\lambda$ (a.u.)")
    ax3.set_ylabel(r"$T_c$ (K)")
    ax3.set_xlabel(r"$U$ (eV)")

    ax1.grid(axis="x")
    ax2.grid(axis="x")
    ax3.grid(axis="x")

    ax3.axhline(5.38, color="red", ls="--", lw=0.5, zorder=0, label="exp")
    ax3.set_ylim(4.5, 6.5)
    ax3.legend(frameon=True)

    if save:
        fig.savefig(FIGS / "V_alat_lambda_tc.png", dpi=900)


def effective_mass_iw(iw, sig_iw):
    """Calculate the effective mass from the imaginary part of the self-energy."""
    deriv = (sig_iw[..., 1] - sig_iw[..., 0]).imag / (iw[1] - iw[0])
    return 1 - deriv


def plot_sigma_iw_temp(save=False):
    s, t2g, eg = 0, 0, 2
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    xlim = 0, 3
    ylim = -0.3, 0

    fig = plt.figure(figsize=[1.2 * 3.375, 2.531])
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.05, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_yticklabels([])

    addtext(ax1, 0.5, 0.9, "t$_{2g}$")
    addtext(ax2, 0.5, 0.9, "e$_{g}$")
    temps = 100, 200, 300, 400
    colors = "C0", "C1", "C3", "C4"
    markers = "o", "s", "D", "v"
    for temp, col, mark in zip(temps, colors, markers):
        folder = EmtoDirectory(ROOT, f"CPA_{temp}K", "u20")
        iw, sig_iw = folder.get_sigma_iw(unit="ev")
        sig_iw = sig_iw[0]

        kwargs = dict(marker=mark, ms=1.5, color=col, lw=0)
        ax1.plot(iw, sig_iw[s, t2g].imag, label=f"{temp}K", **kwargs)
        ax2.plot(iw, sig_iw[s, eg].imag, label=f"{temp}K", **kwargs)

    ax1.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax2.set_xlim(*xlim)
    ax2.set_ylim(*ylim)
    ax1.set_xlabel(r"$i \omega_n$ (eV)")
    ax2.set_xlabel(r"$i \omega_n$ (eV)")
    ax1.set_ylabel(r"Im $\Sigma$ (eV)")
    ax1.legend(frameon=True)
    ax2.legend(frameon=True)
    fig.savefig("sig_iw.png", dpi=900)
    plt.show()


def load_etots(u):
    return np.loadtxt(ROOT / "sws" / f"etot_u{int(u*10):02d}.dat", skiprows=1).T


# noinspection PyTypeChecker,PyUnresolvedReferences
def plot_alat_opt_curves(save=False):
    use_mplstyle("figure", "aps", color_cycle="tableau-colorblind")

    fig, ax = plt.subplots(figsize=[3.375, 0.8 * 2.531])

    u = 2

    sws, alat, etot = load_etots(u)
    poly = Polynomial.fit(alat, etot, deg=3)
    sol = optimize.minimize(poly, x0=3)
    print(sol.x)
    ax.plot(alat, etot, "o", ms=3, label="data", zorder=2)
    ax.plot(*poly.linspace(), color="k", label="poly-fit", zorder=1)
    ax.set_xlabel("$a$ ($Å$)")
    ax.set_ylabel("Total energy (Ry)")
    ax.axvline(sol.x, color="r", ls="--", lw=0.5, label="$a_{eq}$")
    ax.axvline(3.024, color="k", ls="-.", lw=0.5, label="$a$ (exp)")
    ax.legend(loc="upper right", frameon=True)
    if save:
        fig.savefig(FIGS / "vanadium_alat_etot.png", dpi=900)


def main():
    save = True
    # plot_tc_u(save=save)
    # plot_dos(save=save)
    # plot_dos(save=save)
    # plot_dos_3(save=save)
    # plot_dos_4(save=save)
    # plot_sigma_z(save=save)
    # plot_sigma_iw(save=save)
    # plot_sws_lambda(save=save)
    # plot_alat_opt_curves(save=save)
    plt.show()


if __name__ == "__main__":
    main()
