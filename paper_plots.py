# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-10-17

import numpy as np
from numpy.polynomial import Polynomial
from scipy import constants
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mplstyles import use_mplstyle, mplstyle_context, set_colorcycle
import h5py
from emtolib import Path, CONFIG, walk_emtodirs, EmtoDirectory, elements
from emtolib.mcmillan import phonon_coupling, mcmillan
from emtolib.sws import read_data
from emtolib.meff import deriv_iw, effective_mass
from emtolib.utils import fermi_fct, gaussian, lorentzian, convolve_func
import seaborn as sns


stylesheets = plt.style.core.read_style_directory("styles")
plt.style.core.update_nested_dict(plt.style.library, stylesheets)

RY2EV = constants.value("Rydberg constant times hc in eV")  # 13.605693122994
KB = constants.value("Boltzmann constant in eV/K")  # 8.617333262145e-05

ROOT = CONFIG["app"]
FIGS = CONFIG["figs"]
Vanadium = elements["V"]
Titan = elements["Ti"]
DEBYE_TI = Titan.debye_0K
DEBYE_V = Vanadium.debye_0K
MASS_TI = Titan.mass
MASS_V = Vanadium.mass

MARKERS = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "h", "H", "X", "d", "D"]
PALETTE = "colorblind"  # "deep"

use_mplstyle("figure", "paper")
mpl.rcParams["axes.prop_cycle"] = plt.cycler(color=sns.color_palette(PALETTE))


def addtext(ax, x, y, text, ha="center", va="center"):
    ax.text(x, y, text, ha=ha, va=va, transform=ax.transAxes)


def extract_meffs_u(root, s=0, t2g=0, eg=2):
    uu = []
    meffs = []
    for folder in walk_emtodirs(root):
        dat = folder.dat
        u = dat.get_atom("V").u[2]
        try:
            iw, sig_iw = folder.get_sigma_iw(unit="ev")
            sig_iw = sig_iw[0]
            deriv = deriv_iw(iw, sig_iw)
            meff = 1 - deriv
            # print(f"Z_t2g = {1 / meff[s, t2g]:.6f}  m* = {meff[s, t2g]:.3f}")
            # print(f"Z_eg  = {1 / meff[s, eg]:.6f}  m* = {meff[s, eg]:.3f}")
            uu.append(u)
            meffs.append([meff[s, t2g], meff[s, eg]])
        except FileNotFoundError:
            continue

    idx = np.argsort(uu)
    uu = np.array(uu)[idx]
    meffs = np.array(meffs)[idx]
    return uu, meffs


def extract_meffs_c(root, s=0, t2g=0, eg=2):
    cc = []
    meffs = []
    for folder in walk_emtodirs(root):

        dat = folder.dat
        c = dat.get_concentration("Ti")
        try:
            if c > 0.7:
                continue
            iw, sig_iw = folder.get_sigma_iw(unit="ev")
            deriv = deriv_iw(iw, sig_iw)
            meff = 1 - deriv
            cc.append(c)
            if c == 0:
                ti = np.full_like(meff[0], fill_value=np.nan)
                meff = np.array([ti, meff[0]])
            meffs.append([meff[:, s, t2g], meff[:, s, eg]])
        except FileNotFoundError:
            continue

    idx = np.argsort(cc)
    cc = np.array(cc)[idx]
    meffs = np.array(meffs)[idx]
    print(meffs.shape)
    return cc, meffs


def extract_tc_c_tiv(root, ckey):
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


def load_data_p():
    file = ROOT / "V" / "DATA_P.DOS"
    dos_data = np.loadtxt(file)
    c = dos_data.shape[0] // 2
    dos_up, dos_dn = dos_data[:c], dos_data[c:]
    energy, total_md, xy, yz, z2, zx, x2my2, t2g, eg = dos_up.T
    return energy, t2g * RY2EV, eg * RY2EV


# noinspection PyTypeChecker,PyUnresolvedReferences
def plot_alat_opt_curves_v(save=False):
    print("---- V: E(a) ----")

    u = 2
    alat, etot, poly, popt = read_data(
        ROOT / "V" / "sws", key=f"u{int(u*10):02d}", quantity="alat"
    )

    fig, ax = plt.subplots(figsize=[3.375, 0.8 * 2.531])
    print(f"Alat opt: {popt[0]}")
    ax.plot(alat, etot, "o", ms=3, label="data", zorder=2)
    ax.plot(*poly.linspace(), color="k", label="poly-fit", zorder=1)
    ax.set_xlabel("$a$ ($Å$)")
    ax.set_ylabel("Total energy (Ry)")
    ax.axvline(popt, color="r", ls="--", lw=0.5, label="$a_{eq}$")
    ax.axvline(3.024, color="k", ls="-.", lw=0.5, label="$a$ (exp)")
    ax.legend(loc="upper right")
    if save:
        fig.savefig(FIGS / "V_alat_etot.png", dpi=900)


def plot_dos_conv_v(save=False):
    print("---- DOS(z) conv ----")
    root = ROOT / "V" / "sws_opt" / "400K"
    xlim = -8, +8

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
    ax1.plot(x, y, "o", ms=1.0, color="k", label="BIS + XPS (a.u.)")
    addtext(ax1, 0.05, 0.9, "XPS", ha="left", va="top")

    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    x, y = data_bis.T
    ax2.plot(x, y, "o", ms=1.0, color="k")
    addtext(ax2, 0.95, 0.9, "BIS", ha="right", va="top")

    # Panel 2
    u = 2
    folder = EmtoDirectory(root, f"u{u}0")
    dat = folder.dat
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    energy *= RY2EV
    temp = dat.ttt
    beta = 1 / (KB * temp)
    # Plot original DOS
    ax3.plot(energy, dos, "-", color="grey", label="DOS", lw=0.5, zorder=1)
    ax3.fill_between(energy, dos, color="grey", alpha=0.15, zorder=1)

    # Apply Fermi-Dirac smearing
    dos_val = fermi_fct(energy, beta=beta) * dos
    dos_con = fermi_fct(-energy, beta=beta) * dos
    # Convolve with a constant Lorentzian broadening of 0.1eV
    dos_val = convolve_func(energy, dos_val, lorentzian, width=0.1, x0=0)
    dos_con = convolve_func(energy, dos_con, lorentzian, width=0.1, x0=0)
    # Convolve VB (CB) with a Gaussian broadening of 0.55eV (0.7eV)
    dos_val = convolve_func(energy, dos_val, gaussian, width=0.55, x0=0)
    dos_con = convolve_func(energy, dos_con, gaussian, width=0.7, x0=0)
    ax3.plot(energy, dos_val, "-", label="VB", color="C0", zorder=2)
    ax3.plot(energy, dos_con, "-", label="CB", color="C1", zorder=2)

    # Styling
    ax1.set_xlim(xlim[0], 0)
    ax2.set_xlim(0, xlim[1])
    ax3.set_xlim(*xlim)
    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 29)
    ax1.set_ylabel("Exp. (a.u.)")
    ax3.set_xlabel("$E - E_F$ (eV)")
    ax3.set_ylabel("DOS (states/eV)")
    ax3.legend(loc="upper left", fontsize="6",)
    if save:
        fig.savefig(FIGS / "V_dos_conv.png", dpi=900)


def plot_dos_conv_v2(save=False):
    print("---- DOS(z) conv ----")
    root = ROOT / "V" / "sws_opt" / "400K"
    xlim = -8, +8

    fig = plt.figure(figsize=[3.375, 1.2 * 2.531])
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 2])
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2, :])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax2.set_yticklabels([])
    ax1.grid(axis="x")
    ax2.grid(axis="x")
    ax3.set_axisbelow(True)
    ax3.grid(axis="x", zorder=-1)
    ax4.set_axisbelow(True)
    ax4.grid(axis="x", zorder=-1)
    ax3.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax4.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    data_xps = np.loadtxt(Path("exp", "XPS_png.dat"))
    x, y = data_xps.T
    ax1.plot(x, y, "o", ms=1.0, color="k", label="BIS + XPS (a.u.)")
    data_bis = np.loadtxt(Path("exp", "BIS_png.dat"))
    x, y = data_bis.T
    ax2.plot(x, y, "o", ms=1.0, color="k")
    addtext(ax1, 0.05, 0.9, "XPS", ha="left", va="top")
    addtext(ax2, 0.05, 0.9, "BIS", ha="left", va="top")

    # Panel 2
    u = 2
    folder = EmtoDirectory(root, f"u{u}0")
    dat = folder.dat
    dosfile = folder.dos
    energy, dos = dosfile.get_total_dos()
    energy *= RY2EV
    temp = dat.ttt
    beta = 1 / (KB * temp)
    # Apply Fermi-Dirac smearing
    dos_val = fermi_fct(energy, beta=beta) * dos
    dos_con = fermi_fct(-energy, beta=beta) * dos
    # Convolve with a constant Lorentzian broadening of 0.1eV
    dos_val = convolve_func(energy, dos_val, lorentzian, width=0.1, x0=0)
    dos_con = convolve_func(energy, dos_con, lorentzian, width=0.1, x0=0)
    # Convolve VB (CB) with a Gaussian broadening of 0.55eV (0.7eV)
    dos_val = convolve_func(energy, dos_val, gaussian, width=0.55, x0=0)
    dos_con = convolve_func(energy, dos_con, gaussian, width=0.7, x0=0)
    ax3.plot(energy, dos_val, "-", label="VB", color="C2", zorder=2)
    ax3.plot(energy, dos_con, "-", label="CB", color="C3", zorder=2)
    # Plot original DOS
    ax4.plot(energy, dos, "-", color="grey", label="DFT+DMFT", lw=0.5, zorder=1)
    ax4.fill_between(energy, dos, color="grey", alpha=0.15, zorder=1)

    # Read LMTO DOS
    energy, dos_t2g, dos_eg = load_data_p()

    ax4.plot(energy, dos_t2g, "-", label="DFT d-t$_{2g}$", color="C0", zorder=2)
    ax4.plot(energy, dos_eg, "-", label="DFT d-e$_g$", color="C1", zorder=2)

    # Styling
    ax1.set_xlim(xlim[0], 0)
    ax2.set_xlim(0, xlim[1])
    ax3.set_xlim(*xlim)
    ax4.set_xlim(*xlim)
    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 23)
    ax4.set_ylim(0, 33)

    ax1.set_ylabel("Exp. (a.u.)")
    ax3.set_xlabel("$E - E_F$ (eV)")
    ax4.set_ylabel("DOS (states/eV)")
    ax3.legend(loc="upper left", fontsize="5",)
    ax4.legend(loc="upper left", fontsize="5",)

    addtext(ax2, 0.95, 0.9, r"(a)", ha="right", va="top")
    addtext(ax3, 0.975, 0.9, r"(b)", ha="right", va="top")
    addtext(ax4, 0.975, 0.95, r"(c)", ha="right", va="top")

    if save:
        fig.savefig(FIGS / "V_dos_conv2.png", dpi=900)


def plot_sigma_iw_v(save=False):
    print("---- Sigma(iω) ----")
    root = ROOT / "V" / "sws_opt" / "400K"
    xlim = 0, +10

    folder = EmtoDirectory(root, "u20")
    iw, sig_iw = folder.get_sigma_iw(unit="ev")
    sig_iw = sig_iw[0]

    deriv = deriv_iw(iw, sig_iw)
    meff = effective_mass(deriv)
    print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.3f}")
    print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.3f}")

    y0 = sig_iw.imag[:, :, 0] - deriv * iw[0]
    poly_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    poly_eg = Polynomial([y0[0, 2], deriv[0, 2]])

    fig, ax = plt.subplots()

    m1, m2 = MARKERS[:2]
    ax.set_xlabel(r"$i \omega_n$ (eV)")
    ax.set_ylabel(r"Im $\Sigma$ (eV)")
    ax.plot(iw, sig_iw.imag[0, 0], marker=m1, color="C0", label=r"$t_{2g}$", zorder=2)
    ax.plot(iw, sig_iw.imag[0, 2], marker=m2, color="C1", label=r"$e_{g}$", zorder=2)
    ax.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    ax.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.38, 0.01)

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.38, 0.6, 0.5, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    ax2.axhline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)
    z, sig_z = folder.get_sigma_z(unit="ev")
    ax2.plot(z, sig_z.imag[0, 0], color="C0", label=r"$t_{2g}$", zorder=2)
    ax2.plot(z, sig_z.imag[0, 2], color="C1", label=r"$e_{g}$", zorder=2)
    ax2.tick_params(labelsize=5)
    ax2.set_xlabel(r"$E - E_F$ (eV)", fontsize=5)
    ax2.set_xlim(-3, +3)
    ax2.set_ylim(-0.35, 0.03)
    # mask = np.abs(z) < 0.2
    # fit_t2g = Polynomial.fit(z[mask], sig_z.imag[0, 0, mask], 3)
    # fit_eg = Polynomial.fit(z[mask], sig_z.imag[0, 2, mask], 3)
    # ax2.plot(*fit_t2g.linspace(domain=[-0.7, 0.7]), color="C0", ls="--", lw=0.5, zorder=1)
    # ax2.plot(*fit_eg.linspace(domain=[-1.1, 1.1]), color="C1", ls="--", lw=0.5, zorder=1)
    ax.legend(loc="lower left")
    if save:
        fig.savefig(FIGS / "V_selfiw.png", dpi=900)


def plot_meff2_v(save=False):
    print("---- m*(U) ----")
    xlim = 1.4, 4.6
    ylim = 1, 2.12

    root = ROOT / "V" / "sws_opt"
    fig = plt.figure(figsize=[3.375, 0.75 * 2.531])
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.02, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_yticklabels([])
    ax1.text(0.8, 0.93, "t$_{2g}$", transform=ax1.transAxes, ha="center", va="center")
    ax2.text(0.2, 0.93, "e$_{g}$", transform=ax2.transAxes, ha="center", va="center")
    ax1.set_ylabel(r"$m^* / m$")
    ax1.set_xlabel(r"$U$ (eV)")
    ax2.set_xlabel(r"$U$ (eV)")
    ax1.grid(axis="y")
    ax2.grid(axis="y")

    temps = [200, 400, 600]
    lines = ["-", "--", ":", "-."]
    for temp, marker, ls in zip(temps, MARKERS, lines):
        label = f"$T={temp}$K"
        uu, meffs = extract_meffs_u(root / f"{temp}K")
        mask = np.logical_and(1.5 <= uu, uu <= 4.5)
        uu = uu[mask]
        meffs = meffs[mask]
        # meff_total = (3 * meffs[:, 0] + 2 * meffs[:, 1]) / 5
        ax1.plot(uu, meffs[:, 0], ls=ls, marker=marker, ms=1.5, color="C0", label=label)
        ax2.plot(uu, meffs[:, 1], ls=ls, marker=marker, ms=1.5, color="C1", label=label)

    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax2.set_ylim(*ylim)
    ax1.legend(fontsize=5)
    ax2.legend(fontsize=5)

    if save:
        fig.savefig(FIGS / "V_u_meff2.png", dpi=900)


def plot_sws_lambda_tc_v(save=False):
    print("---- SWS, λ, Tc ----")
    root = ROOT / "V" / "sws_opt" / "400K"
    xlim = -0.1, 5.1

    uu, alat, etas = list(), list(), list()
    for folder in walk_emtodirs(root):
        dat = folder.dat
        prn = folder.prn
        if prn is None or not prn.converged:
            continue
        v = dat.get_atom("V")
        u = v.get_u(2)
        if u > 5.0:
            continue
        _, _, _alat = prn.get_lattice_constants()
        hopfield = prn.get_sublat_hopfields(dat).sum(axis=1)[0]
        uu.append(u)
        alat.append(_alat)
        etas.append(hopfield)

    idx = np.argsort(uu)
    uu = np.array(uu)[idx]
    alat = np.array(alat)[idx]
    etas = np.array(etas)[idx]
    mask = np.logical_or(1.5 <= uu, uu == 0)
    uu = uu[mask]
    alat = alat[mask]
    etas = etas[mask]

    mass = Vanadium.mass
    theta = Vanadium.debye_0K
    lambd = phonon_coupling(2 * etas, mass, theta)

    # Plot
    fig = plt.figure(figsize=[3.375, 1.0 * 2.531])
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 2])
    gs.update(left=0.15, bottom=0.10, top=0.97, right=0.97, wspace=0.0, hspace=0.02)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    m1, m2 = MARKERS[:2]
    # ax1.axhline(3.024, color="k", ls="-.", lw=0.5, label="$a$ (exp)")
    ax1.plot(uu, alat, marker=m1, ms=2, lw=0.8, color="grey", )
    ax2.plot(uu, lambd, marker=m1, ms=2, lw=0.8, color="grey", )

    def label(mu):
        return r"$\mu^* = " + f"{mu:.2f}" + "$"

    mustar = 0.13
    tc = mcmillan(theta, lambd, mu_star=mustar)
    ax3.plot(uu, tc, marker=m1, ms=2, lw=0.8, color="C2", label=label(mustar))

    mustar = 0.14
    tc = mcmillan(theta, lambd, mu_star=mustar)
    ax3.plot(uu, tc, ls="--", marker=m2, ms=2, lw=0.8, color="C3", label=label(mustar))

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

    ax3.axhline(5.38, color="C7", ls="--", lw=0.5, zorder=0, label="exp")
    ax3.set_ylim(4.5, 6.3)
    ax3.legend()

    addtext(ax1, 0.975, 0.9, r"(a)", ha="right", va="top")
    addtext(ax2, 0.975, 0.9, r"(b)", ha="right", va="top")
    addtext(ax3, 0.975, 0.95, r"(c)", ha="right", va="top")

    if save:
        fig.savefig(FIGS / "V_alat_lambda_tc.png", dpi=900)


def _callback(folder, ds):
    ds.attrs["conc"] = int(folder.name.replace("Ti", "")) / 100


# noinspection PyTypeChecker,PyUnresolvedReferences,PyCallingNonCallable
def plot_conc_alat_tiv(save=False):
    fig, ax = plt.subplots(figsize=[3.375, 0.8 * 2.531])

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$a_{eq}$ ($\AA$)")
    ax.set_xlim(-0.02, 0.7)
    # ax.set_ylim(2.98, 3.18)
    ax.grid(axis="x")

    u = 2

    root = ROOT / "Ti-V" / "CPA" / f"sws_u{u}_400K"
    # save_alat_etots(root, _callback=callback)
    cc = list()
    alat = list()
    with h5py.File(root / "sws.hdf5", "r") as f:
        for ds in f.values():
            cc.append(ds.attrs["conc"])
            alat.append(ds.attrs["alat_opt"][0])
    idx = np.argsort(cc)
    cc = np.array(cc)[idx]
    alat = np.array(alat)[idx]
    mask = cc < 0.8
    linfit = Polynomial.fit(cc[mask], alat[mask], 1, domain=[0, 0.7])

    m = MARKERS[0]
    ax.plot(cc[mask], alat[mask], marker=m, lw=0, ms=3, label="data", zorder=2)
    ax.plot(*linfit.linspace(100), color="k", label="linear fit", zorder=1)
    print(linfit(1))

    ax.legend()

    if save:
        fig.savefig(FIGS / "TiV_c_alat.png", dpi=900)


def plot_dos_cpa_tiv(save=False):
    xlim = -8, +8
    root = ROOT / "Ti-V" / "CPA" / "u2_400K"

    fig, ax = plt.subplots(figsize=[3.375, 1 * 2.531])  # figsize=[3.375, 0.75 * 2.531])
    ax.grid(axis="x", zorder=-1)
    ax.set_axisbelow(True)
    ax.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    folder = EmtoDirectory(root / "Ti35")
    c = folder.dat.get_concentration("Ti")
    dosfile = folder.dos
    # energy, dos = dosfile.get_total_dos()
    ee, dos = dosfile.get_partial_dos(atom="Ti")
    ax.plot(ee * RY2EV, dos, lw=0.8, zorder=1, ls="-", label="$c_{Ti}$ Ti")
    ee, dos = dosfile.get_partial_dos(atom="V")
    ax.plot(ee * RY2EV, dos, lw=0.8, zorder=1, ls="-", label="$c_{V}$ V")
    ee, dos = dosfile.get_total_dos()
    ax.plot(ee * RY2EV, dos, lw=0.8, zorder=2, color="k", label="Total")

    ax.text(0.53, 0.95, "Ti$_{0.35}$V$_{0.65}$", ha="left", va="top", transform=ax.transAxes)

    root = CONFIG["app"] / "V" / "sws_opt" / "400K_2" / "u20"
    folder = EmtoDirectory(root)
    dosfile = folder.dos
    eev, dosv = dosfile.get_total_dos()
    # ax.plot(ee * ry2ev, dos, lw=0.7, zorder=1, ls="--", label="pure V")

    ax.set_xlabel("$E - E_F$ (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 30)
    ax.legend()

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.18, 0.55, 0.25, 0.35]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(ee * RY2EV, dos, lw=0.7, zorder=1, label="TiV", color="k")
    ax2.plot(eev * RY2EV, dosv, lw=0.7, zorder=1, ls="-", label="V", color="C1")
    ax2.legend(frameon=False, fontsize=5)
    ax2.tick_params(labelsize=5)
    ax2.set_xlim(-1.3, 1.3)
    ax2.set_ylim(0, 30)
    ax2.axvline(0, color="dimgrey", ls="-", lw=0.5, zorder=1)

    if save:
        fig.savefig(FIGS / "TiV_DOS_CPA.png", dpi=900)


def plot_sigma_iw2_tiv(save=False):
    print("---- Sigma(iω) ----")
    inset = False
    u = 2
    c = 0.35
    root = ROOT / "Ti-V" / "CPA"
    xlim = 0, +4.5
    ylim = -0.42, 0.01

    folder = EmtoDirectory(root, f"u{u}_400K", f"Ti{int(c*100)}")
    iw, sig_iw = folder.get_sigma_iw(unit="ev")
    sig_iw = sig_iw[:, 0]
    sig_ti, sig_v = sig_iw[0], sig_iw[1]

    deriv = deriv_iw(iw, sig_iw)
    y0 = sig_iw[..., 0].imag - deriv * iw[0]

    poly_ti_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    poly_ti_eg = Polynomial([y0[0, 2], deriv[0, 2]])
    poly_v_t2g = Polynomial([y0[1, 0], deriv[1, 0]])
    poly_v_eg = Polynomial([y0[1, 2], deriv[1, 2]])
    domain = [0, 5]

    # sig_iw_tot = np.average(sig_iw, weights=[c, 1 - c], axis=0)
    # deriv = deriv_iw(iw, sig_iw_tot)
    # meff = effective_mass(deriv)
    # print(f"Z_t2g = {1 / meff[0, 0]:.6f}  m* = {meff[0, 0]:.3f}")
    # print(f"Z_eg  = {1 / meff[0, 2]:.6f}  m* = {meff[0, 2]:.3f}")
    # y0 = sig_iw_tot.imag[:, :, 0] - deriv * iw[0]
    # poly_t2g = Polynomial([y0[0, 0], deriv[0, 0]])
    # poly_eg = Polynomial([y0[0, 2], deriv[0, 2]])

    fig = plt.figure()  # figsize=[3.375, 1.0 * 2.531])
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.02, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_yticklabels([])
    ax1.text(0.3, 0.95, "t$_{2g}$", transform=ax1.transAxes, ha="center", va="center")
    ax2.text(0.3, 0.95, "e$_{g}$", transform=ax2.transAxes, ha="center", va="center")
    # plot_tc_conc_tiv(save)
    ax1.set_xlabel(r"$i \omega_n$ (eV)")
    ax2.set_xlabel(r"$i \omega_n$ (eV)")
    ax1.set_ylabel(r"Im $\Sigma$ (eV)")
    m1, m2, m3 = MARKERS[:3]
    # ax1.plot(iw, sig_iw_tot.imag[0, 0], marker=m1, color="C0", label=r"Total", zorder=2)
    # ax2.plot(iw, sig_iw_tot.imag[0, 2], marker=m1, color="C1", label=r"Total", zorder=2)
    # ax1.plot(*poly_t2g.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C0", zorder=1)
    # ax2.plot(*poly_eg.linspace(domain=[0, 5]), ls="--", lw=0.7, color="C1", zorder=1)

    ax1.plot(iw, sig_ti.imag[0], marker=m1, ls="--", color="C0", label=r"Ti", zorder=2)
    ax2.plot(iw, sig_ti.imag[2], marker=m1, ls="--", color="C1", label=r"Ti", zorder=2)
    ax1.plot(*poly_ti_t2g.linspace(domain=domain), ls=":", color="C7", zorder=1)
    ax2.plot(*poly_ti_eg.linspace(domain=domain), ls=":", color="C7", zorder=1)

    ax1.plot(iw, sig_v.imag[0], marker=m2, ls="--", color="C0", label=r"V", zorder=2)
    ax2.plot(iw, sig_v.imag[2], marker=m2, ls="--", color="C1", label=r"V", zorder=2)
    ax1.plot(*poly_v_t2g.linspace(domain=domain), ls=":", color="C7", zorder=1)
    ax2.plot(*poly_v_eg.linspace(domain=domain), ls=":", color="C7", zorder=1)

    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax2.set_ylim(*ylim)
    ax1.grid(axis="y")
    ax2.grid(axis="y")
    ax1.legend()
    ax2.legend()
    if save:
        fig.savefig(FIGS / "TiV_selfiw2.png", dpi=900)


def plot_meff2_tiv(save=False):
    print("---- m*(U) ----")
    xlim = -0.03, 0.68
    ylim = 1.05, 1.38

    root = ROOT / "Ti-V" / "CPA"
    fig = plt.figure(figsize=[3.375, 1 * 2.531])
    gs = gridspec.GridSpec(2, 2)
    gs.update(left=0.15, bottom=0.13, top=0.97, right=0.97, wspace=0.02, hspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax2.set_yticklabels([])
    ax4.set_yticklabels([])
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    addtext(ax1, 0.1, 0.92, "t$_{2g}$", "left", "top")
    addtext(ax2, 0.1, 0.92, "e$_{g}$", "left", "top")

    # ax1.text(0.1, 0.92, "t$_{2g}$", transform=ax1.transAxes, ha="center", va="center")
    # ax2.text(0.1, 0.92, "e$_{g}$", transform=ax2.transAxes, ha="center", va="center")
    ax1.set_ylabel(r"Ti $m^* / m$")
    ax3.set_ylabel(r"V $m^* / m$")
    ax3.set_xlabel(r"$x$")
    ax4.set_xlabel(r"$x$")
    u = 2
    temp = 200
    m1, m2 = MARKERS[:2]
    cc, meffs = extract_meffs_c(root / f"u{u}_{temp}K")
    print(meffs.shape)
    ax1.plot(cc, meffs[:, 0, 0], marker=m1, color="C0", label=f"$T={temp}$K")
    ax2.plot(cc, meffs[:, 1, 0], marker=m1, color="C1", label=f"$T={temp}$K")
    ax3.plot(cc, meffs[:, 0, 1], marker=m1, color="C0", label=f"$T={temp}$K")
    ax4.plot(cc, meffs[:, 1, 1], marker=m1, color="C1", label=f"$T={temp}$K")
    # ax2.plot(cc, meffs[:, 1], marker=m1, color="C1", label=f"$T={temp}$K")

    temp = 400
    cc, meffs = extract_meffs_c(root / f"u{u}_{temp}K")
    ax1.plot(cc, meffs[:, 0, 0], marker=m2, color="C0", label=f"$T={temp}$K")
    ax2.plot(cc, meffs[:, 1, 0], marker=m2, color="C1", label=f"$T={temp}$K")
    ax3.plot(cc, meffs[:, 0, 1], marker=m2, color="C0", label=f"$T={temp}$K")
    ax4.plot(cc, meffs[:, 1, 1], marker=m2, color="C1", label=f"$T={temp}$K")
    # ax1.plot(cc, meffs[:, 0, 0], marker=m1, color="C0", label=f"$T={temp}$K")
    # ax2.plot(cc, meffs[:, 0, 1], marker=m1, color="C0", label=f"$T={temp}$K")
    # ax1.plot(cc, meffs[:, 0], marker=m2, ls="--", color="C0", label=f"$T={temp}$K")
    # ax2.plot(cc, meffs[:, 1], marker=m2, ls="--", color="C1", label=f"$T={temp}$K")

    ax1.set_xlim(*xlim)
    ax2.set_xlim(*xlim)
    ax3.set_xlim(*xlim)
    ax4.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax2.set_ylim(*ylim)
    ax3.set_ylim(*ylim)
    ax4.set_ylim(*ylim)
    ax1.grid(axis="y")
    ax2.grid(axis="y")
    ax3.grid(axis="y")
    ax4.grid(axis="y")
    # ax1.legend()
    ax1.legend()
    ax2.legend()
    if save:
        fig.savefig(FIGS / "TiV_c_meff2.png", dpi=900)


def read_tc_data():
    temp = 400
    ckey = "Ti"
    root = ROOT / "Ti-V" / "CPA"
    uu = list()
    data = list()
    for path in root.iterdir():
        name = path.name
        if name.startswith("sws"):
            continue
        if not name.endswith(f"{temp}K"):
            continue
        u = float(name.split("_")[0][1:])
        if u >= 6:
            continue
        uu.append(u)
        cc, tc = extract_tc_c_tiv(path, ckey)
        data.append((cc, tc))
    i = np.argsort(uu)
    uu = np.array(uu)[i]
    data = [np.array(data[j]) for j in i]
    return uu, data


def plot_tc_conc_tiv(save=False):
    fig, ax = plt.subplots()
    ckey = "Ti"
    uu = [0, 2, 5]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5"]
    for u, c, m in zip(uu, colors, MARKERS):
        lab = f"$U={u}$eV" if u > 0 else "LDA"
        root = ROOT / "Ti-V" / "CPA" / f"u{u}_400K"  # / f"nl3_u{u}"
        cc, tc = extract_tc_c_tiv(root, ckey)
        ax.plot(cc, tc, f"--{m}", lw=0.5, ms=1.5, color=c, label=lab)

    regions = [1 - 0.335, 1 - 0.145, 1]
    addtext(ax, 0.1, 0.9, r"Ti$_{x}$V$_{1-x}$")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(0, 9)
    ax.set_xlabel(r"x")
    ax.set_ylabel("$T_c$ (K)")
    ax.fill_between([regions[-2] - 0.02, 1.2], 0, 9, color="k", alpha=0.1, lw=0)
    ax.fill_between([regions[-3], regions[-2] + 0.02], 0, 9, color="red", alpha=0.1, lw=0)
    x1 = regions[-1] + 0.5 * (regions[-2] - regions[-1])
    x2 = regions[-2] + 0.5 * (regions[-3] - regions[-2])
    y = 0.7
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
    ax.plot(x, y, "-", marker=".", lw=0.5, ms=2.5, color="C7", label="exp", zorder=1)
    ax.legend(loc="upper right")

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.2, 0.3, 0.4, 0.35]
    ax2 = fig.add_axes([left, bottom, width, height])
    cc = [0.1, 0.3, 0.6]
    uu, data = read_tc_data()
    for c, m in zip(cc, MARKERS):
        tc = np.zeros(len(uu), dtype=float)
        for i in range(len(uu)):
            cc, ttc = data[i]
            idx = np.where(cc == c)[0][0]
            tc[i] = ttc[idx]
        ax2.plot(uu, tc, marker=m, ms=1, lw=0.5, label=f"$x={c:.1f}$")
    ax2.tick_params(labelsize=5)
    ax2.set_xlabel(r"$U$", fontsize=5)
    ax2.legend(loc="center left", frameon=False, fontsize=5)
    ax2.set_ylim(6.1, None)
    if save:
        fig.savefig(FIGS / "TiV_conc_tc.png", dpi=900)


def main():
    save = True
    # Vanadium
    # plot_alat_opt_curves_v(save)
    # plot_dos_conv_v(save)
    # plot_dos_conv_v2(save)
    # plot_sigma_iw_v(save)
    # plot_meff2_v(save)
    # plot_sws_lambda_tc_v(save)
    # TiV
    # plot_conc_alat_tiv(save)
    # plot_dos_cpa_tiv(save)
    # plot_sigma_iw2_tiv(save)
    plot_meff2_tiv(save)
    # plot_tc_conc_tiv(save)
    plt.show()


if __name__ == "__main__":
    main()
