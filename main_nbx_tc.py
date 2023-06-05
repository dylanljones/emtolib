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
from emtolib.mcmillan import phonon_coupling, mcmillan
from emtolib.xarr import load_dataset, update_datasets
from scipy import interpolate, optimize

FIGDIR = Path("figs")


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


def average(x, ds, axis=-1):
    return np.sum(x * ds.coords["concs"].values, axis=axis)


def read_data(datasets):
    concs = list(sorted(datasets.keys()))
    hopfields = np.zeros((len(concs), 2))
    masses = np.zeros((len(concs)))
    thetas = np.zeros((len(concs)))
    for i, k in enumerate(concs):
        ds = datasets[k]
        mass = ds.attrs["mass"]
        theta = ds.attrs["debye_temp"]
        eta = ds.hopfields.values
        if i == 0:
            eta = [0, eta[0]]
        elif i == len(concs) - 1:
            eta = [eta[0], 0]
        hopfields[i] = eta
        masses[i] = mass
        thetas[i] = theta

    concs = np.array(concs)
    concs = np.array([concs, 1-concs]).T
    return concs, hopfields, masses, thetas


def fit_mustar(lamb, theta_d, tc):

    def cost_funct(mu_star):
        return np.sum((tc - mcmillan(theta_d, lamb, mu_star))**2)

    sol = optimize.minimize(cost_funct, x0=0.13 * np.ones_like(tc))
    return sol


def plot_eta_lambda(x, hopfields_avg, lambdas, atoms, at_hopfields, at_lambdas):
    fig = plt.figure()
    gs = fig.add_gridspec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    gs.update(wspace=0.025, hspace=0.05, bottom=0.15, left=0.15, right=0.95, top=0.95)
    ax1.xaxis.set_tick_params(labelbottom=False)

    # fig, (ax1, ax2) = plt.subplots(nrows=2, sharex="all")
    # fig.subplots_adjust(wspace=0.05, hspace=0.05)
    ax1.axhline(at_hopfields[0], label=atoms[0], color="C0", ls="--", lw=0.8)
    ax1.axhline(at_hopfields[1], label=atoms[1], color="C1", ls="--", lw=0.8)
    # ax1.plot(x, hopfields[:, 0], "o-", label=r"Nb")
    # ax1.plot(x, hopfields[:, 1], "o-", label=r"Ta")
    ax1.plot(x, hopfields_avg, "o-", label=r"$\langle \eta \rangle$", color="k", ms=2)
    ax1.set_xlabel("x")
    # ax.set_xlabel("$c_{Nb}$")
    ax1.set_ylabel(r"$\eta$ (eV/ $\mathrm{\AA}^2$)")
    ax1.legend()
    ax1.set_xmargin(0.02)

    # ax2.plot(x, lambda_nb, "o-", label=r"Nb")
    # ax2.plot(x, lambda_ta, "o-", label=r"Ta")
    ax2.axhline(at_lambdas[0], label=atoms[0], color="C0", ls="--", lw=0.8)
    ax2.axhline(at_lambdas[1], label=atoms[1], color="C1", ls="--", lw=0.8)
    ax2.plot(x, lambdas, "o-", label=r"$\langle \lambda \rangle$", color="k", ms=2)
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$\lambda$ (eV)")
    ax2.legend()
    return fig, (ax1, ax2)


def plot_tc(x, tc, x_exp, tc_exp):
    fig, ax = plt.subplots()
    ax.plot(x, tc, "o-", label="McMillan", ms=2)
    # ax.plot(x, tc_var, "o-", label="McMillan (var)", ms=2)
    ax.plot(x_exp, tc_exp, "o-", color="k", label="Exp", ms=2)
    ax.set_xlabel("x")
    # ax.set_xlabel("$c_{Nb}$")
    ax.set_ylabel(r"$T_C$ (K)")
    ax.set_ylim(0, None)
    ax.legend()
    ax.set_xmargin(0.02)
    return fig, ax


def main():
    mu_star_ta = 0.1027
    mu_star_v = 0.1763
    mu_star_nb = 0.1953

    at1, at2 = "Nb", "V"
    root = Path("app") / f"{at1}-{at2}" / "CPA2"
    update_datasets(root, force=False)

    atoms = (at1, at2)
    datasets = load_datasets(root)
    for k, ds in datasets.items():
        if k not in (0, 1):
            assert tuple(ds.attrs["atoms"]) == atoms
        else:
            assert ds.attrs["atoms"][0] in atoms

    exp = Path("exp") / f"{at1.lower()}_{at2.lower()}_exp.dat"
    data = np.loadtxt(exp)
    x_exp = 1 - data[:, 0]
    tc_exp = data[:, 1]

    concs, hopfields, masses, thetas = read_data(datasets)
    x = concs[:, 0]
    hopfields_avg = np.sum(concs * hopfields, axis=-1)
    lambdas = phonon_coupling(hopfields_avg, masses, thetas)

    # Comput McMillan Tc with fixed mu_star
    mu_star = x * mu_star_nb + (1 - x) * mu_star_ta
    tc = mcmillan(thetas, lambdas, mu_star=0.13)
    # print(f"T_C  ", "  ".join(f"{t:.2f}" for t in tc))

    # Compute pure quantities
    # hopfields[0, 0] = np.nan
    # hopfields[-1, 1] = np.nan
    at_masses = datasets[0.5].coords["masses"].values
    at_thetas = datasets[0.5].coords["debye"].values
    at_hopfields = [hopfields[-1, 0], hopfields[0, 1]]
    at_lambdas = [
        phonon_coupling(at_hopfields[0], at_masses[0], at_thetas[0]),
        phonon_coupling(at_hopfields[1], at_masses[1], at_thetas[1])
    ]

    # Fit mu_star to (linear interpolated) experimental data
    tc_inter = interpolate.interp1d(x_exp, tc_exp, kind="linear")(x)
    sol = fit_mustar(lambdas, thetas, tc_inter)
    mu_star = sol.x
    print(mu_star[0], mu_star[-1])
    tc_var = mcmillan(thetas, lambdas, mu_star=mu_star_nb)

    fig, ax = plt.subplots()
    ax.plot(x, mu_star, "o-")
    ax.axhline(0.13, lw=.5)
    # ax.set_ylim(0.1, 0.2)

    use_mplstyle("figure", "aps", color_cycle="seaborn-colorblind")

    fig, axs = plot_eta_lambda(x, hopfields_avg, lambdas, atoms, at_hopfields, at_lambdas)
    fig.savefig(FIGDIR / f"{at1}{at2}_eta_lambda_cpa.png")

    fig, ax = plot_tc(x, tc, x_exp, tc_exp)
    fig.savefig(FIGDIR / f"{at1}{at2}_tc_cpa.png")
    plt.show()


if __name__ == "__main__":
    main()
