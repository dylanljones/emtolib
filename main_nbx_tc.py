# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones


from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mplstyles import use_mplstyle
from emtolib import CONFIG
from emtolib.mcmillan import phonon_coupling, mcmillan
from emtolib.xarr import load_datasets, update_datasets
from scipy import optimize

FIGDIR = Path(CONFIG["fig_dir"])
EXPDIR = Path(CONFIG["exp_dir"])
DATADIR = Path(CONFIG["data_dir"])
XARRDIR = Path(CONFIG["xarr_dir"])


def plot_eta_lambda(x, hopfields_avg, lambdas, grid=False):
    fig = plt.figure()
    gs = fig.add_gridspec(2, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    gs.update(wspace=0.025, hspace=0.05, bottom=0.15, left=0.15, right=0.95, top=0.95)
    ax1.xaxis.set_tick_params(labelbottom=False)
    ax1.plot(x, hopfields_avg, "o-", label=r"$\langle \eta \rangle$", color="C0", ms=2)
    ax1.set_xlabel("x")
    ax1.set_ylabel(r"$\eta$ (eV/ $\mathrm{\AA}^2$)")
    ax1.legend()
    ax1.set_xmargin(0.02)

    ax2.plot(x, lambdas, "o-", label=r"$\langle \lambda \rangle$", color="C0", ms=2)
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$\lambda$ (eV)")
    ax2.legend()

    if grid:
        ax1.grid()
        ax2.grid()
    return fig, (ax1, ax2)


def plot_tc(x, tc, x_exp, tc_exp, grid=False):
    fig, ax = plt.subplots()
    ax.plot(x, tc, "o-", label="McMillan", ms=2)
    ax.plot(x_exp, tc_exp, "o-", color="k", label="Exp", ms=2)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$T_C$ (K)")
    ax.set_ylim(0, None)
    ax.legend()
    ax.set_xmargin(0.02)
    if grid:
        ax.grid()
    return fig, ax


def keyfunc(ds):
    try:
        i = ds.attrs["atoms"].index("Nb")
        return ds.coords["concs"].values[i]
    except ValueError:
        return 0


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


def main():
    at1, at2 = "Nb", "V"
    dirname = "CPA_4_275"
    mu_star = 0.13
    c_eta = 1
    savefig = False

    # Update and read datasets
    root = DATADIR / f"{at1}-{at2}" / dirname
    xarr_root = XARRDIR / f"{at1}-{at2}" / dirname
    update_datasets(root, location=xarr_root, force=False)
    datasets = load_datasets(xarr_root, keyfunc)

    # Read experimental data
    exp = EXPDIR / f"{at1.lower()}_{at2.lower()}_exp.dat"
    data = np.loadtxt(exp)
    x_exp = 1 - data[:, 0]
    tc_exp = data[:, 1]

    # Compute McMillan Tc with fixed mu_star
    concs, hopfields, masses, thetas = read_data(datasets)
    x = concs[:, 0]
    hopfields *= c_eta
    hopfields_avg = np.sum(concs * hopfields, axis=-1)
    lambdas = phonon_coupling(hopfields_avg, masses, thetas)
    tc = mcmillan(thetas, lambdas, mu_star=mu_star)

    # ===================== #
    #         Plots         #
    # ===================== #
    title = at1 + r"$_{x}$" + at2 + r"$_{1-x}$"
    use_mplstyle("figure", "aps", color_cycle="seaborn-colorblind")
    fig1, axs = plot_eta_lambda(x, hopfields_avg, lambdas, grid=True)
    axs[0].set_title(f"{title} {dirname} $C_\eta$={c_eta}")
    fig2, ax = plot_tc(x, tc, x_exp, tc_exp, grid=True)
    ax.set_title(f"{title} {dirname} $C_\eta$={c_eta}")
    if savefig:
        fig1.savefig(FIGDIR / f"{at1}{at2}_{dirname}_eta_lambda.png")
        fig2.savefig(FIGDIR / f"{at1}{at2}_{dirname}_tc.png")
    plt.show()


if __name__ == "__main__":
    main()
