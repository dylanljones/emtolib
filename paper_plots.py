# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-08-21

import numpy as np
import matplotlib.pyplot as plt
from mplstyles import use_mplstyle
from emtolib import CONFIG, Path, walk_emtodirs, EmtoDirectory, elements
from emtolib.mcmillan import phonon_coupling, mcmillan

Vanadium = elements["V"]


def plot_tc_u_vanadium():

    use_mplstyle("figure", "aps")
    root = CONFIG["app"] / "V" / "nl3_j06"
    folders = list(walk_emtodirs(root))
    uu, etas = [], []
    for folder in folders:
        dat = folder.dat
        prn = folder.prn
        u = dat.get_atom("V").u[2]
        if u <= 5:
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

    fig, ax = plt.subplots()
    ax.plot(uu, tc, "-o")
    ax.set_xlabel("U (eV)")
    ax.set_ylabel("Tc (K)")
    ax.set_ylim(3.15, 4.05)
    ax.set_xlim(0.9, 5.1)
    ax.grid()
    fig.savefig("vanadium_tc_u.png")
    plt.show()


def plot_dos_vanadium():
    pass



def main():
    plot_tc_u_vanadium()



if __name__ == "__main__":
    main()
