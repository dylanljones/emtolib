# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2024-11-07

"""Methods for calculating McMillan's critical temperature.

Methods for calculating the critical temperature of a superconductor
using McMillan's equation via the electron-phonon coupling constant computed
from the Hopfield parameters output by the EMTO code.
"""

import numpy as np
from scipy import constants as const


def deriv_iw(iw, sig_iw):
    """Calculate the derivative from the imaginary part of the self-energy."""
    return (sig_iw[..., 1] - sig_iw[..., 0]).imag / (iw[1] - iw[0])


def effective_mass(deriv):
    """Calculate the effective mass from the derivative of the self-energy."""
    return 1 - deriv


def effective_mass_iw(iw, sig_iw):
    """Calculate the effective mass from the imaginary part of the self-energy."""
    return effective_mass(deriv_iw(iw, sig_iw))


def quasiparticle_weight(deriv):
    """Calculate the quasiparticle weight from the derivative of the self-energy."""
    return 1 / (1 - deriv)


def quasiparticle_weight_iw(iw, sig_iw):
    """Calculate the quasiparticle weight from the derivative of the self-energy."""
    return quasiparticle_weight(deriv_iw(iw, sig_iw))


def quasiparticle_damping_iw(iw, sig_iw):
    """Calculate the quasiparticle damping from the Matsubara self-energy."""
    return -quasiparticle_weight_iw(iw, sig_iw) * sig_iw[..., 0].imag


def phonon_coupling(eta, mass, theta_d):
    """Compute the phonon coupling constant .math:`λ` used in McMillan's equation.

    Parameters
    ----------
    eta : float or np.ndarray
        The Hopfield of the atom (in .math:`eV/Å²`).
    mass : float or np.ndarray
         The atomic mass of the atom (in .math:`u`).
    theta_d : float or np.ndarray
         The Debye temperature at 0K of the atom (in .math:`K`).

    Returns
    -------
    lamb : float or np.ndarray
        The dimensionless phonon coupling constant.
    """
    kelvin_to_ev = const.k / const.eV

    factor = const.eV / (const.angstrom**2 * const.m_u)
    x = factor * eta / mass  # SI Unit: 1/s²
    x /= (2 * np.pi * const.c) ** 2  # SI Unit: 1/m²
    theta = kelvin_to_ev * theta_d  # SI Unit: eV
    theta = theta * const.eV / (const.h * const.c)  # SI Unit: 1/m
    lamb = x / (theta**2)  # Dimension-less
    return lamb


def mcmillan(theta_d, lamb, mu_star):
    r"""Compute the critical temperature using the McMillan equation

    .. math::
        T_C = \frac{θ_D}{1.45} \exp[ -\frac{1.04(1 + λ)}{λ - μ⁺(1 + 0.62λ)} ]

    Parameters
    ----------
    theta_d : array_like
        The Debye temperature .math:`θ_D`.
    lamb : array_like
         The electron—phonon coupling term .math:`λ`.
    mu_star : float or np.ndarray
        The Coulomb pseudo-potential .math:`μ⁺`.

    Returns
    -------
    tc : float or np.ndarray
        The critical temperature .math:`T_C` of the superconductor.
    """
    arg = -1.04 * (1.0 + lamb) / (lamb - mu_star * (1.0 + 0.62 * lamb))
    return theta_d / 1.45 * np.exp(arg)


def compute_mcmillan_tc(concs, hopfields, masses, thetas, mu_star):
    r"""Compute the McMillan critical temperature for disordered superconductors.

    Parameters
    ----------
    concs : (..., N) np.ndarray
        The concentrations for each of the `N` atom types.
    hopfields : (..., N) np.ndarray
        The hopfields for each of the `N` atom types in eV/AA^2.
    masses : (N,) np.ndarray
        The masses for each of the `N` atom types in Kg.
    thetas : (N,) np.ndarray
        The Debye temperatures of the `N` atom types in Kelvin.
    mu_star : float
        The pseudo-potential .math:`\mu^*`.

    Returns
    -------
    tc : (...) np.ndarray
        The critical temperature for each concentration.
    """
    masses_avg = np.sum(concs * masses, axis=1)
    thetas_avg = np.sum(concs * thetas, axis=1)
    hopfields_avg = np.sum(concs * hopfields, axis=-1)

    lambdas = phonon_coupling(hopfields_avg, masses_avg, thetas_avg)
    return mcmillan(thetas_avg, lambdas, mu_star)
