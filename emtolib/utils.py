# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-10-11

import numpy as np
from scipy.special import expit


def fermi_fct(eps, beta):
    r"""
    Return the Fermi function `1/(exp(βϵ)+1)`.

    For complex inputs the function is not as accurate as for real inputs.

    Parameters
    ----------
    eps : complex or float or np.ndarray
        The energy at which the Fermi function is evaluated.
    beta : float
        The inverse temperature :math:`beta = 1/k_B T`.

    Returns
    -------
    complex of float or np.ndarray
        The Fermi function, same type as eps.

    See Also
    --------
    fermi_fct_inv : The inverse of the Fermi function for real arguments.
    """
    z = eps * beta
    try:
        return expit(-z)  # = 0.5 * (1. + tanh(-0.5 * beta * eps))
    except TypeError:
        pass  # complex arguments not handled by expit
    z = np.asanyarray(z)
    pos = z.real > 0
    res = np.empty_like(z)
    res[~pos] = 1.0 / (np.exp(z[~pos]) + 1)
    exp_m = np.exp(-z[pos])
    res[pos] = exp_m / (1 + exp_m)
    return res


def lorentzian(x, x0, width):
    return width / (np.pi * (width**2 + (x - x0) ** 2))


def gaussian(x, x0, width):
    return np.exp(-((x - x0) ** 2) / (2 * width**2)) / (width * np.sqrt(2 * np.pi))


def convolve_func(x, y, func, normalize=True, **kwargs):
    kernel = func(x, **kwargs)
    yconv = np.convolve(y, kernel, mode="same")
    if normalize:
        yconv /= np.sum(kernel)
    return yconv
