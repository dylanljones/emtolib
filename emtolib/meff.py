# -*- coding: utf-8 -*-
# Author: Dylan Jones
# Date:   2023-10-09


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
