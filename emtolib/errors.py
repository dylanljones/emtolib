# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-26

"""Defines all exceptions used in emtolib."""


class EmtoError(Exception):
    """Base class for all exceptions in emtolib."""

    pass


class KGRNError(EmtoError):
    """Base class for all errors related to KGRN files."""

    pass


class KGRNReadError(KGRNError):
    """Errors related to reading KGRN files."""

    pass


class KGRNWriteError(KGRNError):
    """Errors related to writing KGRN files."""

    pass


class DMFTError(EmtoError):
    """Base class for all errors related to DMFT files."""

    pass


class DMFTReadError(DMFTError):
    """Errors related to reading DMFT files."""

    pass


class DMFTWriteError(DMFTError):
    """Errors related to writing DMFT files."""

    pass


class DOSError(EmtoError):
    """Base class for all errors related to DOS files."""

    pass


class DOSReadError(DOSError):
    """Errors related to reading DOS files."""

    pass
