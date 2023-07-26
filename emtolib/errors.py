# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-26

"""This module defines all exceptions used in emtolib."""


class EmtoException(Exception):
    """Base class for all exceptions in emtolib."""
    pass


# Input and Output file errors

class KGRNError(EmtoException):
    """Base class for all errors related to KGRN files."""
    pass


class KGRNReadError(KGRNError):
    pass


class KGRNWriteError(KGRNError):
    pass


class DMFTError(EmtoException):
    """Base class for all errors related to DMFT files."""
    pass


class BMDLError(EmtoException):
    """Base class for all errors related to BMDL files."""
    pass


class PRNError(EmtoException):
    """Base class for all errors related to PRN files."""
    pass


class DOSError(EmtoException):
    """Base class for all errors related to DOS files."""
    pass


class SlurmError(EmtoException):
    """Base class for all errors related to SLURM files."""
    pass


class MakefileError(EmtoException):
    """Base class for all errors related to Makefiles."""
    pass


# Directory errors

class EmtoDirectoryError(EmtoException):
    """Base class for all errors related to directories."""
    pass
