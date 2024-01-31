# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-03

from pathlib import Path  # noqa: F401

from .files import (
    Atom,
    DosFile,
    PrnFile,
    BmdlFile,
    KgrnFile,
    Makefile,
    SlurmScript,
    load_dos,
    read_dos,
    generate_makefile,
)
from .common import logger, elements
from .config import CONFIG, read_config, update_emto_paths, update_slurm_settings
from .errors import (
    DOSError,
    DMFTError,
    EmtoError,
    KGRNError,
    DOSReadError,
    DMFTReadError,
    KGRNReadError,
    DMFTWriteError,
    KGRNWriteError,
)
from .directory import EmtoDirectory, is_emtodir, walk_emtodirs

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
